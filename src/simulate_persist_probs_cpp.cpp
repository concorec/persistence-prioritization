

// src/simulate_persist_probs_cpp.cpp
//
// CRN-first persistence simulation kernels.
//
// Design contract (CRNs as default):
//   - A CRN context defines the experiment: {seed, n_draws, reps, years, chunk_size, seed_stride}.
//   - Reuse the same context across:
//       (Level 1) scanning K for a persistence curve,
//       (Level 2) paired sensitivity changes (e.g., ext_thr / rm_fraction),
//     so differences are driven by parameters, not random noise.
//   - The context does NOT store normals (no giant Z). Instead, we deterministically
//     regenerate the same Normal(0,1) stream per chunk on every call.
//
// Exports:
//   1) crn_context_create()
//   2) simulate_persist_probs_cpp()   // full p_hat vector
//   3) simulate_persist_qvec_cpp()    // q50/q16/q025 only
//   4) [optional] simulate_persist_qvec_manyK_cpp()  // vectorized K in one call
//
// -----------------------------------------------------------------------------

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

//==============================================================================
// 0) Tiny utilities (no heavy validation; just safety / determinism)
//==============================================================================

// splitmix64: seed mixer (stable, fast)
static inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// R type-7 quantile, computed in-place via nth_element
static inline double quantile_type7_inplace(std::vector<double>& xs, double p) {
  const std::size_t n = xs.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return xs[0];
  
  const double h  = 1.0 + (static_cast<double>(n) - 1.0) * p;
  const double hf = std::floor(h);
  const double hc = std::ceil(h);
  
  const std::size_t lo = static_cast<std::size_t>(hf - 1.0);
  const std::size_t hi = static_cast<std::size_t>(hc - 1.0);
  
  std::nth_element(xs.begin(), xs.begin() + lo, xs.end());
  const double xlo = xs[lo];
  if (hi == lo) return xlo;
  
  std::nth_element(xs.begin(), xs.begin() + hi, xs.end());
  const double xhi = xs[hi];
  
  const double frac = h - hf;
  return xlo + frac * (xhi - xlo);
}

//==============================================================================
// 1) CRNContext (seed-only; no stored normals)
//==============================================================================
//
// Chunking defines the deterministic RNG “substreams”:
//   chunk c uses rng seeded from base_seed64 + (c+1)*seed_stride
// and is restarted every call so the same (z) stream repeats across K or scenario.
//
struct CRNContext {
  // experiment definition
  std::uint64_t base_seed64;  // mixed seed for robustness
  int seed_stride;
  
  int n_draws;
  int reps;
  int years;
  int chunk_size;
  
  // derived layout
  int n_chunks;
  std::vector<int> start_draw;  // start index in [0, n_draws)
  std::vector<int> d_draws;     // draws in this chunk
  
  CRNContext(int seed,
             int n_draws_,
             int reps_,
             int years_,
             int chunk_size_,
             int seed_stride_)
    : base_seed64(splitmix64(static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed)))),
      seed_stride(seed_stride_),
      n_draws(n_draws_),
      reps(reps_),
      years(years_),
      chunk_size(chunk_size_) {
    
    if (n_draws <= 0 || reps <= 0 || years <= 0 || chunk_size <= 0) {
      throw std::runtime_error("CRNContext: n_draws/reps/years/chunk_size must be positive");
    }
    if (seed_stride <= 0) {
      throw std::runtime_error("CRNContext: seed_stride must be positive");
    }
    
    n_chunks = (n_draws + chunk_size - 1) / chunk_size;
    
    start_draw.resize(n_chunks);
    d_draws.resize(n_chunks);
    
    for (int c = 0; c < n_chunks; ++c) {
      const int s = c * chunk_size;
      const int d = std::min(chunk_size, n_draws - s);
      start_draw[c] = s;
      d_draws[c]    = d;
    }
  }
  
  inline std::uint64_t chunk_seed64(int c) const {
    // stable, collision-resistant mapping: (base_seed64, c) -> seed
    const std::uint64_t x =
      base_seed64 + static_cast<std::uint64_t>(c + 1) * static_cast<std::uint64_t>(seed_stride);
    return splitmix64(x);
  }
};

static inline Rcpp::XPtr<CRNContext> get_ctx(SEXP crn_ctx) {
  if (crn_ctx == R_NilValue) Rcpp::stop("crn_ctx is NULL");
  Rcpp::XPtr<CRNContext> ctx(crn_ctx);
  if (ctx.get() == nullptr) Rcpp::stop("crn_ctx pointer is invalid");
  return ctx;
}

//==============================================================================
// 2) Export: create context
//==============================================================================

// [[Rcpp::export]]
SEXP crn_context_create(const int seed,
                        const int n_draws,
                        const int reps,
                        const int years,
                        const int chunk_size,
                        const int seed_stride = 10007) {
  CRNContext* raw = nullptr;
  try {
    raw = new CRNContext(seed, n_draws, reps, years, chunk_size, seed_stride);
  } catch (const std::exception& e) {
    delete raw;
    Rcpp::stop("crn_context_create failed: %s", e.what());
  }
  return Rcpp::XPtr<CRNContext>(raw, true);
}

//==============================================================================
// 3) Internal simulation core (single source of truth)
//==============================================================================
//
// Writes p_hat (length n_draws) into out[].
//
// CRN behavior:
//   - RNG is restarted per chunk using ctx->chunk_seed64(c)
//   - For extinct trajectories we STILL draw and discard z to preserve stream alignment
//   - Therefore: repeated calls with the same ctx/r/sigma produce paired outcomes
//     across K or scenario changes.
//
static inline void simulate_core_crn(
    const double* rptr,
    const double* sptr,
    const CRNContext* ctx,
    const double K,
    const int ext_thr,
    const double cap_factor,
    double* out) {
  
  const int reps    = ctx->reps;
  const int years   = ctx->years;
  
  const double invK = 1.0 / K;
  const double capN = cap_factor * K;
  
#ifdef _OPENMP
#pragma omp parallel
{
  std::vector<double> N;  // reused per thread
  
#pragma omp for schedule(static)
  for (int c = 0; c < ctx->n_chunks; ++c) {
    const int s  = ctx->start_draw[c];
    const int d  = ctx->d_draws[c];
    const int ns = d * reps;
    
    // deterministic RNG for this chunk (restart every call)
    std::mt19937_64 rng(ctx->chunk_seed64(c));
    std::normal_distribution<double> norm01(0.0, 1.0);
    
    N.assign(static_cast<std::size_t>(ns), K);
    
    for (int t = 0; t < years; ++t) {
      for (int j = 0; j < d; ++j) {
        const double rj = rptr[s + j];
        const double sj = sptr[s + j];
        const int base_idx = j * reps;
        
        for (int k = 0; k < reps; ++k) {
          const int idx = base_idx + k;
          
          // Always advance RNG to preserve stream alignment
          const double z = norm01(rng);
          
          if (N[idx] == 0.0) continue;
          
          const double G = rj * (1.0 - N[idx] * invK) + sj * z;
          double nextN   = N[idx] * std::exp(G);
          
          if (nextN < ext_thr) nextN = 0.0;
          if (nextN > capN)    nextN = capN;
          
          N[idx] = nextN;
        }
      }
    }
    
    for (int j = 0; j < d; ++j) {
      int alive = 0;
      const int base_idx = j * reps;
      for (int k = 0; k < reps; ++k) {
        alive += (N[base_idx + k] > 0.0);
      }
      out[s + j] = static_cast<double>(alive) / static_cast<double>(reps);
    }
  }
}
#else
for (int c = 0; c < ctx->n_chunks; ++c) {
  const int s  = ctx->start_draw[c];
  const int d  = ctx->d_draws[c];
  const int ns = d * reps;
  
  std::mt19937_64 rng(ctx->chunk_seed64(c));
  std::normal_distribution<double> norm01(0.0, 1.0);
  
  std::vector<double> N(static_cast<std::size_t>(ns), K);
  
  for (int t = 0; t < years; ++t) {
    for (int j = 0; j < d; ++j) {
      const double rj = rptr[s + j];
      const double sj = sptr[s + j];
      const int base_idx = j * reps;
      
      for (int k = 0; k < reps; ++k) {
        const int idx = base_idx + k;
        
        const double z = norm01(rng);
        if (N[idx] == 0.0) continue;
        
        const double G = rj * (1.0 - N[idx] * invK) + sj * z;
        double nextN   = N[idx] * std::exp(G);
        
        if (nextN < ext_thr) nextN = 0.0;
        if (nextN > capN)    nextN = capN;
        
        N[idx] = nextN;
      }
    }
  }
  
  for (int j = 0; j < d; ++j) {
    int alive = 0;
    const int base_idx = j * reps;
    for (int k = 0; k < reps; ++k) {
      alive += (N[base_idx + k] > 0.0);
    }
    out[s + j] = static_cast<double>(alive) / static_cast<double>(reps);
  }
}
#endif
}

//==============================================================================
// 4) Export: full p_hat (CRN default)
//==============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector simulate_persist_probs_cpp(
    const Rcpp::NumericVector& r,
    const Rcpp::NumericVector& sigma,
    const double K,
    const int ext_thr,
    const double cap_factor,
    SEXP crn_ctx) {
  
  const int n_draws = r.size();
  if (sigma.size() != n_draws) Rcpp::stop("r and sigma length mismatch");
  if (!(K > 0.0))              Rcpp::stop("K must be > 0");
  
  Rcpp::XPtr<CRNContext> ctx = get_ctx(crn_ctx);
  if (ctx->n_draws != n_draws) Rcpp::stop("crn_ctx n_draws mismatch");
  
  Rcpp::NumericVector p_hat(n_draws);
  
  simulate_core_crn(
    r.begin(), sigma.begin(),
    ctx.get(),
    K, ext_thr, cap_factor,
    p_hat.begin()
  );
  
  return p_hat;
}

//==============================================================================
// 5) Export: quantiles only (CRN default)
//==============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector simulate_persist_qvec_cpp(
    const Rcpp::NumericVector& r,
    const Rcpp::NumericVector& sigma,
    const double K,
    const int ext_thr,
    const double cap_factor,
    SEXP crn_ctx) {
  
  const int n_draws = r.size();
  if (sigma.size() != n_draws) Rcpp::stop("r and sigma length mismatch");
  if (!(K > 0.0))              Rcpp::stop("K must be > 0");
  
  Rcpp::XPtr<CRNContext> ctx = get_ctx(crn_ctx);
  if (ctx->n_draws != n_draws) Rcpp::stop("crn_ctx n_draws mismatch");
  
  std::vector<double> xs(static_cast<std::size_t>(n_draws));
  
  simulate_core_crn(
    r.begin(), sigma.begin(),
    ctx.get(),
    K, ext_thr, cap_factor,
    xs.data()
  );
  
  const double q50  = quantile_type7_inplace(xs, 0.50);
  const double q16  = quantile_type7_inplace(xs, 0.16);
  const double q025 = quantile_type7_inplace(xs, 0.025);
  
  return Rcpp::NumericVector::create(
    Rcpp::Named("q50")  = q50,
    Rcpp::Named("q16")  = q16,
    Rcpp::Named("q025") = q025
  );
}