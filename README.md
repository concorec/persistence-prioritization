# Persistence-Based Spatial Conservation Prioritization

This repository contains the analysis workflow for a conservation spatial-prioritization project that links trait-informed population persistence models to spatial prioritization of habitat in Madagascar. The workflow estimates species-specific persistence curves, converts distribution and habitat rasters into connected population units, and ranks habitat cells by their marginal contribution to multispecies persistence.

The current application focuses on Malagasy endemic mammals and birds and uses a persistence objective based on the probability that at least one connected population unit for a species remains above a quasi-extinction threshold over a 100-year time horizon.

## Repository notes

Several large raster and intermediate-output directories are intentionally omitted in this GitHub repository because of file-size limitations. To run the full pipeline, these files/directories must be restored locally:

- `bird_ppm_bin/`
- `bird_rangebag_bin/`
- `mammal_ppm_bin/`
- `mammal_rangebag_bin/`
- `Data/Clean/PriorityInputs/`
- `Data/Raw/esacci_2022_pfts.tif`
- Within each run directory under `Data/Results/PriorityRuns/`:
  - `patch_lookup_tables/`
  - `ana/persist_cmp_abf/`

The scripts use relative paths, so they should be run from the project root.

## Workflow overview

The scripts are numbered in approximate execution order. Scripts `1` through `4` build demographic and species-level inputs, scripts `5` and `5.1` build or illustrate spatial population-unit inputs, script `6` runs the persistence-based prioritization, and scripts `7.1` through `7.3` generate derived outputs and benchmark comparisons.

### `1_rm_sigma_models.Rmd`

**Purpose:** Fits trait-based demographic calibration models for intrinsic growth rate (`rm`) and environmental variability (`sigma_r`) for mammals and birds. Mammal `rm` and `sigma_r` are modeled as functions of body mass. Bird `rm` is modeled as a function of body mass, while bird `sigma_r` is modeled using body mass plus ILR-transformed diet composition. The script also creates diagnostic/manuscript figures in the knitted output.

**Required inputs:**

- `Data/Raw/mammal_rmax.txt`
- `Data/Raw/bird_data.txt`
- `Data/Raw/sigma.csv`
- `Data/Raw/bird_synonyms.csv`
- JAGS/rjags installation if `RUN_FITS = TRUE`

**Outputs:**

- `Data/Clean/post_mammal_rm_coefs.csv`
- `Data/Clean/post_mammal_sigma_coefs.csv`
- `Data/Clean/post_bird_r_coefs.csv`
- `Data/Clean/post_bird_sigma_ilr_coefs.csv`
- `Data/Clean/bird_ilr_combo_coords.csv`
- `Data/Clean/bird_mass_grid.csv`
- `Data/Clean/mammal_mass_grid.csv`

### `2_persist_points_crn.Rmd`

**Purpose:** Uses posterior draws from script 1 and a C++ stochastic population simulator to estimate persistence probability across equilibrium abundances (`K`). It uses common random numbers and writes tidy persistence-point tables for later Gompertz fitting. The simulation uses a 100-year horizon, quasi-extinction threshold of 500 individuals, 760 demographic draws, and 2,500 replicate trajectories per draw.

**Required inputs:**

- `Data/Clean/mammal_mass_grid.csv`
- `Data/Clean/bird_mass_grid.csv`
- `Data/Clean/bird_ilr_combo_coords.csv`
- `Data/Clean/post_mammal_rm_coefs.csv`
- `Data/Clean/post_mammal_sigma_coefs.csv`
- `Data/Clean/post_bird_r_coefs.csv`
- `Data/Clean/post_bird_sigma_ilr_coefs.csv`
- `src/simulate_persist_probs_cpp.cpp`

**Outputs:**

- `Data/Results/persist_points_mammals.csv`
- `Data/Results/persist_points_birds.csv`
- Temporary compiled-code cache in `Data/.rcpp_cache/`

### `3_gompertz_loess_models.Rmd`

**Purpose:** Fits shifted two-parameter Gompertz curves to the simulated persistence points from script 2, then fits LOESS models that interpolate Gompertz parameters over body mass for mammals and over body mass within bird diet-composition classes. Diagnostic figures are printed in the knitted document, including comparisons with the Wolff-style Gompertz specification.

**Required inputs:**

- `Data/Results/persist_points_mammals.csv`
- `Data/Results/persist_points_birds.csv`
- `Data/Clean/bird_ilr_combo_coords.csv` for some diagnostic bird plots
- `Data/Clean/post_bird_sigma_ilr_coefs.csv` for some diagnostic bird plots

**Outputs:**

- `Data/Clean/gompertz_loess_models.rds`

### `4_build_species_table.Rmd`

**Purpose:** Builds the canonical species-level table used by downstream spatial scripts. It joins species metadata, traits, raster-source information, density/home-range/dispersal quantities, minimum patch and population-unit size thresholds, Gompertz persistence parameters, and IUCN habitat information.

**Required inputs:**

- `Data/Raw/simple_summary.csv`
- `Data/Raw/synonyms.csv`
- `Data/Raw/mammal_data.txt`
- `Data/Raw/bird_data.txt`
- `Data/Raw/random_effects.csv`
- `Data/Clean/gompertz_loess_models.rds`
- Presence raster directories:
  - `mammal_ppm_bin/`
  - `mammal_rangebag_bin/`
  - `bird_ppm_bin/`
  - `bird_rangebag_bin/`
- `IUCN_REDLIST_KEY` environment variable for IUCN habitat queries

**Outputs:**

- `Data/Clean/species_table.csv`

### `5_build_patches_and_connectivity.Rmd`

**Purpose:** Converts species distributions and habitat masks into spatial population-unit inputs. For each selected species, the script intersects suitable habitat with the species distribution raster, delineates AOH patches using rook adjacency, removes patches below the species-specific minimum patch area, groups remaining patches into dispersal-connected population units, removes units below the minimum population-unit area, and writes patch rasters plus lookup/connectivity objects.

**Required inputs:**

- `Data/Clean/species_table.csv`
- `Data/Raw/esacci_2022_pfts.tif`
- Species distribution rasters referenced in `species_table.csv`, from the PPM and/or RangeBag raster directories
- GRASS GIS / `fasterRaster` for clumping; the script currently uses `C:/Program Files/GRASS GIS 8.4` as the GRASS path

**Outputs:**

- `Data/Clean/Patches/<species>.tif`
- `Data/Clean/all_patch_lookup.rds`
- `Data/Clean/all_connectivity.rds`
- Optional/convenience output from the later block:
  - `Data/Clean/Patches_binary/<species>.tif`

### `5.1_single_species_aoh_patches_pu_process.Rmd`

**Purpose:** Generates a single-species visual walkthrough of the AOH, patch, and population-unit construction workflow. The provided target species is `Cryptoprocta ferox`. This script is mainly for figure production and process illustration, not for generating canonical pipeline inputs.

**Required inputs:**

- `Data/Clean/species_table.csv`
- `Data/Raw/esacci_2022_pfts.tif`
- `mammal_ppm_bin/<target_raster>.tif`
- GRASS GIS / `fasterRaster`

**Outputs:**

- Printed figure in the knitted RMarkdown output
- Optional commented-out exports:
  - `Figures/cf_process_figure.png`
  - `Figures/cf_process_figure.pdf`

### `6_spatial_prioritization_pipeline.Rmd`

**Purpose:** Runs the persistence-based reverse-removal spatial prioritization. The script loads the species table, patch rasters, patch lookup table, and connectivity object; builds or loads a pruning-input bundle; scores frontier cells by marginal persistence loss; removes low-loss cells in batches; updates patch and population-unit state; and writes a nested removal-order raster and stage-level lookup tables.

**Required inputs:**

- `Data/Clean/species_table.csv`
- `Data/Clean/Patches/<species>.tif`
- `Data/Clean/all_patch_lookup.rds`
- `Data/Clean/all_connectivity.rds`
- `Data/Clean/PriorityInputs/pruning_inputs_curve_<curve>_taxa_<taxa>_sdm_<sdm>_redlist_<filter>.rds`

**Outputs:**

- `Data/Clean/PriorityInputs/pruning_inputs_curve_<curve>_taxa_<taxa>_sdm_<sdm>_redlist_<filter>.rds`
- `Data/Results/PriorityRuns/<run_id>/removal_order.tif`
- `Data/Results/PriorityRuns/<run_id>/removal_events.csv`
- `Data/Results/PriorityRuns/<run_id>/patch_lookup_tables/stage_patch_lookup_stage_####.csv`

**Important note:** the initialization/bundle-building section at the top of this RMarkdown file is marked `eval=FALSE` in the provided script. If the pruning-input bundle is not already present, run that setup block manually or change the chunk behavior before running the production pipeline section.

### `7.1_stage_meta_and_priority_surface.Rmd`

**Purpose:** Converts raw priority-run outputs into analysis-ready stage metadata and a normalized priority surface. `stage_meta.csv` summarizes cumulative cells/area removed and retained by stage. `priority_surface_prop_removed.tif` converts the removal-order raster into a 0-1 surface, where cells removed later have higher values and cells retained to the end receive value 1.

**Required inputs:**

- `Data/Clean/PriorityInputs/pruning_inputs_curve_q025_taxa_mammals_birds_sdm_ppm_redlist_include_lc.rds` or another matching pruning-input bundle
- `Data/Results/PriorityRuns/<run_id>/removal_order.tif`
- `Data/Results/PriorityRuns/<run_id>/removal_events.csv`
- `Data/Results/PriorityRuns/<run_id>/patch_lookup_tables/`

**Outputs:**

- `Data/Results/PriorityRuns/<run_id>/ana/stage_meta.csv`
- `Data/Results/PriorityRuns/<run_id>/ana/priority_surface_prop_removed.tif`

The YAML parameter `fp_name: stage_meta_fp.rds` is defined, but the current script does not appear to write that file.

### `7.2_build_rank_lut.Rmd`

**Purpose:** Builds benchmark lookup tables for a selected rank-map method, controlled by the YAML parameter `rank_method` (`"abf"`, `"caz"`, or `"cazmax"`). It matches the selected benchmark rank map to the same retained-area or retained-cell targets as the persistence-based stages, reconstructs retained landscapes at each stage, recalculates patches and population units, and writes stage-specific lookup tables for persistence evaluation.

**Required inputs:**

- Matching pruning-input bundle from `Data/Clean/PriorityInputs/`
- `Data/Results/PriorityRuns/<run_id>/ana/stage_meta.csv`
- `Data/Results/PriorityRuns/<run_id>/removal_order.tif` as the template grid
- `Data/Results/PriorityRuns/<run_id>/benchmark_rank_maps/rankmap_<rank_method>.tif`
  - for ABF: `Data/Results/PriorityRuns/<run_id>/benchmark_rank_maps/rankmap_abf.tif`
  - for CAZ: `Data/Results/PriorityRuns/<run_id>/benchmark_rank_maps/rankmap_caz.tif`- `Data/Clean/species_table.csv`
- `Data/Clean/Patches/<species>.tif`
- GRASS GIS / `fasterRaster`, unless using the `terra` clump backend

**Outputs:**

- `Data/Results/PriorityRuns/<run_id>/ana/rank_lut_<rank_method>/rank_stage_targets_<rank_method>.csv`
- `Data/Results/PriorityRuns/<run_id>/ana/rank_lut_<rank_method>/lut_stage_0000.rds`
- `Data/Results/PriorityRuns/<run_id>/ana/rank_lut_<rank_method>/lut_stage_####.rds`
- `Data/Results/PriorityRuns/<run_id>/ana/rank_lut_<rank_method>/rank_lut_inventory_<rank_method>.csv`

### `7.3_persist_cmp.Rmd`

**Purpose:** Evaluates and compares modeled persistence outcomes for the persistence-based priority run and the ABF benchmark. It reads the persistence-run patch lookup tables and the benchmark LUTs from script 7.2, standardizes them, computes PU-level persistence, species-level persistence, stage summaries, and method differences, and produces manuscript-style comparison plots in the knitted output.

**Required inputs:**

- Matching pruning-input bundle from `Data/Clean/PriorityInputs/`
- `Data/Clean/species_table.csv`
- `Data/Results/PriorityRuns/<run_id>/ana/stage_meta.csv`
- `Data/Results/PriorityRuns/<run_id>/patch_lookup_tables/stage_patch_lookup_stage_####.csv`
- `Data/Results/PriorityRuns/<run_id>/ana/rank_lut_abf/lut_stage_####.rds`

**Outputs:**

- `Data/Results/PriorityRuns/<run_id>/ana/persist_cmp_abf/pu_long_abf.csv`
- `Data/Results/PriorityRuns/<run_id>/ana/persist_cmp_abf/sp_long_abf.csv`
- `Data/Results/PriorityRuns/<run_id>/ana/persist_cmp_abf/stage_sum_abf.csv`
- `Data/Results/PriorityRuns/<run_id>/ana/persist_cmp_abf/stage_compare_abf.csv`
- Printed comparison figures in the knitted RMarkdown output

## Expected working directory structure

```text
.
├── 1_rm_sigma_models.Rmd
├── 2_persist_points_crn.Rmd
├── 3_gompertz_loess_models.Rmd
├── 4_build_species_table.Rmd
├── 5_build_patches_and_connectivity.Rmd
├── 5.1_single_species_aoh_patches_pu_process.Rmd
├── 6_spatial_prioritization_pipeline.Rmd
├── 7.1_stage_meta_and_priority_surface.Rmd
├── 7.2_build_rank_lut.Rmd
├── 7.3_persist_cmp.Rmd
├── src/
│   └── simulate_persist_probs_cpp.cpp
├── bird_ppm_bin/                  # large; omitted
├── bird_rangebag_bin/             # large; omitted
├── mammal_ppm_bin/                # large; omitted
├── mammal_rangebag_bin/           # large; omitted
├── Data/
│   ├── Raw/
│   │   ├── bird_data.txt
│   │   ├── bird_synonyms.csv
│   │   ├── mammal_data.txt
│   │   ├── mammal_rmax.txt
│   │   ├── random_effects.csv
│   │   ├── sigma.csv
│   │   ├── simple_summary.csv
│   │   ├── synonyms.csv
│   │   └── esacci_2022_pfts.tif   # large; omitted
│   ├── Clean/
│   │   ├── bird_ilr_combo_coords.csv
│   │   ├── bird_mass_grid.csv
│   │   ├── mammal_mass_grid.csv
│   │   ├── post_*_coefs.csv
│   │   ├── gompertz_loess_models.rds
│   │   ├── species_table.csv
│   │   ├── all_patch_lookup.rds
│   │   ├── all_connectivity.rds
│   │   ├── Patches/
│   │   │   └── <species>.tif
│   │   ├── Patches_binary/
│   │   │   └── <species>.tif
│   │   └── PriorityInputs/         # large; omitted
│   │       └── pruning_inputs_*.rds
│   └── Results/
│       ├── persist_points_mammals.csv
│       ├── persist_points_birds.csv
│       └── PriorityRuns/
│           └── <run_id>/
│               ├── benchmark_rank_maps/
│               │   ├── rankmap_abf.tif
│               │   └── rankmap_caz.tif
│               ├── removal_order.tif
│               ├── removal_events.csv
│               ├── patch_lookup_tables/       # large; omitted
│               │   └── stage_patch_lookup_stage_####.csv
│               └── ana/
│                   ├── stage_meta.csv
│                   ├── priority_surface_prop_removed.tif
│                   ├── rank_lut_abf/
│                   │   ├── rank_stage_targets_abf.csv
│                   │   ├── lut_stage_####.rds
│                   │   └── rank_lut_inventory_abf.csv
│                   └── persist_cmp_abf/        # large; omitted
│                       ├── pu_long_abf.csv
│                       ├── sp_long_abf.csv
│                       ├── stage_sum_abf.csv
│                       └── stage_compare_abf.csv
└── Figures/
    └── optional exported figures
```

## Typical execution order

1. Run `1_rm_sigma_models.Rmd` to fit demographic allometries and produce posterior coefficient draws.
2. Run `2_persist_points_crn.Rmd` to simulate persistence points for selected taxa/diet combinations.
3. Run `3_gompertz_loess_models.Rmd` to fit shifted-Gompertz curves and LOESS interpolation objects.
4. Run `4_build_species_table.Rmd` to assemble `Data/Clean/species_table.csv`.
5. Run `5_build_patches_and_connectivity.Rmd` to create patch rasters, patch lookup tables, and connectivity objects.
6. Run `6_spatial_prioritization_pipeline.Rmd` to create/load the pruning-input bundle and execute the priority-ranking algorithm.
7. Run `7.1_stage_meta_and_priority_surface.Rmd` to build analysis metadata and the priority surface.
8. Run `7.2_build_rank_lut.Rmd` to reconstruct benchmark ABF retained landscapes at matched stages.
9. Run `7.3_persist_cmp.Rmd` to compare persistence outcomes between the persistence-based prioritization and ABF benchmark.
10. Use `5.1_single_species_aoh_patches_pu_process.Rmd` independently when a visual example of the AOH-to-population-unit workflow is needed.

## Software and environment requirements

The workflow is written in R/RMarkdown and uses several spatial, Bayesian, and simulation packages. Required or frequently used packages include:

- `readr`, `data.table`, `dplyr`, `tibble`, `stringr`, `tools`
- `ggplot2`, `cowplot`, `scales`, `viridisLite`, `tidyterra`
- `terra`, `sf`, `igraph`, `units`, `rnaturalearth`
- `rjags` and a working JAGS installation
- `Rcpp`
- `fasterRaster` and GRASS GIS
- `rredlist`
- `zCompositions` for some diet-composition diagnostics

Before running `4_build_species_table.Rmd`, set an IUCN Red List API key in the environment:

```r
Sys.setenv(IUCN_REDLIST_KEY = "your_key_here")
```

On non-Windows systems, or if GRASS is installed somewhere else, update the hard-coded `grassDir` value in scripts `5`, `5.1`, and `7.2`.

## Output interpretation

The main prioritization output is `removal_order.tif`, which records when each initially occupied habitat cell was removed by the reverse-removal algorithm. `priority_surface_prop_removed.tif` rescales that removal order to the cumulative proportion of area or cells removed at the time of loss. Higher values correspond to cells retained until later in the sequence and therefore treated as higher priority under the fitted persistence objective.

The comparison outputs in `persist_cmp_abf/` summarize persistence through the removal sequence for both the persistence-based run and the ABF benchmark. `sp_long_abf.csv` stores species-level persistence trajectories, `pu_long_abf.csv` stores population-unit-level summaries, `stage_sum_abf.csv` stores method/stage summaries, and `stage_compare_abf.csv` stores side-by-side method differences.
