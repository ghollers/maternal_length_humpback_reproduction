# Full pipeline: from anonymized tag data to the figures

This folder reproduces the entire chain behind the figures

```
anonymized tag tables              data/model_input_tables/*.csv   [provided]
        |  1_refit_models_from_tag_data.R
        v
7 re-fitted sub-models             (in memory: `models`)
        |  2_run_simulation.R  (engine/02 + engine/03)
        v
processed datasets                 data/processed/simulation_results.*,
        |                          full_table.*, lactation_estimates.*
        v
figures 1-4                        ../figure_1_.. etc.
```

The anonymized tag tables in `../data/model_input_tables/` are **provided with
this package**; they were produced once from the raw data in the authors' full
analysis project and carry no individual-animal information. Nothing in this
package needs the raw recordings.

## Steps (optional - only to verify the pipeline from the tables)

1. **Re-fit the sub-models**: `source("1_refit_models_from_tag_data.R")`
   rebuilds the seven fitted models from `../data/model_input_tables/` with the
   study's formulas.

2. **Run the simulation**: `source("2_run_simulation.R")` runs the Monte Carlo
   energy simulation (`engine/02_dee_simulation.Rmd`) and the lactation step
   (`engine/03_maternal_supply.Rmd`), regenerating the processed datasets in
   `../data/processed/` that Figures 1 and 3 read. (This is slow.)

The processed datasets are already provided in `../data/processed/`, so Figures
1-4 reproduce without re-running any of this. Steps 1-2 are here so the whole
pipeline, including the empirical model fitting and the simulation, can be
verified end-to-end from the non-sensitive tables.

## Contents

- `engine/energetics_functions.R` - all model/simulation helper functions
- `engine/02_dee_simulation.Rmd` - thrust/metabolic power + Monte Carlo
- `engine/03_maternal_supply.Rmd` - lactation cost from the simulation output
- `1_refit_models_from_tag_data.R` - re-fit the 7 sub-models from the tables
- `2_run_simulation.R` - orchestrates steps 1-2 and writes the processed data
