# Reviewer reproducibility package

Streamlined data and code to reproduce the main figures of the manuscript on
maternal body length and the energetics of reproduction in humpback whales. Each
main figure has its own folder with a single, self-contained script. No raw
biologging data is required: figures reproduce from processed data, and the
model-fitting step reproduces from anonymized per-animal summaries.

## Layout

```
reviewer_package/
+-- reviewer_package.Rproj         <- open this first (sets the working directory)
+-- reproduce_all_figures.R        <- runs all four figure scripts
+-- data/
|   +-- processed/                 <- derived datasets the figures read
|   |   +-- lactation_estimates.csv
|   |   +-- departure_df.csv
|   |   +-- yates_CT.csv, yates_uprop.csv   <- foil lookup tables (used by full_pipeline)
|   +-- model_input_tables/        <- anonymized tag-derived tables (see its README)
+-- figure_1_lactation_cost/       <- make_figure.R + README
+-- figure_2_departure_and_reserves/
+-- figure_3_annual_time_budget/
+-- figure_4_population_growth/     <- make_figure.R + fecundity_function.R + README
+-- full_pipeline/                 <- tag tables -> refit models -> simulation -> figures
```

## Quick start

1. Open `reviewer_package.Rproj` in RStudio.
2. Install packages (once):
   `install.packages(c("dplyr","tidyr","ggplot2","tibble","truncnorm","patchwork","lme4","glmmTMB"))`
3. Reproduce everything: `source("reproduce_all_figures.R")`, or run any single
   figure by opening its folder and `source("make_figure.R")`.

Each figure writes a `figure_*.png` into its own folder. Figures 1 and 3 read a
processed CSV; Figures 2 and 4 are self-contained (deterministic growth,
allometry and demographic models from published/hardcoded coefficients).

## Two levels of reproducibility

- **Figures from processed data (one click):** the four `figure_*/make_figure.R`
  scripts. This is all a reviewer needs to reproduce the published figures.
- **The full chain from the tag data:** `full_pipeline/` reproduces how the
  processed datasets are produced, entirely within this package - anonymized tag
  tables -> re-fit the seven empirical sub-models (`1_refit_models_from_tag_data.R`)
  -> Monte Carlo simulation (`2_run_simulation.R`, using `full_pipeline/engine/`).
  See `full_pipeline/README.md`.


## Data and privacy

All data here are processed or aggregated. The `data/model_input_tables/` tables
are per-animal summaries with whale IDs replaced by anonymous codes - no raw
sensor streams, no locations, no dates (see `data/model_input_tables/README.md`).
The raw biologging recordings are owned by the Marine Mammal Research Program
(PI L. Bejder) and are not required to reproduce any result in the paper.

## Packages

Figures: dplyr, tidyr, ggplot2, tibble, truncnorm, patchwork.
Full pipeline (`full_pipeline/`): additionally purrr, MASS, broom.mixed, lme4,
glmmTMB, furrr, future, knitr.

## Licensing

Code (the `.R` and `.Rmd` files) is released under the MIT License (`LICENSE`).
Data (the `.csv` files under `data/`) are released under Creative Commons
Attribution 4.0 International (CC-BY 4.0; `LICENSE-data.txt`). Fill in the author
name in both files before depositing.
