# ============================================================
# Full pipeline, step 2: run the Monte Carlo energy simulation and regenerate the
# processed datasets used by Figures 1 and 3 - entirely within reviewer_package,
# from the anonymized tag tables. No raw sensor data, no external repository.
#
# This:
#   1. sources the simulation engine (full_pipeline/engine/energetics_functions.R)
#   2. re-fits the 7 sub-models from data/model_input_tables/ (step 1)
#   3. builds the coefficient/covariance objects the simulation needs
#   4. runs the simulation (engine/02) and lactation step (engine/03), which write
#      data/processed/simulation_results.* , full_table.* and lactation_estimates.*
#
# NOTE: step 1 (the Monte Carlo) is slow. Open reviewer_package.Rproj first, then
#   source("full_pipeline/2_run_simulation.R")
# ============================================================

suppressMessages({
  library(dplyr); library(purrr); library(tidyr)
  library(lme4); library(glmmTMB); library(broom.mixed)
  library(furrr); library(future); library(knitr)
})

PKG <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
engine <- file.path(PKG, "full_pipeline", "engine")

# 1. engine functions
source(file.path(engine, "energetics_functions.R"))

# 2. re-fit the seven sub-models from the anonymized tag tables -> `models`
source(file.path(PKG, "full_pipeline", "1_refit_models_from_tag_data.R"))

# 3. derived objects the simulation expects in the global environment
out                 <- purrr::imap(models, extract_model)
fixed_effects_table <- bind_rows(lapply(out, `[[`, "fixed"))
vcov_list           <- lapply(out, `[[`, "vcov"); names(vcov_list) <- names(models)
calf_propstroke     <- models$propstroke_calf
mom_propstroke      <- models$propstroke_mom

# 4. run the simulation and lactation modules (they read/write data/processed/)
run_module <- function(rmd) {
  rf <- knitr::purl(file.path(engine, rmd), output = tempfile(fileext = ".R"), quiet = TRUE)
  sys.source(rf, envir = globalenv())
}
run_module("02_dee_simulation.Rmd")    # -> simulation_results.*, full_table.*
run_module("03_maternal_supply.Rmd")   # -> lactation_estimates.*

message("Done. Regenerated data/processed/ from the anonymized tag tables. ",
        "Figures 1 and 3 now reproduce from these freshly simulated datasets.")
