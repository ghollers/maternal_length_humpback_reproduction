# ============================================================
# Full pipeline, step 1: re-fit the seven empirical sub-models from the
# anonymized tag-derived tables in ../data/model_input_tables/ - no raw data.
#
# This reproduces the `models` object (equivalent to the fitted-model object
# used by the simulation). After this, step 2 runs the Monte Carlo simulation
# that produces the processed datasets used by Figures 1 and 3.
#
# Run: source("1_refit_models_from_tag_data.R")
# ============================================================

suppressMessages({ library(dplyr); library(lme4); library(glmmTMB) })

d  <- file.path("..", "data", "model_input_tables")
rd <- function(f) {
  p <- file.path(d, f)
  if (!file.exists(p))
    stop("Missing ", f, " in data/model_input_tables/. See ",
         "data/model_input_tables/README.md for what these tables are.")
  read.csv(p)
}

df_chord           <- rd("fluke_chord.csv")
df_area            <- rd("fluke_area.csv")
df_freq            <- rd("tailbeat_frequency.csv")
calf_speed_summary <- rd("swim_speed_calf.csv")
mom_speed_summary  <- rd("swim_speed_mother.csv")
mom_long_stroke    <- rd("prop_fluking_mother.csv")
calf_long_stroke   <- rd("prop_fluking_calf.csv")

# re-fit (identical formulas to the original analysis)
model_chord      <- lm(ChordLength ~ BodyLength, data = df_chord)
model_area       <- lm(FlukePlanarArea ~ BodyLength, data = df_area)
model_freq       <- lm(log10(mean_freq) ~ log10(BodyLength), data = df_freq)
model_speed_calf <- lmer(mean_speed ~ Behavior + (1 | ID), data = calf_speed_summary, REML = TRUE)
model_speed_mom  <- lmer(mean_speed ~ Behavior + (1 | ID), data = mom_speed_summary,  REML = TRUE)
mom_propstroke   <- glmmTMB(prop_sec ~ Behavior + (1 | ID), data = mom_long_stroke,  family = beta_family())
calf_propstroke  <- glmmTMB(prop_sec ~ Behavior + BodyLength + (1 | ID), data = calf_long_stroke, family = beta_family())

models <- list(
  propstroke_calf = calf_propstroke, propstroke_mom = mom_propstroke,
  speed_calf = model_speed_calf, speed_mom = model_speed_mom,
  freq_body = model_freq, chord_body = model_chord, area_body = model_area
)
assign("models", models, envir = globalenv())

message("Re-fit 7 sub-models from the anonymized tag tables into `models`.\n",
        "Next: step 2 (full_pipeline/2_run_simulation.R) runs the Monte Carlo ",
        "simulation that regenerates the processed datasets used by Figures 1 and 3.")
