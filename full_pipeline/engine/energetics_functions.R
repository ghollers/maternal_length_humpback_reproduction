# energetics_functions.R
# Single source of truth for all custom functions in the maternal-energetics pipeline.
# GENERATED from analysis/*.Rmd -- if you change a function, change it here.
# Sourced at the top of every module (and 00_setup) so functions are available
# without re-running the heavy model-fitting or simulation stages.
#
# EXTERNAL MODEL PROVENANCE
# Several growth and body-composition relationships used below come from two
# companion studies by van Aswegen et al.:
#
#   [1] Age-specific growth costs
#       van Aswegen M, Szabo A, Currie JJ, Stack SH, Straley J, Neilson J,
#       Gabriele C, Cates K, Baker CS, Steel D, Pack AA, West KL, Hofmann NJ,
#       Toonen RJ, Loockerman C, Sullivan FA, Cheeseman T, Sharpe F, Bejder L
#       (2025) Age-specific body length, mass, and energetic cost of growth in
#       humpback whales. Mar Ecol Prog Ser 770:171-194. doi:10.3354/meps14937
#       -> calf and maternal length->mass and length-at-age relationships
#          (e.g. predict_calf_mass, predict_mother_mass, growth terms) and the
#          calf energy density E_es.
#
#   [2] Maternal investment
#       van Aswegen M, Szabo A, Currie JJ, Stack SH, Evans L, Straley J,
#       Neilson J, Gabriele C, Cates K, Steel D, ... Bejder L (2025) Maternal
#       investment, body condition and calf growth in humpback whales.
#       J Physiol 603(2):551-578. doi:10.1113/JP287379
#       -> daily body-volume loss and its conversion to energy
#          (predict_daily_bv_loss, predict_mother_energy_loss[_v2]).
#
# Hardcoded body-composition and tissue-density coefficients also follow these
# two studies, with Nishiwaki (1959) and Lockyer (1976) for tissue proportions
# and body-fluid fraction.

suppressMessages({
  library(dplyr); library(purrr); library(tidyr); library(MASS)
  library(glmmTMB); library(lme4); library(broom.mixed)
  library(performance); library(truncnorm)
})


# ================= from 01_empirical_models =================

get_rest_time_exact <- function(dives, start_index, end_index) {
  
  # create full index range
  idx <- rep(FALSE, end_index - start_index + 1)
  
  for (i in 1:nrow(dives)) {
    
    s <- dives$StartIndex[i]
    e <- dives$EndIndex[i]
    
    # clip just in case
    s <- max(s, start_index)
    e <- min(e, end_index)
    
    if (e >= s) {
      idx[(s - start_index + 1):(e - start_index + 1)] <- TRUE
    }
  }
  
  sum(idx) / 10
}

get_sigma <- function(model) {
  if (inherits(model, c("lmerMod", "lmerModLmerTest", "lm"))) {
    return(sigma(model))
  } else {
    return(NA)
  }
}

get_vcov_matrix <- function(model) {
  tryCatch({
    vcov(model)
  }, error = function(e) NA)
}

get_random_var <- function(model) {
  if (!inherits(model, c("lmerMod", "lmerModLmerTest"))) {
    return(NA)
  }
  
  vc <- as.data.frame(VarCorr(model))
  
  # assumes first grouping factor = ID
  data.frame(
    var_ID = vc$vcov[1],
    sd_ID  = sqrt(vc$vcov[1])
  )
}

get_icc <- function(model) {
  
  if (!inherits(model, c("lmerMod", "lmerModLmerTest"))) {
    return(NA)  # lm models
  }
  
  icc_val <- tryCatch(performance::icc(model), error = function(e) NA)
  
  if (is.na(icc_val)[1]) return(NA)
  
  if (is.data.frame(icc_val)) {
    return(icc_val$ICC_adjusted[1])
  }
  
  if (is.numeric(icc_val)) {
    return(unname(icc_val[1]))
  }
  
  NA
}

make_equation <- function(model_name) {
  
  eq <- case_when(
    
    model_name == "stroke_calf" ~ 
      "StrokeRate ~ Behavior * BodyLength + (1 | ID)",
    
    model_name == "stroke_mom" ~ 
      "StrokeRate ~ Behavior + (1 | ID)",
    
    model_name == "speed_calf" ~ 
      "MeanSpeed ~ Behavior + (1 | ID)",
    
    model_name == "speed_mom" ~ 
      "MeanSpeed ~ Behavior + (1 | ID)",
    
    model_name == "freq_body" ~ 
      "log10(MeanFreq) ~ log10(BodyLength)",
    
    model_name == "chord_body" ~ 
      "ChordLength ~ BodyLength",
    
    model_name == "area_body" ~ 
      "FlukePlanarArea ~ BodyLength",
    
    TRUE ~ NA_character_
  )
  
  return(eq)
}

extract_model <- function(model, model_name) {
  
  is_lmer <- inherits(model, c("lmerMod", "lmerModLmerTest"))
  
  fixed <- if (is_lmer) {
    broom.mixed::tidy(model, effects = "fixed")
  } else {
    broom::tidy(model)
  }
  
  vcov_mat <- get_vcov_matrix(model)
  
  stats <- data.frame(
    model = model_name,
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    n = nobs(model),
    ICC = get_icc(model),
    sigma = get_sigma(model),
    equation = make_equation(model_name),
    rand_SD_ID = if (is_lmer) sqrt(as.data.frame(VarCorr(model))$vcov[1]) else NA,
    var_ID = if (is_lmer) as.data.frame(VarCorr(model))$vcov[1] else NA
  )
  
  fixed$model <- model_name
  
  list(
    fixed = fixed,
    stats = stats,
    vcov = vcov_mat
  )
}

# Residual SD source. In the full pipeline this is NULL, so get_resid_sigma()
# reads the residual SD from the fitted model object. In parameters-only mode,
# load_model_parameters.R populates model_sigmas and the model objects are not
# needed. (See R/export_model_parameters.R / R/load_model_parameters.R.)
model_sigmas <- NULL

get_resid_sigma <- function(name) {
  if (!is.null(model_sigmas) && !is.null(model_sigmas[[name]])) return(model_sigmas[[name]])
  sigma(models[[name]])
}

draw_coefs <- function(model_name, n) {
  
  beta <- fixed_effects_table %>%
    dplyr::filter(model == model_name) %>%
    dplyr::arrange(match(term, names(vcov_list[[model_name]]))) %>%
    dplyr::pull(estimate)
  
  vcov_mat <- vcov_list[[model_name]]
  
  MASS::mvrnorm(n = n, mu = beta, Sigma = vcov_mat)
}

predict_chord <- function(body_length, n = 1000) {
  
  coefs <- draw_coefs("chord_body", n)
  
  intercept <- coefs[, "(Intercept)"]
  slope     <- coefs[, "BodyLength"]
  
  mu <- intercept + slope * body_length
  
  sigma <- get_resid_sigma("chord_body")
  
  rnorm(n, mean = mu, sd = sigma)
}

predict_area <- function(body_length, n = 1000) {
  
  coefs <- draw_coefs("area_body", n)
  
  intercept <- coefs[, "(Intercept)"]
  slope     <- coefs[, "BodyLength"]
  
  mu <- intercept + slope * body_length
  
  sigma <- get_resid_sigma("area_body")
  
  rnorm(n, mean = mu, sd = sigma)
}

predict_frequency <- function(body_length, n = 1000) {
  
  coefs <- draw_coefs("freq_body", n)
  
  intercept <- coefs[, "(Intercept)"]
  slope     <- coefs[, "log10(BodyLength)"]
  
  log_mu <- intercept + slope * log10(body_length)
  
  sigma <- get_resid_sigma("freq_body")
  
  log_freq <- rnorm(n, mean = log_mu, sd = sigma)
  
  10^log_freq
}

predict_speed_calf <- function(behavior, n = 1000) {
  
  coefs <- draw_coefs("speed_calf", n)
  
  intercept <- coefs[, "(Intercept)"]
  beta_rest <- coefs[, "Behaviorrest"]
  
  mu <- ifelse(
    behavior == "rest",
    intercept + beta_rest,
    intercept
  )
  
  sigma <- get_resid_sigma("speed_calf")
  
  rnorm(n, mean = mu, sd = sigma)
}

predict_speed_mom <- function(behavior, n = 1000) {
  
  coefs <- draw_coefs("speed_mom", n)
  
  intercept <- coefs[, "(Intercept)"]
  beta_rest <- coefs[, "Behaviorrest"]
  
  mu <- ifelse(
    behavior == "rest",
    intercept + beta_rest,
    intercept
  )
  
  sigma <- get_resid_sigma("speed_mom")
  
  rnorm(n, mean = mu, sd = sigma)
}

predict_speed <- function(role, behavior, n = 1000) {
  
  if (role == "calf") {
    predict_speed_calf(behavior, n)
  } else if (role == "mother") {
    predict_speed_mom(behavior, n)
  } else {
    stop("role must be 'calf' or 'mother'")
  }
}

predict_prop_fluking <- function(model, behavior, n = 200) {
  
  if (inherits(model, "prop_params")) {
    beta <- model$beta; vcov_mat <- model$vcov
  } else {
    beta    <- fixef(model)$cond
    vcov_mat <- vcov(model)$cond
  }
  
  coef_draws <- MASS::mvrnorm(n, mu = beta, Sigma = vcov_mat)
  
  intercept  <- coef_draws[, "(Intercept)"]
  beta_rest  <- if ("Behaviorrest" %in% colnames(coef_draws)) 
                  coef_draws[, "Behaviorrest"] else 0
  
  eta <- if (behavior == "rest") intercept + beta_rest else intercept
  
  # inverse logit (beta family uses logit link)
  plogis(eta)
}

predict_prop_fluking_calf <- function(model, behavior, body_length, n = 200) {
  
  if (inherits(model, "prop_params")) {
    beta <- model$beta; vcov_mat <- model$vcov
  } else {
    beta     <- fixef(model)$cond
    vcov_mat <- vcov(model)$cond
  }
  
  coef_draws <- MASS::mvrnorm(n, mu = beta, Sigma = vcov_mat)
  
  intercept   <- coef_draws[, "(Intercept)"]
  beta_rest   <- if ("Behaviorrest" %in% colnames(coef_draws)) 
                   coef_draws[, "Behaviorrest"] else 0
  beta_BL     <- if ("BodyLength" %in% colnames(coef_draws)) 
                   coef_draws[, "BodyLength"] else 0
  
  eta <- intercept + beta_BL * body_length +
         if (behavior == "rest") beta_rest else 0
  
  plogis(eta)
}


# ================= from 02_dee_simulation =================

calc_reduced_freq <- function(f, chord, U) {
  omega <- 2 * pi * f
  (omega * chord) / U
}

calc_feathering <- function(U, f, body_length,
                            alpha_deg = 20,   # max pitch angle (deg); baseline 20
                            amp_frac  = 0.1) { # fluke amplitude as fraction of body length; baseline 0.1
  alpha <- alpha_deg * pi / 180
  h <- amp_frac * body_length
  omega <- 2 * pi * f
  (U * alpha) / (omega * h)
}

get_ct_coeffs <- function(feathering, yatesCT) {
  
  rounded <- round(feathering, 2)
  
  # exact match
  idx <- which(round(yatesCT$Feathering, 2) == rounded)
  
  if (length(idx) > 0) {
    return(yatesCT[idx[1], ])
  }
  
  # fallback if feathering >= 0.80
  if (rounded >= 0.80) {
    return(yatesCT[nrow(yatesCT), ])
  }
  
  # last resort: nearest bin (rare)
  idx <- which.min(abs(yatesCT$Feathering - rounded))
  yatesCT[idx, ]
}

calc_Ct <- function(x, coeffs) {
  coeffs$a * x^3 + coeffs$b * x^2 + coeffs$c * x + coeffs$d
}

calc_Ct_vectorized <- function(U, f, body_length, chord, yatesCT,
                               alpha_deg = 20,    # max pitch angle (deg); baseline 20
                               amp_frac  = 0.1) { # fluke amplitude as fraction of body length; baseline 0.1

  n <- length(U)

  # ensure vectors
  body_length <- rep(body_length, length.out = n)
  chord       <- rep(chord, length.out = n)

  alpha <- alpha_deg * pi / 180
  h <- amp_frac * body_length
  omega <- 2 * pi * f
  feather <- (U * alpha) / (omega * h)
  
  # select Ct coefficients per MATLAB rules
  coeffs <- lapply(feather, get_ct_coeffs, yatesCT = yatesCT)
  coeffs <- do.call(rbind, coeffs)
  
  x <- (omega * chord) / U
  
  Ct <- coeffs$a * x^3 + coeffs$b * x^2 + coeffs$c * x + coeffs$d
  
  return(Ct)
}

calc_thrust_power_full <- function(U, f, body_length, chord, area, yatesCT,
                                   water_density = 1025, # seawater density (kg/m3); baseline 1025
                                   alpha_deg     = 20,   # max pitch angle (deg); baseline 20
                                   amp_frac      = 0.1) { # fluke amplitude as fraction of body length; baseline 0.1

  Ct <- calc_Ct_vectorized(U, f, body_length, chord, yatesCT,
                           alpha_deg = alpha_deg, amp_frac = amp_frac)
  h  <- amp_frac * body_length
  PT <- 0.5 * water_density * Ct * (U^3) * area * (h / chord)^2

  return(PT)
}

get_uprop_coeffs <- function(feathering, yatesUProp) {
  
  rounded <- round(feathering, 2)
  
  # exact match (MATLAB behavior)
  idx <- which(round(yatesUProp$Feathering, 2) == rounded)
  
  if (length(idx) > 0) {
    return(yatesUProp[idx[1], ])
  }
  
  # fallback regime (high feathering)
  if (rounded >= 0.80) {
    return(yatesUProp[nrow(yatesUProp), ])
  }
  
  # fallback nearest bin (rare edge case)
  idx <- which.min(abs(yatesUProp$Feathering - rounded))
  yatesUProp[idx, ]
}

calc_u_prop_vectorized <- function(U, f, body_length, chord, yatesUProp,
                                   alpha_deg = 20,    # max pitch angle (deg); baseline 20
                                   amp_frac  = 0.1) { # fluke amplitude as fraction of body length; baseline 0.1

  n <- length(U)

  body_length <- rep(body_length, length.out = n)
  chord       <- rep(chord, length.out = n)

  alpha <- alpha_deg * pi / 180
  h <- amp_frac * body_length

  omega <- 2 * pi * f
  
  feather <- (U * alpha) / (omega * h)
  
  # MATLAB-style coefficient selection
  coeffs <- lapply(feather, get_uprop_coeffs, yatesUProp = yatesUProp)
  coeffs <- do.call(rbind, coeffs)
  
  x <- (omega * chord) / U
  
  u_prop <- coeffs$a * x^3 + coeffs$b * x^2 + coeffs$c * x + coeffs$d
  
  # clamp to physical bounds
  u_prop <- pmin(pmax(u_prop, 0.01), 0.99)
  
  return(u_prop)
}

calc_metabolic_power_full <- function(U, f, body_length, chord, area,
                                      yatesCT, yatesUProp, n,
                                      u_met         = 0.25, # muscle metabolic efficiency; baseline 0.25 (Potvin et al. 2021)
                                      water_density = 1025, # seawater density (kg/m3); baseline 1025
                                      alpha_deg     = 20,   # max pitch angle (deg); baseline 20
                                      amp_frac      = 0.1) { # fluke amplitude as fraction of body length; baseline 0.1

  PT     <- calc_thrust_power_full(U, f, body_length, chord, area, yatesCT,
                                   water_density = water_density,
                                   alpha_deg = alpha_deg, amp_frac = amp_frac)
  u_prop <- calc_u_prop_vectorized(U, f, body_length, chord, yatesUProp,
                                   alpha_deg = alpha_deg, amp_frac = amp_frac)
  u_met  <- rep(u_met, n)   # baseline 0.25 (Potvin et al. 2021)

  Pmet <- PT / (u_met * u_prop)

  return(Pmet)
}

calculate_relative_length <- function(age) {
  tibble(
    Age = age,
    Calf_RL = 0.737106 - (0.737106 - 0.3504304) * exp(-0.006740519 * age),
    Calf_RL_Lwr = 0.7250429 - (0.7250429 - 0.3471785) * exp(-0.006297628 * age),
    Calf_RL_Upr = 0.7503627 - (0.7503627 - 0.3536075) * exp(-0.007190727 * age)
  )
}

convert_relative_to_total_length <- function(data) {
  data %>%
    mutate(
      Calf_TL = Calf_RL * ML,
      Calf_TL_Lwr = Calf_RL_Lwr * ML,
      Calf_TL_Upr = Calf_RL_Upr * ML
    )
}

calculate_body_volume <- function(data) {
  data$BV <- exp(-3.991 + 2.881 * log(data$Calf_TL))
  data$BV_Lwr <- exp(-3.991 + 2.881 * log(data$Calf_TL_Lwr))
  data$BV_Upr <- exp(-3.991 + 2.881 * log(data$Calf_TL_Upr))
  
  return(data)
}

calculate_tissue_masses <- function(data) {
  # Mass calculations (mean, lower, upper CI)
  data$Musc_Mass <- data$BV * (exp(1.934 + 0.657 * log(data$Calf_TL)) / 100) * 972.6259
  data$Bone_Mass <- data$BV * (exp(4.441 - 0.705 * log(data$Calf_TL)) / 100) * 1123.534
  data$Visc_Mass <- data$BV * (exp(3.276 - 0.369 * log(data$Calf_TL)) / 100) * 953.707
  data$Blood_Mass <- data$BV * 0.06 * 1048.657
  data$Blubber_Mass <-((1 + 0) * data$BV - data$Musc_Mass/972.6259 - data$Bone_Mass/1123.534 - data$Visc_Mass/953.707-data$Blood_Mass /1048.657)*986
  # lower CI.
  data$Musc_Mass_Lwr<-(data$BV_Lwr*(exp(1.934+0.657*log(data$Calf_TL_Lwr))/100))*972.6259 
  data$Bone_Mass_Lwr<-(data$BV_Lwr*(exp(4.441-0.705*log(data$Calf_TL_Lwr))/100))*1123.534  
  data$Visc_Mass_Lwr<-(data$BV_Lwr*(exp(3.276-0.369*log(data$Calf_TL_Lwr))/100))*953.707
  data$Blood_Mass_Lwr<-(data$BV_Lwr*0.06)*1048.657 
  data$Blubber_Mass_Lwr<-((1+0)*data$BV_Lwr-data$Musc_Mass_Lwr/972.6259-data$Bone_Mass_Lwr/1123.534-data$Visc_Mass_Lwr/953.707-data$Blood_Mass_Lwr/1048.657)*986
  # Upper CI
  data$Musc_Mass_Upr<-(data$BV_Upr*(exp(1.934+0.657*log(data$Calf_TL_Upr))/100))*972.6259 
  data$Bone_Mass_Upr<-(data$BV_Upr*(exp(4.441-0.705*log(data$Calf_TL_Upr))/100))*1123.534  
  data$Visc_Mass_Upr<-(data$BV_Upr*(exp(3.276-0.369*log(data$Calf_TL_Upr))/100))*953.707
  data$Blood_Mass_Upr<-(data$BV_Upr*0.06)*1048.657 
  data$Blubber_Mass_Upr<-((1+0)*data$BV_Upr-data$Musc_Mass_Upr/972.6259-data$Bone_Mass_Upr/1123.534-data$Visc_Mass_Upr/953.707-data$Blood_Mass_Upr/1048.657)*986
  
  return(data)
}

calculate_total_mass <- function(data) {
  data$Mass_kg <- data$Musc_Mass + data$Bone_Mass + data$Visc_Mass + data$Blubber_Mass + data$Blood_Mass
  data$Mass_kg_Lwr <- data$Musc_Mass_Lwr + data$Bone_Mass_Lwr + data$Visc_Mass_Lwr + data$Blubber_Mass_Lwr + data$Blood_Mass_Lwr
  data$Mass_kg_Upr <- data$Musc_Mass_Upr + data$Bone_Mass_Upr + data$Visc_Mass_Upr + data$Blubber_Mass_Upr + data$Blood_Mass_Upr
  
  return(data)
}

predict_mother_mass <- function(BL) {
  BM <- exp(3.069) * BL^2.825
  return(BM)
}

draw_lmer_coefs <- function(model, n) {
  
  beta <- lme4::fixef(model)
  vcov_mat <- vcov(model)
  
  draws <- MASS::mvrnorm(n = n, mu = beta, Sigma = vcov_mat)
  
  #  FIX: force matrix
  if (n == 1) {
    draws <- matrix(draws, nrow = 1)
    colnames(draws) <- names(beta)
  }
  
  return(draws)
}

simulate_daily_energy <- function(row, p_rest, n = 500,
                                  yatesCT, yatesUProp,
                                  calf_propstroke,
                                  mom_propstroke) {

  calf_BL <- as.numeric(row$Calf_TL)
  mom_BL  <- as.numeric(row$MotherLength)

  # --- TIME ALLOCATION ---
  n_rest   <- 1440 * p_rest
  n_active <- 1440 - n_rest

  # =========================
  # MORPHOLOGY
  # =========================
  calf_chord <- pmax(as.vector(predict_chord(calf_BL, n)), 0.05)
  calf_area  <- pmax(as.vector(predict_area(calf_BL, n)),  0.05)
  mom_chord  <- pmax(as.vector(predict_chord(mom_BL, n)),  0.05)
  mom_area   <- pmax(as.vector(predict_area(mom_BL, n)),   0.05)

  # =========================
  # LOCOMOTION (your unique data)
  # =========================
  U_r_calf <- predict_speed("calf",   "rest",  n)
  U_r_mom  <- predict_speed("mother", "rest",  n)
  U_a_calf <- predict_speed("calf",   "other", n)
  U_a_mom  <- predict_speed("mother", "other", n)

  f_r_calf <- pmax(as.vector(predict_frequency(calf_BL, n)), 0.01)
  f_r_mom  <- pmax(as.vector(predict_frequency(mom_BL,  n)), 0.01)
  f_a_calf <- pmax(as.vector(predict_frequency(calf_BL, n)), 0.01)
  f_a_mom  <- pmax(as.vector(predict_frequency(mom_BL,  n)), 0.01)

  P_r_calf <- calc_metabolic_power_full(U_r_calf, f_r_calf,
                                        calf_BL, calf_chord, calf_area,
                                        yatesCT, yatesUProp, n)
  P_r_mom  <- calc_metabolic_power_full(U_r_mom, f_r_mom,
                                        mom_BL, mom_chord, mom_area,
                                        yatesCT, yatesUProp, n)
  P_a_calf <- calc_metabolic_power_full(U_a_calf, f_a_calf,
                                        calf_BL, calf_chord, calf_area,
                                        yatesCT, yatesUProp, n)
  P_a_mom  <- calc_metabolic_power_full(U_a_mom, f_a_mom,
                                        mom_BL, mom_chord, mom_area,
                                        yatesCT, yatesUProp, n)

  prop_fluke_r_calf <- predict_prop_fluking_calf(calf_propstroke, "rest",  calf_BL, n)
  prop_fluke_a_calf <- predict_prop_fluking_calf(calf_propstroke, "other", calf_BL, n)
  prop_fluke_r_mom  <- predict_prop_fluking(mom_propstroke, "rest",  n)
  prop_fluke_a_mom  <- predict_prop_fluking(mom_propstroke, "other", n)

  # MJ/day
  loco_calf <- pmax(
    ((n_rest   * 60 * prop_fluke_r_calf * P_r_calf) +
     (n_active * 60 * prop_fluke_a_calf * P_a_calf)) / 1e6, 0)

  loco_mom <- pmax(
    ((n_rest   * 60 * prop_fluke_r_mom * P_r_mom) +
     (n_active * 60 * prop_fluke_a_mom * P_a_mom)) / 1e6, 0)

  # =========================
  # BMR - Kleiber point estimate (FIXED; no random +-10%)
  # The +-10% BMR uncertainty is no longer sampled in the Monte Carlo; it is
  # evaluated separately in the sensitivity analysis (bmr_mult), consistent with
  # the treatment of the other fixed literature constants.
  # =========================
  bmr_calf <- 0.04843 * 70 * (row$Mass_kg      ^ 0.75) * 86400 / 1e6
  bmr_mom  <- 0.04843 * 70 * (row$MotherMass_kg ^ 0.75) * 86400 / 1e6

  # =========================
  # HIF - calf only, as proportion of BMR
  # FIXED at 0.08 (midpoint of 0.04-0.12; ~Fahlman et al. 2025 dolphin ~8.2%).
  # The 0.04-0.12 range is no longer sampled in the Monte Carlo; it is evaluated
  # in the sensitivity analysis, consistent with the other fixed literature
  # constants (pinniped 4-10%; extended to 4-12% for mysticete uncertainty).
  # =========================
  hif_coef <- 0.08
  hif_calf <- hif_coef * bmr_calf

  # =========================
  # TOTALS
  # note: growth costs added separately from the van Aswegen et al. growth models (age-specific growth costs)
  # =========================
  # calf: locomotion + BMR + HIF (growth added post-hoc)
  total_calf <- loco_calf + bmr_calf + hif_calf

  # mother: locomotion + BMR
  # lactation costs added post-hoc once supply side is built
  total_mom  <- loco_mom + bmr_mom

  # =========================
  # RETURN - summarize from joint distribution
  # keep components for energy partitioning figures
  # =========================
  return(data.frame(
    MotherLength  = mom_BL,
    MotherMass_kg = row$MotherMass_kg,
    CalfLength    = calf_BL,
    CalfMass_kg   = row$Mass_kg,
    Age           = row$Age,
    p_rest        = p_rest,

    # calf components
    calf_loco_mean = mean(loco_calf),
    calf_loco_lwr  = quantile(loco_calf, 0.025),
    calf_loco_upr  = quantile(loco_calf, 0.975),

    calf_BMR_mean  = mean(bmr_calf),
    calf_BMR_lwr   = quantile(bmr_calf, 0.025),
    calf_BMR_upr   = quantile(bmr_calf, 0.975),

    calf_HIF_mean  = mean(hif_calf),
    calf_HIF_lwr   = quantile(hif_calf, 0.025),
    calf_HIF_upr   = quantile(hif_calf, 0.975),

    # calf total excluding growth
    calf_total_mean = mean(total_calf),
    calf_total_lwr  = quantile(total_calf, 0.025),
    calf_total_upr  = quantile(total_calf, 0.975),

    # mother components
    mom_loco_mean  = mean(loco_mom),
    mom_loco_lwr   = quantile(loco_mom, 0.025),
    mom_loco_upr   = quantile(loco_mom, 0.975),

    mom_BMR_mean   = mean(bmr_mom),
    mom_BMR_lwr    = quantile(bmr_mom, 0.025),
    mom_BMR_upr    = quantile(bmr_mom, 0.975),

    mom_total_mean = mean(total_mom),
    mom_total_lwr  = quantile(total_mom, 0.025),
    mom_total_upr  = quantile(total_mom, 0.975)
  ))
}


# ================= from 03_maternal_supply =================

Delta_MV_HI <- function(L,alpha = 0.056197, beta = -0.01262){
  delta_MV = alpha + L*beta
  return(delta_MV) #m3 day-1
}

get_mother_energy_loss_from_BC <- function(
    mother_lengths,      # vector of maternal lengths (m)
    doy_arrival  = 1,  # day of year arriving Hawaii
    doy_departure = 100, # day of year leaving Hawaii
    rho_b  = 986,        # blubber density kg/m3 (van Aswegen post-mortem calves, 986 +- 24.3)
    E_es   = 23761.12,   # energy density storage tissue kJ/kg
    zeta_s = 0.7,        # storage efficiency
    n      = 500         # Monte Carlo draws
  ) {
  
  days_in_hawaii <- doy_departure - doy_arrival
  
  purrr::map_dfr(mother_lengths, function(L) {
    
    # --- BC at arrival and departure (non-lactating baseline) ---
    BC_arrival   <- BC_pred(doy = doy_arrival,   rest = 0)
    BC_departure <- BC_pred(doy = doy_departure,  rest = 0)
    
    # --- body volumes ---
    vol_arrival   <- get_volume(BC = BC_arrival,   L = L)
    vol_departure <- get_volume(BC = BC_departure, L = L)
    
    # --- total volume lost over stay ---
    delta_vol <- vol_arrival - vol_departure  # m3
    
    # --- convert to energy lost per day (kJ/day -> MJ/day) ---
    # instantaneous rate from Delta_MV_HI
    daily_vol_loss <- Delta_MV_HI(L)          # m3/day
    
    # energy mobilized per day
    EMR_daily_kJ <- (daily_vol_loss * rho_b * E_es) / zeta_s
    EMR_daily_MJ <- EMR_daily_kJ / 1000
    
    # total energy lost over full stay
    total_energy_loss_MJ <- (delta_vol * rho_b * E_es) / (zeta_s * 1000)
    
    data.frame(
      MotherLength         = L,
      days_in_hawaii       = days_in_hawaii,
      BC_arrival           = BC_arrival,
      BC_departure         = BC_departure,
      delta_vol_m3         = delta_vol,
      EMR_daily_MJ         = EMR_daily_MJ,       # energy loss rate per day
      total_energy_loss_MJ = total_energy_loss_MJ # total over stay
    )
  })
}

predict_daily_bv_loss <- function(TL, n = 500) {

  # Model 2: coefficient uncertainty (independent normal draws from SEs)
  int_draws <- rnorm(n, mean =  0.056197, sd = 0.014901)
  btl_draws <- rnorm(n, mean = -0.012621, sd = 0.001159)

  mu_slope <- int_draws + btl_draws * TL   # expected daily BV change for this TL

  resid_draws <- rnorm(n, mean = 0, sd = 0.008191)  # Model 2 residual: individual variation not explained by length

  mu_slope + resid_draws   # m3/day (negative = losing volume)
}

predict_mother_energy_loss <- function(
    mother_lengths,
    rho_b = 986,   # blubber density kg/m3 (van Aswegen post-mortem calves, 986 +- 24.3)
    E_es  = 26990,   # kJ/kg - largest Hawaii humpback calf (van Aswegen et al., age-specific growth costs)
    n     = 5000
) {
  purrr::map_dfr(mother_lengths, function(L) {
    
    bv_loss_draws      <- predict_daily_bv_loss(L, n)
    daily_energy_draws <- (abs(bv_loss_draws) * rho_b * E_es) / 1000
    
    data.frame(
      MotherLength        = L,
      energy_loss_mean_MJ = mean(daily_energy_draws),
      energy_loss_lwr_MJ  = quantile(daily_energy_draws, 0.025),
      energy_loss_upr_MJ  = quantile(daily_energy_draws, 0.975)
    )
  })
}

predict_mother_energy_loss_v2 <- function(
    mother_lengths,
    E_blubber_draws,   # pre-drawn, shared across lengths
    rho_b = 986,   # blubber density kg/m3 (van Aswegen post-mortem calves, 986 +- 24.3)
    n     = 5000
) {
  purrr::map_dfr(mother_lengths, function(L) {
    bv_loss_draws      <- predict_daily_bv_loss(L, n)
    daily_energy_draws <- (abs(bv_loss_draws) * rho_b * E_blubber_draws) / 1000
    data.frame(
      MotherLength        = L,
      energy_loss_mean_MJ = mean(daily_energy_draws),
      energy_loss_lwr_MJ  = quantile(daily_energy_draws, 0.025),
      energy_loss_upr_MJ  = quantile(daily_energy_draws, 0.975)
    )
  })
}

smooth_vec <- function(x, y, span = 0.5) predict(loess(y ~ x, span = span))

to_left  <- function(x) x * scale_factor + offset

to_right <- function(y) (y - offset) / scale_factor


# ================= from 04_growth_departure_reserves =================

age_at_threshold <- function(ML, threshold_length = 6.0) {
  RL_target <- threshold_length / ML
  # check if threshold is reachable (RL_target must be < asymptote)
  if (RL_target >= 0.737106) return(NA_real_)
  age <- -log((0.737106 - RL_target) / (0.737106 - 0.3504304)) / 0.006740519
  return(age)
}

age_at_threshold_lwr <- function(ML, threshold_length = 6.0) {
  RL_target <- threshold_length / ML
  if (RL_target >= 0.7250429) return(NA_real_)
  age <- -log((0.7250429 - RL_target) / (0.7250429 - 0.3471785)) / 0.006297628
  return(age)
}

age_at_threshold_upr <- function(ML, threshold_length = 6.0) {
  RL_target <- threshold_length / ML
  if (RL_target >= 0.7503627) return(NA_real_)
  age <- -log((0.7503627 - RL_target) / (0.7503627 - 0.3536075)) / 0.007190727
  return(age)
}


# ================= from 05_figures_tables =================

s2p <- function(x) {
  (x - secondary_range[1]) / diff(secondary_range) * diff(primary_range) + primary_range[1]
}

p2s <- function(x) {
  (x - primary_range[1]) / diff(primary_range) * diff(secondary_range) + secondary_range[1]
}

fmt_ci <- function(mn, lo, hi, d = 1) {
  sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), mn, lo, hi)
}


# ================= from 06_sensitivity =================

simulate_daily_energy_SA <- function(row, p_rest, n = 200,
                                     yatesCT, yatesUProp,
                                     mom_propstroke,
                                     u_met = 0.25, alpha_deg = 20,
                                     amp_frac = 0.10, water_density = 1025,
                                     bmr_mult = 1.0) {

  mom_BL  <- as.numeric(row$MotherLength)

  n_rest   <- 1440 * p_rest
  n_active <- 1440 - n_rest

  # morphology
  mom_chord  <- pmax(as.vector(predict_chord(mom_BL, n)),  0.05)
  mom_area   <- pmax(as.vector(predict_area(mom_BL, n)),   0.05)

  # locomotion kinematics
  U_r_mom  <- predict_speed("mother", "rest",  n)
  U_a_mom  <- predict_speed("mother", "other", n)

  f_r_mom  <- pmax(as.vector(predict_frequency(mom_BL,  n)), 0.01)
  f_a_mom  <- f_r_mom

  # metabolic power with the sensitivity parameters threaded through
  mp <- function(U, f, BL, chord, area)
    calc_metabolic_power_full(U, f, BL, chord, area, yatesCT, yatesUProp, n,
                              u_met = u_met, water_density = water_density,
                              alpha_deg = alpha_deg, amp_frac = amp_frac)

  P_r_mom  <- mp(U_r_mom,  f_r_mom,  mom_BL,  mom_chord,  mom_area)
  P_a_mom  <- mp(U_a_mom,  f_a_mom,  mom_BL,  mom_chord,  mom_area)

  prop_fluke_r_mom  <- predict_prop_fluking(mom_propstroke, "rest",  n)
  prop_fluke_a_mom  <- predict_prop_fluking(mom_propstroke, "other", n)

  loco_mom  <- pmax(((n_rest   * 60 * prop_fluke_r_mom * P_r_mom) +
                     (n_active * 60 * prop_fluke_a_mom * P_a_mom)) / 1e6, 0)

  # BMR (Kleiber point estimate x bmr_mult; baseline bmr_mult = 1, i.e. fixed)
  bmr_mom   <- 0.04843 * 70 * (row$MotherMass_kg ^ 0.75) * bmr_mult * 86400 / 1e6

  total_mom  <- loco_mom + bmr_mom

  data.frame(
    MotherLength    = mom_BL,
    Age             = row$Age,
    p_rest          = p_rest,
    mom_total_mean  = mean(total_mom)
  )
}

run_sim_grid_SA <- function(params, focal_lengths = seq(10, 15, 0.5),
                            focal_age = 50, p_rest = 0.4, n = 200) {
  rows <- full_table %>%
    dplyr::mutate(MLr = round(MotherLength, 1)) %>%
    dplyr::filter(MLr %in% round(focal_lengths, 1), Age == focal_age) %>%
    dplyr::distinct(MLr, .keep_all = TRUE) %>%
    dplyr::arrange(MLr)

  purrr::map_dfr(seq_len(nrow(rows)), function(i)
    simulate_daily_energy_SA(
      rows[i, ], p_rest = p_rest, n = n,
      yatesCT = YatesCT, yatesUProp = YatesUProp,
      mom_propstroke = mom_propstroke,
      u_met         = params$u_met,
      alpha_deg     = params$alpha_deg,
      amp_frac      = params$amp_frac,
      water_density = params$water_density,
      bmr_mult      = params$bmr_mult
    ))
}

compute_pattern_SA <- function(sim_grid,
                               E_es = 26990, rho_blubber = 986,
                               focal_lengths = seq(10, 15, 0.5)) {
  L <- focal_lengths

  bvloss      <- abs(0.056197 - 0.012621 * L)
  E_blub_loss <- bvloss * rho_blubber * E_es / 1000

  mom_DEE  <- approx(sim_grid$MotherLength, sim_grid$mom_total_mean,  xout = L, rule = 2)$y

  # absolute lactation supply (MJ/day)
  lactation        <- E_blub_loss - mom_DEE

  # mass-specific lactation cost (MJ/day/kg) - paper's Figure 1 red curve
  mother_mass      <- predict_mother_mass(L)
  lactation_per_kg <- lactation / mother_mass

  # NEW: % of total daily energy loss going to calf - paper's Figure 1 blue curve
  pct_to_calf      <- lactation / E_blub_loss * 100

  # reserves remaining (kept for the table but not the figure - see text)
  expected_BV     <- exp(-3.891 + 2.845 * log(L))
  body_mass       <- exp( 3.069 + 2.825 * log(L))
  # structural tissue MASS -> VOLUME via tissue-specific densities (7 post-mortem calves; van Aswegen 2025b)
  # bone density uncertain: 7-calf mean 1110 (range 963-1257); toggle for sensitivity
  musc_vol        <- (body_mass * exp(1.934 + 0.657 * log(L)) / 100) / 971
  bone_vol        <- (body_mass * exp(4.441 - 0.705 * log(L)) / 100) / 1110
  visc_vol        <- (body_mass * exp(3.276 - 0.369 * log(L)) / 100) / 954
  blood_vol       <- expected_BV * 0.06                    # body fluids, 6% of BV (Lockyer 1976)
  blubb_vol_birth <- expected_BV - musc_vol - bone_vol - visc_vol - blood_vol
  tank_MJ         <- blubb_vol_birth * rho_blubber * E_es / 1000
  dep_age         <- approx(departure_df$MotherLength, departure_df$age_mean,
                            xout = L, rule = 2)$y
  expended        <- E_blub_loss * dep_age
  reserves_pct    <- (tank_MJ - expended) / tank_MJ * 100

  data.frame(
    MotherLength     = L,
    lactation        = lactation,
    lactation_per_kg = lactation_per_kg,
    pct_to_calf      = pct_to_calf,           # NEW
    reserves_pct     = reserves_pct
  )
}

pattern_metrics <- function(df, outcome) {
  y <- df[[outcome]]; x <- df$MotherLength
  fit <- lm(y ~ x)
  data.frame(
    slope       = unname(coef(fit)[2]),
    level_at_12 = unname(predict(fit, newdata = data.frame(x = 12))),
    spearman    = suppressWarnings(cor(x, y, method = "spearman"))
  )
}

eval_scenario <- function(mods = list()) {
  p     <- modifyList(baseline_params, mods)
  rerun <- any(vapply(sim_param_names,
                      function(nm) !isTRUE(all.equal(p[[nm]], baseline_params[[nm]])),
                      logical(1)))
  grid  <- if (rerun) run_sim_grid_SA(p, focal_lengths = focal_lengths) else baseline_grid
  compute_pattern_SA(grid, E_es = p$E_es, rho_blubber = p$rho_blubber,
                     focal_lengths = focal_lengths)
}

add_scenario <- function(label, parameter, mods) {
  pat <- eval_scenario(mods)
  # store the full pattern curve for the figure
  sa_patterns[[label]] <<- pat %>%
    dplyr::mutate(scenario = label, parameter = parameter)

  # store the summary metrics for the table
  for (oc in outcomes) {
    m <- pattern_metrics(pat, oc)
    sa_rows[[length(sa_rows) + 1]] <<- data.frame(
      scenario = label, parameter = parameter, outcome = oc,
      slope = m$slope, level_at_12 = m$level_at_12, spearman = m$spearman
    )
  }
  invisible(pat)
}

