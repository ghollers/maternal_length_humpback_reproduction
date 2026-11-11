# ============================================================
# Figure 4 - Leslie-matrix population growth rate (lambda) vs female
#            asymptotic length, under length-decline scenarios.
#
# Self-contained: sources fecundity_function.R (in this folder) and uses
# published/hardcoded demographic parameters. No raw data needed.
#
# Run: source("make_figure.R")  ->  figure_4_population_growth.png
#
# (original script: leslie_population_growth.R)
#
# Intrinsic population growth rate (lambda) of Hawaiian humpback
# whales as a function of FEMALE ASYMPTOTIC LENGTH (L_inf).
#
# lambda = the dominant eigenvalue of an age-structured Leslie
# matrix whose fertility schedule is driven by a length -> calving
# probability relationship (Pirotta et al. 2024, re-anchored to
# Hawaii calving rates). Uncertainty in the demographic parameters
# (survival, maximum age, age at first reproduction) is propagated
# by Monte Carlo to give a DISTRIBUTION of lambda for each
# asymptotic-length scenario.
#
# ------------------------------------------------------------
# TERMINOLOGY (to avoid the earlier confusion):
#   lambda        = INTRINSIC POPULATION GROWTH RATE (this script).
#                   lambda > 1 growing, = 1 stable, < 1 declining.
#   The object called "lambda" INSIDE humpback_fecundity_function.R
#   is a DIFFERENT quantity: the intercept of the Pirotta calving-
#   probability model. It is internal to that function and never
#   used here. (Recommend renaming it there, e.g. to `beta0`.)
# ------------------------------------------------------------
#
# Requires: humpback_fecundity_function.R in the SAME folder.
# Run from that folder (setwd or open the .Rproj first).
# ============================================================

library(truncnorm)   # truncated-normal draws for parameters
library(ggplot2)
library(dplyr)
library(patchwork)
# Locate the helper regardless of working directory (works from the package root
# when opened via humpback_maternal_energetics.Rproj, or from the R/ folder directly).
source("fecundity_function.R")   # provides humpback_fecundity(L, anchor)

out_dir <- "."   # write the figure next to this script


# ============================================================
# 1. CORE DEMOGRAPHIC FUNCTIONS
# ============================================================

# --- age -> length: female von Bertalanffy growth --------------
# L_0 and k are fixed to the project's female values so the whole
# analysis stays consistent with the growth model used elsewhere.
vb_length_at_age <- function(age, L_inf, L_0 = 8.42, k = 0.15) {
  L_inf - (L_inf - L_0) * exp(-k * age)
}

# --- age -> calving probability (via length) -------------------
# Immature females (age < A_r) do not calve.
calving_prob_at_age <- function(age, L_inf, anchor = "herman",
                                L_0 = 8.42, k = 0.15, A_r = 6) {
  L <- vb_length_at_age(age, L_inf, L_0, k)
  P <- humpback_fecundity(L, anchor = anchor)$P_calving
  P[age < A_r] <- 0
  P
}

# --- intrinsic population growth rate = dominant eigenvalue -----
lambda_dom <- function(M) max(Re(eigen(M, only.values = TRUE)$values))

# --- build ONE Leslie matrix and return its lambda -------------
# Survival schedule (matches Alec's construction):
#   P_sc first-year, P_sj juvenile (ages 2-9), P_sa adult (10+).
leslie_lambda <- function(L_inf, anchor = "herman",
                          omega = 70, A_r = 6,
                          P_sc = 0.85, P_sj = 0.93, P_sa = 0.97,
                          sex_ratio = 0.5,
                          L_0 = 8.42, k = 0.15) {

  # age-specific survival on the sub-diagonal (age i -> i+1)
  surv <- numeric(omega - 1)
  surv[1]              <- P_sc
  surv[2:9]            <- P_sj
  surv[10:(omega - 1)] <- P_sa

  survival  <- diag(surv, omega - 1, omega - 1)
  survival2 <- cbind(survival, 0)                # add terminal-age column

  # length-driven calving probability at each age (the fertility driver)
  br <- calving_prob_at_age(1:(omega - 1), L_inf, anchor, L_0, k, A_r)

  # --- fertility (top) row of the Leslie matrix ---
  # NOTE / VERIFY WITH ALEC: this reproduces his construction, which
  # multiplies calving probability by the *age-specific* survival
  # (surv) and then the female offspring fraction. Depending on the
  # census timing (pre- vs post-breeding) you may instead want to
  # multiply by first-year calf survival (P_sc) only. Confirm before
  # finalizing; it shifts the level of lambda but not its ranking
  # across length scenarios.
  f <- c(br * surv, 0) * sex_ratio

  M <- rbind(f, survival2)                        # omega x omega matrix
  lambda_dom(M)
}


# ============================================================
# 2. PARAMETER RANGES   <-- EDIT THESE (placeholders)
# ============================================================
# Each demographic parameter is drawn from a plausible distribution
# so its uncertainty propagates into lambda. Values below are
# literature-grounded (see comments); survival and A_r use truncated
# normals, omega a uniform range.
#
# LITERATURE BENCHMARKS (humpback whale demography):
#   Calf (first-year) survival : 0.875 (SE 0.047)   Barlow & Clapham (1997), Gulf of Maine
#   Non-calf survival          : 0.960 (SE 0.008)   Barlow & Clapham (1997)
#                                range 0.95-0.985 across ocean basins; females
#                                higher (0.98-0.99); Zerbini et al. (2010) used ~0.985
#                                as an upper bound.
#   Age at first parturition   : youngest 5 yr, mean ~6 yr (Barlow & Clapham 1997;
#                                Clapham 1992). N. Pacific may run higher - confirm for Hawaii.
#   Maximum age (omega)        : oldest confirmed ~95 yr; ~90+ typical. (Little effect on lambda.)
PARAMS <- list(
  # first-year calf survival - Barlow & Clapham (1997): 0.875 (SE 0.047)
  P_sc = list(mean = 0.875, sd = 0.047, lo = 0.75, hi = 0.95),
  # juvenile survival - no humpback-specific juvenile/adult split exists;
  # set just below adult (literature lumps these as "non-calf" survival)
  P_sj = list(mean = 0.95,  sd = 0.02,  lo = 0.90, hi = 0.98),
  # adult (non-calf) survival - B&C 0.960 (SE 0.008); across studies 0.95-0.985,
  # females at the high end. Center 0.96, allow up to ~0.99.
  P_sa = list(mean = 0.96,  sd = 0.015, lo = 0.94, hi = 0.99),
  # age at first reproduction - youngest 5 yr, mean ~6 (B&C 1997; Clapham 1992).
  # NOTE: North Pacific / Hawaii estimates may differ; treated as a sensitivity parameter.
  A_r  = list(mean = 6.0,   sd = 1.0,   lo = 5,    hi = 9),
  # maximum age - oldest confirmed ~95 yr (earplug); negligible effect on lambda
  omega_lo = 85, omega_hi = 95,
  sex_ratio = 0.5                                             # female offspring fraction
)

# draw n independent demographic parameter sets
draw_params <- function(n) {
  data.frame(
    P_sc      = rtruncnorm(n, PARAMS$P_sc$lo, PARAMS$P_sc$hi, PARAMS$P_sc$mean, PARAMS$P_sc$sd),
    P_sj      = rtruncnorm(n, PARAMS$P_sj$lo, PARAMS$P_sj$hi, PARAMS$P_sj$mean, PARAMS$P_sj$sd),
    P_sa      = rtruncnorm(n, PARAMS$P_sa$lo, PARAMS$P_sa$hi, PARAMS$P_sa$mean, PARAMS$P_sa$sd),
    A_r       = round(rtruncnorm(n, PARAMS$A_r$lo, PARAMS$A_r$hi, PARAMS$A_r$mean, PARAMS$A_r$sd)),
    omega     = round(runif(n, PARAMS$omega_lo, PARAMS$omega_hi)),
    sex_ratio = PARAMS$sex_ratio
  )
}


# ============================================================
# 3. MONTE CARLO: distribution of lambda for a given L_inf
# ============================================================
# Hold asymptotic length fixed (it is the treatment) and draw the
# demographic parameters n_sim times, computing lambda for each.
simulate_lambda <- function(L_inf, anchor = "herman", n_sim = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- draw_params(n_sim)
  vapply(seq_len(n_sim), function(i)
    leslie_lambda(L_inf, anchor,
                  omega = p$omega[i], A_r = p$A_r[i],
                  P_sc = p$P_sc[i], P_sj = p$P_sj[i], P_sa = p$P_sa[i],
                  sex_ratio = p$sex_ratio[i]),
    numeric(1))
}

# mean and 95% interval of a lambda sample
summarise_lambda <- function(lam) {
  c(mean = mean(lam),
    lwr  = unname(quantile(lam, 0.025)),
    upr  = unname(quantile(lam, 0.975)))
}


# ============================================================
# 4. FOUR MANUSCRIPT SCENARIOS (asymptotic-length gradient)
# ============================================================
scenarios <- data.frame(
  name  = c("Historical", "Current", "Moderate decline", "Severe decline"),
  L_inf = c(13.7,          13.0,      12.5,               12.0)
)

set.seed(1)
scenario_lambda <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  L  <- scenarios$L_inf[i]
  sh <- summarise_lambda(simulate_lambda(L, "herman", n_sim = 2000))  # primary anchor
  sb <- summarise_lambda(simulate_lambda(L, "baker",  n_sim = 2000))  # sensitivity anchor
  data.frame(
    scenario   = scenarios$name[i],
    L_inf      = L,
    lambda_herman     = sh["mean"],
    lambda_herman_lwr = sh["lwr"],
    lambda_herman_upr = sh["upr"],
    lambda_baker      = sb["mean"],
    lambda_baker_lwr  = sb["lwr"],
    lambda_baker_upr  = sb["upr"]
  )
}))
row.names(scenario_lambda) <- NULL
print(scenario_lambda, digits = 4)
# scenario lambda values are also provided in ../data/processed if needed
write.csv(scenario_lambda, file.path(out_dir, "figure_4_scenarios_lambda.csv"), row.names = FALSE)

# --- HEADLINE: RELATIVE effect of asymptotic length on growth potential ---
# We do NOT interpret the absolute lambda as the realized growth rate:
# the demographic parameters are from earlier recovery-phase studies and
# the North Pacific / Hawaii population has since approached carrying
# capacity and declined (Cheeseman et al. 2024). Because all demographic
# parameters are held identical across scenarios, DIFFERENCES in lambda
# reflect the effect of asymptotic length alone. We express each scenario
# as the % change in the per-capita growth increment r = lambda - 1
# relative to the current-length baseline.
scenario_lambda$r              <- scenario_lambda$lambda_herman - 1
r_cur                          <- scenario_lambda$r[scenario_lambda$scenario == "Current"]
scenario_lambda$pct_r_vs_current <- (scenario_lambda$r / r_cur - 1) * 100
print(scenario_lambda[, c("scenario", "L_inf", "lambda_herman", "r", "pct_r_vs_current")],
      digits = 4)

# structural sanity check only (not a prediction): lambda should fall
# within the biologically plausible envelope for the species
# (<= ~1.118; Zerbini et al. 2010 maximum plausible rate of increase).
cat(sprintf("\nAll scenario lambda within plausible envelope (<= 1.118): %s\n",
            all(scenario_lambda$lambda_herman <= 1.118)))


# ============================================================
# 5. CONTINUOUS SWEEP + FIGURE (lambda vs asymptotic length)
# ============================================================
L_seq <- seq(12.0, 14.0, by = 0.1)

set.seed(2)
sweep <- do.call(rbind, lapply(L_seq, function(L) {
  s <- summarise_lambda(simulate_lambda(L, "herman", n_sim = 1000))
  data.frame(L_inf = L, lambda_mean = s["mean"], lambda_lwr = s["lwr"], lambda_upr = s["upr"])
}))
row.names(sweep) <- NULL

fig_lambda <- ggplot(sweep, aes(L_inf, lambda_mean)) +
  # Monte Carlo 95% band
  geom_ribbon(aes(ymin = lambda_lwr, ymax = lambda_upr),
              fill = "#2c7bb6", alpha = 0.20) +
  geom_line(color = "#2c7bb6", linewidth = 1.1) +
  # replacement level (stable population)
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  # the four named scenarios with their CIs
  geom_errorbar(data = scenario_lambda,
                aes(x = L_inf, y = lambda_herman,
                    ymin = lambda_herman_lwr, ymax = lambda_herman_upr),
                width = 0.05, color = "#d7191c") +
  geom_point(data = scenario_lambda,
             aes(x = L_inf, y = lambda_herman),
             color = "#d7191c", size = 2.6) +
  scale_x_continuous("Female asymptotic length (m)", breaks = seq(12, 14, 0.5)) +
  scale_y_continuous(expression("Intrinsic population growth rate ("*lambda*")")) +
  theme_classic(base_size = 13)

fig_lambda
# (extra exploratory figure not part of the manuscript; omitted)


# ============================================================
# 6. (OPTIONAL) sensitivity of lambda to each parameter
# ============================================================
# Which demographic parameter drives most of the spread in lambda?
# Correlate lambda against each drawn parameter at a fixed length.
# set.seed(3)
# p  <- draw_params(2000)
# lam <- vapply(seq_len(nrow(p)), function(i)
#   leslie_lambda(13.0, "herman", p$omega[i], p$A_r[i],
#                 p$P_sc[i], p$P_sj[i], p$P_sa[i], p$sex_ratio[i]), numeric(1))
# sapply(p[c("P_sc","P_sj","P_sa","A_r","omega")],
#        function(x) cor(x, lam))


# ============================================================
# 7. RELATIVE-EFFECT FIGURES (candidate manuscript figures)
# ============================================================
# The absolute lambda values are not predictions of realized growth, so the
# manuscript figure should report the RELATIVE effect of asymptotic length.
# Both options below put a RELATIVE metric on the y-axis, so absolute lambda
# never appears in the figure.
#
# KEY DESIGN POINT - PAIRED DRAWS:
# Sections 4-5 draw INDEPENDENT demographic parameter sets for each length, so
# a ratio between two lengths would inherit both scenarios' full uncertainty
# and the relative interval would look needlessly wide. Below we instead draw
# ONE parameter set and evaluate EVERY length with that SAME set, then take
# ratios within each draw. That isolates the effect of length and gives an
# honest (and much tighter) interval on the RELATIVE difference.

library(patchwork)   # for combining panels

simulate_relative_lambda <- function(lengths, baseline_L = 13.0,
                                     anchor = "herman", n_sim = 2000, seed = 42) {
  set.seed(seed)
  p <- draw_params(n_sim)        # ONE set of draws, reused for every length

  # lambda at every length under each parameter draw -> [n_lengths x n_sim]
  lam <- vapply(seq_len(n_sim), function(i)
    vapply(lengths, function(L)
      leslie_lambda(L, anchor,
                    omega = p$omega[i], A_r = p$A_r[i],
                    P_sc = p$P_sc[i], P_sj = p$P_sj[i], P_sa = p$P_sa[i],
                    sex_ratio = p$sex_ratio[i]),
      numeric(1)),
    numeric(length(lengths)))

  r        <- lam - 1                            # per-capita growth increment
  base_idx <- which.min(abs(lengths - baseline_L))
  r_base   <- r[base_idx, ]                      # baseline increment, per draw
  pct      <- t(t(r) / r_base) * 100 - 100       # % change vs baseline, per draw

  data.frame(
    L_inf    = lengths,
    pct_mean = rowMeans(pct),
    pct_lwr  = apply(pct, 1, quantile, 0.025),
    pct_upr  = apply(pct, 1, quantile, 0.975)
  )
}

# NOTE: ~2000 draws x 21 lengths can take ~1 min; lower n_sim while exploring.
rel_df <- simulate_relative_lambda(L_seq, baseline_L = 13.0, n_sim = 2000)

# scenario points read off the same curve
rel_scen <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  j <- which.min(abs(rel_df$L_inf - scenarios$L_inf[i]))
  data.frame(scenario = scenarios$name[i], rel_df[j, ])
}))
row.names(rel_scen) <- NULL
print(rel_scen, digits = 3)

# Colour-blind-safe palette (Okabe-Ito). Ordered cool -> warm so larger
# (historical) lengths read cool and declining lengths read warm. Avoids the
# red/green pair, which is indistinguishable under deuteranopia/protanopia.
scen_cols <- c("Historical"       = "#0072B2",  # blue
               "Current"          = "#56B4E9",  # sky blue
               "Moderate decline" = "#E69F00",  # orange
               "Severe decline"   = "#D55E00")  # vermillion

# Redundant linetype encoding: the panel stays readable in greyscale and for
# any form of colour vision deficiency (never rely on hue alone).
scen_lty <- c("Historical"       = "longdash",
              "Current"          = "solid",
              "Moderate decline" = "dashed",
              "Severe decline"   = "dotted")


# ------------------------------------------------------------
# OPTION A - two panels: mechanism (a) + relative consequence (b)
# ------------------------------------------------------------

# (a) the mechanism: the whole fecundity schedule shifts down with length.
# NOTE: plotted from the baseline age at first reproduction (A_r = 6) rather
# than age 1 - below A_r calving probability is structurally zero, so starting
# at 1 produced a flat zero segment and an artificial vertical jump at age 6.
# (A_r is drawn from 5-9 in the Monte Carlo; panel (a) shows the baseline
# value of 6 - state this in the figure caption.)
A_r_plot  <- 6
ages_plot <- A_r_plot:40
calv_df <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i)
  data.frame(age      = ages_plot,
             scenario = scenarios$name[i],
             P        = calving_prob_at_age(ages_plot, scenarios$L_inf[i], "herman",
                                            A_r = A_r_plot))))
calv_df$scenario <- factor(calv_df$scenario, levels = scenarios$name)

# Direct labels at the right-hand end of each curve instead of a legend:
# cleaner than a detached bottom legend, and it removes colour-matching from
# the reader entirely (the most accessible option). The x-axis is padded to
# 52 to make room for the text, but ticks stop at 40.
lab_df <- calv_df[calv_df$age == max(calv_df$age), ]

p_a <- ggplot(calv_df, aes(age, P, colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 1) +
  geom_text(data = lab_df, aes(x = age + 1, y = P, label = scenario),
            hjust = 0, size = 3.8, show.legend = FALSE) +
  scale_colour_manual(values = scen_cols, guide = "none") +
  scale_linetype_manual(values = scen_lty, guide = "none") +
  scale_x_continuous("Female age (years)",
                     breaks = seq(10, 40, 10), limits = c(A_r_plot, 52)) +
  scale_y_continuous("Annual calving probability", limits = c(0, NA)) +
  theme_classic(base_size = 14)

# (b) the consequence, expressed only in relative terms
p_b <- ggplot(rel_df, aes(L_inf, pct_mean)) +
  geom_ribbon(aes(ymin = pct_lwr, ymax = pct_upr), fill = "grey70", alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) +
  geom_point(data = rel_scen, aes(L_inf, pct_mean, colour = scenario), size = 3) +
  scale_colour_manual(values = scen_cols, guide = "none") +
  scale_x_continuous("Female asymptotic length (m)", breaks = seq(12, 14, 0.5)) +
  scale_y_continuous(
    expression(atop("Change in growth increment ("*lambda*" - 1)",
                    "relative to current length (%)"))) +
  theme_classic(base_size = 14)

fig_rel <- (p_a | p_b) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

fig_rel
ggsave(file.path(out_dir, "figure_4_population_growth.png"), fig_rel,
       width = 10, height = 4.5, units = "in", dpi = 300)


# ------------------------------------------------------------
# OPTION B - simpler: scenario lollipop of the relative effect
# ------------------------------------------------------------
rel_lp <- rel_scen
rel_lp$scenario <- factor(rel_lp$scenario, levels = rev(scenarios$name))

fig_rel_lollipop <- ggplot(rel_lp, aes(pct_mean, scenario, colour = scenario)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_segment(aes(x = 0, xend = pct_mean, yend = scenario), linewidth = 0.9) +
  geom_errorbar(aes(xmin = pct_lwr, xmax = pct_upr),
                orientation = "y", width = 0.12) +
  geom_point(size = 4) +
  scale_colour_manual(values = scen_cols, guide = "none") +
  scale_x_continuous(
    expression("Change in growth increment ("*lambda*" - 1) relative to current length (%)")) +
  ylab(NULL) +
  theme_classic(base_size = 12)

fig_rel_lollipop
if (FALSE) ggsave(file.path(out_dir, "fig_lambda_relative_lollipop.png"), fig_rel_lollipop,
       width = 7, height = 3.8, units = "in", dpi = 300)
