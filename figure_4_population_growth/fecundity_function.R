# ============================================================
# humpback_fecundity_function.R
#
# Provides a function to predict calving probability and
# interbirth interval (IBI) for female humpback whales as a
# function of body length.
#
# USAGE
#   source("humpback_fecundity_function.R")
#
#   # Point estimate at a single length
#   humpback_fecundity(L = 13.0)
#
#   # Point estimates at a vector of lengths
#   humpback_fecundity(L = c(11, 12, 13, 14))
#
#   # Monte Carlo draws (returns a list of matrices)
#   humpback_fecundity(L = 13.0, n_draws = 1000)
#
# OUTPUT
#   A data frame (point estimate) or list of data frames
#   (Monte Carlo), each with columns:
#     L           : input body length (m)
#     P_calving   : calving probability asymptote
#     IBI         : interbirth interval (years) = 1 / P_calving
#
# MODEL
#   Adapts the logistic L^3 model of Pirotta et al. (2024),
#   which relates female body length to the asymptote of
#   calving probability (i.e. the maximum calving probability
#   achievable at a given length, at peak health). The slope
#   and predictor standardization are taken directly from
#   Pirotta et al.; only the intercept is re-anchored to
#   Hawaii humpback calving rates.
#
#   logit(m) = lambda + zeta * (L^3 - mu_L3) / sigma_L3
#
#   where:
#     m        = calving probability asymptote
#     lambda   = intercept (population-specific; see ANCHOR below)
#     zeta     = 0.75 [0.52, 1.00] (Pirotta et al. 2024, 95% CrI)
#     mu_L3    = 1987.2 m^3  (mean L^3 in Pirotta fitting dataset)
#     sigma_L3 = 515.5  m^3  (SD   L^3 in Pirotta fitting dataset)
#
#   IBI = 1 / m  (raw calving probability, not conditioned on
#   availability; appropriate because Hawaii humpbacks include
#   documented annual calvings, so minimum IBI may be 1 yr).
#
# ANCHOR
#   The intercept lambda is solved so that m(L_anchor) = P_anchor.
#   Two anchor scenarios are available via the `anchor` argument:
#
#   "herman"  Herman et al. (2011), Hawaii, n = 39 females
#             P_anchor = 0.48; preferred anchor (Hawaii-specific,
#             larger sample, more recent)
#
#   "baker"   Baker et al. (1987), Central/Eastern North Pacific,
#             n = 18 females
#             P_anchor = 0.58; wider uncertainty; not exclusively
#             Hawaii residents - use for sensitivity analysis only
#
#   Anchor length: 13.0 m (mean adult female, North Pacific
#   humpback; van Aswegen et al. (2025) J Physiol 603:551-578,).
#   Update ANCHOR_LENGTH below if a better estimate is available.
#
# UNCERTAINTY (Monte Carlo mode, n_draws > 0)
#   Samples jointly from:
#   (1) zeta  ~ truncated Normal(0.75, sd), where sd is derived
#       from the Pirotta 95% CrI [0.52, 1.00]; truncated at 0
#       to enforce positive effect of length on calving probability
#   (2) P_anchor ~ logit-Normal centred on the reported calving
#       rate, with SD on logit scale derived from binomial
#       uncertainty (Wilson CI) given the study sample size:
#         herman: SD_logit ~ 0.306  (n = 39)
#         baker:  SD_logit ~ 0.435  (n = 18)
#   These two sources of uncertainty are treated as independent.
#   The anchor length is held fixed (uncertainty is second-order).
#
# REFERENCES
#   Pirotta et al. (2024) R. Soc. Open Sci. 11:240050
#     https://doi.org/10.1098/rsos.240050
#   Herman et al. (2011) Mar. Mamm. Sci.
#   Baker et al. (1987)  Mar. Ecol. Prog. Ser. 41:103-114
#   van Aswegen et al. (2025) J Physiol 603:551-578, doi:10.1113/JP287379
# ============================================================

# ---- Fixed parameters from Pirotta et al. (2024) ----------------

.PIROTTA <- list(
  zeta_mean  =  0.75,    # posterior mean, slope on standardized L^3
  zeta_lo    =  0.52,    # 95% CrI lower
  zeta_hi    =  1.00,    # 95% CrI upper
  mu_L3      =  1987.2,  # mean(L^3) in right whale fitting dataset (m^3)
  sigma_L3   =  515.5    # SD(L^3)   in right whale fitting dataset (m^3)
)

# ---- Anchor scenarios -------------------------------------------

.ANCHORS <- list(
  herman = list(
    P       = 0.48,    # Herman et al. (2011) Hawaii calving rate
    n       = 39,      # number of females
    sd_logit = 0.306,  # approx SD on logit scale (from binomial Wilson CI)
    cite    = "Herman et al. (2011) Mar. Mamm. Sci."
  ),
  baker = list(
    P       = 0.58,    # Baker et al. (1987) calving rate
    n       = 18,
    sd_logit = 0.435,  # wider - smaller sample, mixed population
    cite    = "Baker et al. (1987) Mar. Ecol. Prog. Ser. 41:103-114"
  )
)

# ---- Anchor length ----------------------------------------------
# Mean adult female body length, North Pacific humpback whales.
# Source: van Aswegen et al. (2025) J Physiol 603:551-578, doi:10.1113/JP287379.
# Update this value if a better population-specific estimate
# becomes available.

ANCHOR_LENGTH <- 13.0   # metres

# ---- Internal helpers -------------------------------------------

.solve_lambda <- function(P_anchor, L_anchor,
                           zeta, mu_L3, sigma_L3) {
  qlogis(P_anchor) - zeta * (L_anchor^3 - mu_L3) / sigma_L3
}

.predict_m <- function(L, lambda, zeta, mu_L3, sigma_L3) {
  plogis(lambda + zeta * (L^3 - mu_L3) / sigma_L3)
}

# ---- Main function ----------------------------------------------

#' Predict humpback whale calving probability and IBI from body length
#'
#' @param L         Numeric vector of female body lengths (metres)
#' @param anchor    Character: "herman" (default) or "baker"
#' @param n_draws   Integer: number of Monte Carlo draws (0 = point estimate)
#' @param seed      Random seed for reproducibility (used when n_draws > 0)
#'
#' @return
#'   If n_draws == 0: a data frame with columns L, P_calving, IBI
#'   If n_draws  > 0: a list with elements
#'     $summary  - data frame with columns L, P_mean, P_lo, P_hi,
#'                 IBI_mean, IBI_lo, IBI_hi  (95% uncertainty interval)
#'     $draws    - matrix [n_draws x length(L)] of P_calving samples
#'     $IBI_draws - matrix [n_draws x length(L)] of IBI samples

humpback_fecundity <- function(L,
                                anchor  = "herman",
                                n_draws = 0,
                                seed    = 42) {

  # Validate inputs
  anchor <- match.arg(anchor, c("herman", "baker"))
  stopifnot(is.numeric(L), all(L > 0))
  stopifnot(is.numeric(n_draws), n_draws >= 0)

  anc  <- .ANCHORS[[anchor]]
  pir  <- .PIROTTA

  # ---- Point estimate -------------------------------------------
  if (n_draws == 0) {
    lambda <- .solve_lambda(anc$P, ANCHOR_LENGTH,
                             pir$zeta_mean, pir$mu_L3, pir$sigma_L3)
    P      <- .predict_m(L, lambda, pir$zeta_mean, pir$mu_L3, pir$sigma_L3)
    return(data.frame(L = L, P_calving = P, IBI = 1 / P))
  }

  # ---- Monte Carlo ----------------------------------------------
  set.seed(seed)

  # Approximate SD for zeta from 95% CrI (assume Normal)
  zeta_sd <- (pir$zeta_hi - pir$zeta_lo) / (2 * 1.96)

  # Draw zeta > 0 (truncated Normal; rejection sample - fast given
  # that P(zeta <= 0) is negligibly small here)
  zeta_draws <- rnorm(n_draws, pir$zeta_mean, zeta_sd)
  zeta_draws <- pmax(zeta_draws, 1e-6)   # enforce positive

  # Draw P_anchor on logit scale, then back-transform
  P_logit_draws <- rnorm(n_draws, qlogis(anc$P), anc$sd_logit)
  P_draws       <- plogis(P_logit_draws)

  # For each draw, solve lambda and predict P at each length
  P_mat <- matrix(NA_real_, nrow = n_draws, ncol = length(L))

  for (k in seq_len(n_draws)) {
    lambda_k    <- .solve_lambda(P_draws[k], ANCHOR_LENGTH,
                                  zeta_draws[k], pir$mu_L3, pir$sigma_L3)
    P_mat[k, ]  <- .predict_m(L, lambda_k, zeta_draws[k],
                                pir$mu_L3, pir$sigma_L3)
  }

  IBI_mat <- 1 / P_mat

  # Summary across draws
  summary_df <- data.frame(
    L        = L,
    P_mean   = colMeans(P_mat),
    P_lo     = apply(P_mat,   2, quantile, 0.025),
    P_hi     = apply(P_mat,   2, quantile, 0.975),
    IBI_mean = colMeans(IBI_mat),
    IBI_lo   = apply(IBI_mat, 2, quantile, 0.025),
    IBI_hi   = apply(IBI_mat, 2, quantile, 0.975)
  )

  list(summary   = summary_df,
       draws     = P_mat,
       IBI_draws = IBI_mat)
}

# ---- Convenience: print anchor info -----------------------------

fecundity_anchor_info <- function(anchor = "herman") {
  anchor <- match.arg(anchor, c("herman", "baker"))
  anc    <- .ANCHORS[[anchor]]
  cat(sprintf("Anchor scenario: %s\n", anchor))
  cat(sprintf("  Source:        %s\n", anc$cite))
  cat(sprintf("  P_anchor:      %.2f  (calving rate)\n", anc$P))
  cat(sprintf("  n females:     %d\n", anc$n))
  cat(sprintf("  SD(logit P):   %.3f  (binomial uncertainty)\n", anc$sd_logit))
  cat(sprintf("  Anchor length: %.1f m  (van Aswegen et al.)\n", ANCHOR_LENGTH))
  cat(sprintf("  zeta:          %.2f [%.2f, %.2f]  (Pirotta et al. 2024)\n",
      .PIROTTA$zeta_mean, .PIROTTA$zeta_lo, .PIROTTA$zeta_hi))
}
