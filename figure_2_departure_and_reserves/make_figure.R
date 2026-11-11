# ============================================================
# Figure 2 - Predicted calf departure age (blue, left axis) and maternal
#            blubber reserves remaining at departure (orange, right axis),
#            vs maternal length.
#
# The blue curve is read from the processed departure ages
#   ../data/processed/departure_df.csv
# The reserves curve is recomputed from deterministic body-volume allometry
# (hardcoded published coefficients; van Aswegen et al. 2025). No raw data needed.
#
# Run: source("make_figure.R")  ->  figure_2_departure_and_reserves.png
# ============================================================

library(dplyr)
library(ggplot2)

find_data <- function(fname) {
  for (p in c(file.path("..", "data", "processed", fname),
              file.path("data", "processed", fname), fname))
    if (file.exists(p)) return(p)
  stop("Could not find ", fname, " in ../data/processed/")
}
departure_df <- read.csv(find_data("departure_df.csv"))   # MotherLength, age_mean/lwr/upr (0.1 m grid)

blue_color <- "#0072B2"; red_color <- "#D55E00"

# --- deterministic age at the 6 m threshold (for the focal integer lengths) --
age_at_threshold <- function(ML, RL_inf, RL0, k, threshold = 6.0) {
  RL_target <- threshold / ML
  if (RL_target >= RL_inf) return(NA_real_)
  -log((RL_inf - RL_target) / (RL_inf - RL0)) / k
}
age_mean <- function(ML) age_at_threshold(ML, 0.737106, 0.3504304, 0.006740519)
age_lwr  <- function(ML) age_at_threshold(ML, 0.7250429, 0.3471785, 0.006297628)
age_upr  <- function(ML) age_at_threshold(ML, 0.7503627, 0.3536075, 0.007190727)

# --- deterministic initial blubber volume (m^3) by maternal length ----------
blubber_vol_fn <- function(L) {
  expected_BV <- exp(-3.891 + 2.845 * log(L))
  body_mass   <- exp( 3.069 + 2.825 * log(L))
  m_musc <- body_mass * (exp(1.934 + 0.657 * log(L)) / 100)
  m_bone <- body_mass * (exp(4.441 - 0.705 * log(L)) / 100)
  m_visc <- body_mass * (exp(3.276 - 0.369 * log(L)) / 100)
  v_struct <- m_musc/971 + m_bone/1110 + m_visc/954
  blood <- expected_BV * 0.06
  expected_BV - v_struct - blood
}

focal <- seq(10, 15, by = 1)
surplus_B <- tibble::tibble(MotherLength = focal) %>%
  mutate(
    age_mean = vapply(MotherLength, age_mean, numeric(1)),
    age_lwr  = vapply(MotherLength, age_lwr,  numeric(1)),
    age_upr  = vapply(MotherLength, age_upr,  numeric(1)),
    blubb_vol_birth = blubber_vol_fn(MotherLength),
    vol_loss = abs(0.056197 - 0.012621 * MotherLength),
    pct_remaining_mean = (blubb_vol_birth - vol_loss * age_mean) / blubb_vol_birth * 100,
    pct_remaining_lwr  = (blubb_vol_birth - vol_loss * age_upr)  / blubb_vol_birth * 100,
    pct_remaining_upr  = (blubb_vol_birth - vol_loss * age_lwr)  / blubb_vol_birth * 100
  )

# --- dual-axis mapping (reserves plotted onto the age scale) ----------------
primary_range   <- c(0, 180)
secondary_range <- range(c(surplus_B$pct_remaining_lwr, surplus_B$pct_remaining_upr), na.rm = TRUE)
s2p <- function(x) (x - secondary_range[1]) / diff(secondary_range) * diff(primary_range) + primary_range[1]
p2s <- function(x) (x - primary_range[1]) / diff(primary_range) * diff(secondary_range) + secondary_range[1]
zero_cross <- approx(x = surplus_B$pct_remaining_mean, y = surplus_B$MotherLength, xout = 0)$y
maui_mean  <- 12.61   # van Aswegen et al. 2025

fig2 <- ggplot(mapping = aes(x = MotherLength)) +
  annotate("rect", xmin = 10, xmax = zero_cross, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.5) +
  annotate("text", x = mean(c(10, zero_cross)), y = 165, label = "Depleted\nreserves",
           color = "grey30", size = 3.2, fontface = "italic", lineheight = 0.9) +
  annotate("segment", x = maui_mean, xend = maui_mean, y = 0, yend = 180,
           linetype = "dotted", color = "grey20", linewidth = 0.7) +
  annotate("text", x = maui_mean + 0.06, y = 6, label = "Mean maternal length (12.6 m)",
           color = "grey20", size = 2.9, hjust = 0, vjust = 0) +
  # reserves remaining (right axis)
  geom_ribbon(data = surplus_B, aes(ymin = s2p(pct_remaining_lwr), ymax = s2p(pct_remaining_upr)),
              fill = red_color, alpha = 0.2) +
  geom_line(data = surplus_B, aes(y = s2p(pct_remaining_mean)), color = red_color, linewidth = 1.2) +
  geom_point(data = surplus_B, aes(y = s2p(pct_remaining_mean)), color = red_color, size = 2.5) +
  # departure age (left axis), from the processed data
  geom_ribbon(data = departure_df, aes(ymin = age_lwr, ymax = age_upr), fill = blue_color, alpha = 0.2) +
  geom_line(data = departure_df, aes(y = age_mean), color = blue_color, linewidth = 1.2) +
  scale_x_continuous(name = "Maternal length (m)", breaks = seq(10, 15, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(name = "Predicted age at calf departure (days)",
                     breaks = seq(0, 180, 20), limits = primary_range,
                     sec.axis = sec_axis(~ p2s(.),
                       name = "Reserves remaining at departure (% of initial)",
                       breaks = pretty(secondary_range, n = 5))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(axis.title.y.left  = element_text(color = blue_color, size = 12),
        axis.text.y.left   = element_text(color = blue_color, size = 10),
        axis.line.y.left   = element_line(color = blue_color),
        axis.ticks.y.left  = element_line(color = blue_color),
        axis.title.y.right = element_text(color = red_color, size = 12),
        axis.text.y.right  = element_text(color = red_color, size = 10),
        axis.line.y.right  = element_line(color = red_color),
        axis.ticks.y.right = element_line(color = red_color),
        panel.grid = element_blank(),
        plot.margin = margin(10, 20, 10, 10))

ggsave("figure_2_departure_and_reserves.png", fig2, width = 7, height = 5, units = "in", dpi = 300)
message("Wrote figure_2_departure_and_reserves.png")
