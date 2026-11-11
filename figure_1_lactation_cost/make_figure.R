# ============================================================
# Figure 1 - Mass-specific lactation cost and % of daily energy
#            loss transferred to the calf, vs maternal length.
#
# Reproduces the figure from the processed dataset
#   ../data/processed/lactation_estimates.csv
# (the Monte Carlo energy-budget output). No raw data needed.
#
# Run: open reviewer_package in RStudio (or setwd here), then
#   source("make_figure.R")
# Output: figure_1_lactation_cost.png (written next to this script).
# ============================================================

library(dplyr)
library(ggplot2)

# --- locate the processed data regardless of working directory --------------
find_data <- function(fname) {
  for (p in c(file.path("..", "data", "processed", fname),
              file.path("data", "processed", fname),
              fname)) if (file.exists(p)) return(p)
  stop("Could not find ", fname, " in ../data/processed/")
}
lactation <- read.csv(find_data("lactation_estimates.csv"))

# --- colours and dual-axis mapping ------------------------------------------
blue_color  <- "#0072B2"   # % of daily energy loss to calf (left axis)
red_color   <- "#D55E00"   # mass-specific lactation cost   (right axis)
y_left_min  <- 35;   y_left_max  <- 70
y_right_min <- 0.02; y_right_max <- 0.12
scale_factor <- (y_left_max - y_left_min) / (y_right_max - y_right_min)
offset       <- y_left_min - y_right_min * scale_factor
to_left  <- function(x) x * scale_factor + offset   # right-axis value -> left scale
to_right <- function(y) (y - offset) / scale_factor
smooth_vec <- function(x, y, span = 0.5) predict(loess(y ~ x, span = span))

# --- build the two outcomes at three resting proportions --------------------
p_rest_show <- c(0.2, 0.4, 0.8)   # solid = empirical 0.4; dashed = 0.2 and 0.8

lac_lines <- lactation %>%
  filter(Age == 50, p_rest %in% p_rest_show) %>%
  mutate(
    pct_to_calf = (lactation_mean / energy_loss_mean_MJ) * 100,
    lac_per_kg  = lactation_mean / MotherMass_kg
  ) %>%
  group_by(p_rest) %>% arrange(MotherLength) %>%
  mutate(pct_s = smooth_vec(MotherLength, pct_to_calf),
         kg_s  = smooth_vec(MotherLength, lac_per_kg)) %>%
  ungroup()

lines_long <- bind_rows(
  lac_lines %>% transmute(MotherLength, p_rest,
                          outcome = "% of total daily energy loss to calf", y = pct_s),
  lac_lines %>% transmute(MotherLength, p_rest,
                          outcome = "Mass-specific lactation cost",          y = to_left(kg_s))
) %>% mutate(emph = p_rest == 0.4)

line_ends <- lines_long %>% group_by(outcome, p_rest) %>%
  slice_max(MotherLength, n = 1) %>% ungroup()

# --- figure -----------------------------------------------------------------
fig1 <- ggplot(lines_long,
               aes(MotherLength, y, group = interaction(outcome, p_rest),
                   color = outcome, linetype = emph, linewidth = emph)) +
  geom_line() +
  geom_text(data = line_ends, aes(label = sprintf("%.1f", p_rest)),
            hjust = 0, nudge_x = 0.07, size = 3, show.legend = FALSE) +
  scale_linetype_manual(values = c(`FALSE` = "dashed", `TRUE` = "solid"), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.7, `TRUE` = 1.5), guide = "none") +
  scale_color_manual(name = NULL,
    values = c("% of total daily energy loss to calf" = blue_color,
               "Mass-specific lactation cost"          = red_color)) +
  scale_y_continuous(
    name = "% of total daily energy loss to calf",
    sec.axis = sec_axis(~ to_right(.),
      name   = expression("Mass-specific lactation cost (MJ day"^{-1}*" kg"^{-1}*")"),
      breaks = seq(y_right_min, y_right_max, by = 0.02))) +
  scale_x_continuous(name = "Maternal length (m)", breaks = seq(10, 15, 1),
                     expand = expansion(mult = c(0.02, 0.10))) +
  coord_cartesian(ylim = c(y_left_min, y_left_max)) +
  labs(caption = "Solid = empirical resting proportion (p_rest = 0.4); dashed = p_rest 0.2 and 0.8. Labels give p_rest.") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        axis.title.y.left  = element_text(color = blue_color, size = 12),
        axis.text.y.left   = element_text(color = blue_color),
        axis.title.y.right = element_text(color = red_color,  size = 12),
        axis.text.y.right  = element_text(color = red_color),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10))

ggsave("figure_1_lactation_cost.png", fig1, width = 7, height = 5, units = "in", dpi = 300)
message("Wrote figure_1_lactation_cost.png")
