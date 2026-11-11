# ============================================================
# Figure 3 - Maternal annual time budget by maternal length: days on the
#            breeding ground (Hawaii), in migration, and on the feeding
#            ground (Alaska).
#
# Reproduced from the processed departure ages
#   ../data/processed/departure_df.csv
# plus fixed migration parameters (Mate et al. 1998). No raw data needed.
#
# Run: source("make_figure.R")  ->  figure_3_annual_time_budget.png
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

find_data <- function(fname) {
  for (p in c(file.path("..", "data", "processed", fname),
              file.path("data", "processed", fname), fname))
    if (file.exists(p)) return(p)
  stop("Could not find ", fname, " in ../data/processed/")
}
departure_df <- read.csv(find_data("departure_df.csv"))

# --- parameters -------------------------------------------------------------
mig_mean  <- 39    # one-way migration days (Mate et al. 1998)
year_days <- 365

time_budget <- departure_df %>%
  filter(MotherLength %in% seq(10, 15, by = 1)) %>%
  mutate(
    hawaii_mean    = age_mean,
    hawaii_lwr     = age_lwr,
    hawaii_upr     = age_upr,
    migration_mean = mig_mean * 2,
    alaska_mean    = pmax(year_days - hawaii_mean - migration_mean, 0)
  )

budget_long <- time_budget %>%
  dplyr::select(MotherLength, hawaii_mean, migration_mean, alaska_mean) %>%
  pivot_longer(-MotherLength, names_to = "phase", values_to = "days") %>%
  mutate(phase = factor(phase,
                        levels = c("alaska_mean", "migration_mean", "hawaii_mean"),
                        labels = c("Feeding ground\n(Alaska)", "Migration",
                                   "Breeding ground\n(Hawaii)")))

uncertainty_rect <- time_budget %>% mutate(rect_xmin = hawaii_lwr, rect_xmax = hawaii_upr)

phase_colors <- c("Breeding ground\n(Hawaii)" = "#d95f02",
                  "Migration"                 = "#9ecae1",
                  "Feeding ground\n(Alaska)"  = "#2171b5")

fig3 <- ggplot() +
  geom_col(data = budget_long,
           aes(x = days, y = factor(MotherLength), fill = phase),
           position = "stack", width = 0.65, color = "white", linewidth = 0.4) +
  geom_rect(data = uncertainty_rect,
            aes(xmin = rect_xmin, xmax = rect_xmax,
                ymin = as.numeric(factor(MotherLength)) - 0.35,
                ymax = as.numeric(factor(MotherLength)) + 0.35),
            fill = "grey20", alpha = 0.2) +
  scale_fill_manual(values = phase_colors, name = NULL,
                    guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(name = "Days per year", breaks = seq(0, 360, by = 60),
                     expand = expansion(mult = c(0.01, 0.06))) +
  scale_y_discrete(name = "Maternal length (m)") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
        plot.margin = margin(10, 40, 10, 10))

ggsave("figure_3_annual_time_budget.png", fig3, width = 7, height = 5, units = "in", dpi = 300)
message("Wrote figure_3_annual_time_budget.png")
