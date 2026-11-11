# Figure 2 - Calf departure age and maternal reserves

Reproduces the combined figure: predicted calf age at the 6 m departure size
(blue, left axis) and maternal blubber reserves remaining at departure (orange,
right axis), versus maternal length, with the "depleted reserves" zone and the
mean maternal length (12.6 m) marked.

**Run:** `source("make_figure.R")`.

**Input:** `../data/processed/departure_df.csv` (blue curve). The reserves curve is
recomputed from deterministic body-volume allometry (hardcoded published
coefficients; van Aswegen et al. 2025), so no other data are needed.
**Output:** `figure_2_departure_and_reserves.png` (written here)
**Packages:** dplyr, ggplot2, tibble
