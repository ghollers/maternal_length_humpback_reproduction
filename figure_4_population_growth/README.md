# Figure 4 - Population growth rate vs maternal asymptotic length

Reproduces the Leslie-matrix figure: intrinsic population growth rate (lambda) as
a function of female asymptotic length, under length-decline scenarios, with
Monte Carlo uncertainty on the demographic parameters.

**Run:** `source("make_figure.R")`.

**Input:** none. Self-contained - it sources `fecundity_function.R` (in this
folder) and uses published/hardcoded demographic parameters.
**Output:** `figure_4_population_growth.png` and `figure_4_scenarios_lambda.csv`
(written here)
**Packages:** truncnorm, ggplot2, dplyr, patchwork
