# reproduce_all_figures.R
# Runs the four figure scripts in turn. Each writes its PNG into its own folder.
# Open reviewer_package.Rproj first (or setwd() to this folder), then:
#   source("reproduce_all_figures.R")

root <- getwd()
figs <- c("figure_1_lactation_cost",
          "figure_2_departure_and_reserves",
          "figure_3_annual_time_budget",
          "figure_4_population_growth")

for (f in figs) {
  message("\n=========== ", f, " ===========")
  old <- setwd(file.path(root, f))
  tryCatch(source("make_figure.R"), finally = setwd(old))
}

message("\nDone. Each figure_*.png is in its own folder.")
