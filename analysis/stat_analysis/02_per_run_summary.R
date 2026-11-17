# 2. per-run summary table
# Section 2 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/02_per_run_summary.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
source("analysis/stat_analysis/00_load_data.R")

# Collapse each run's trajectory to single numbers: from here on the run is the
# unit of analysis, never the generation.
summarize_run <- function(d) {
  data.frame(
    model            = as.character(d$model[1]),
    run              = d$run[1],
    mean_sd_all      = mean(d[[sd_col]],      na.rm = TRUE),
    mean_mean_all    = mean(d[[mean_col]],    na.rm = TRUE),
    psi              = d$psi[1],
    n                = d$n[1],
    rep              = d$rep[1],
    param_set        = d$param_set[1],
    n_timepoints     = nrow(d)
  )
}
run_df <- do.call(rbind, lapply(split(burn_df, burn_df$run), summarize_run))
run_df <- run_df[order(run_df$psi, run_df$n, match(run_df$model, models), run_df$rep), ]
rownames(run_df) <- NULL

# Sanity check: every run should have 50 recorded time points.
if (any(run_df$n_timepoints != 50)) {
  warning("Some runs do not have 50 time points: ",
          paste(unique(run_df$n_timepoints), collapse = ", "))
}

write.csv(run_df, out_path("per_run_summary.csv"), row.names = FALSE)
print(run_df, row.names = FALSE)
