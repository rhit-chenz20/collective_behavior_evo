#!/usr/bin/env Rscript
# Trajectories behind the tracking classification (classify_tracking.R, 10% cutoff).
#   page 1    a clear-cut sim the metadata calls failed (mean fitness dipped below 10% after the
#             shift) that the equilibrium classification calls tracking
#   close calls: every sim whose fit_500 is within +-0.01 of the cutoff, one per page
# Each page, from tick 9900: (top) mean phenotype z_mean and the optimum, (middle) the lag,
# optimum - z_mean, (bottom) log10 mean fitness / fitness at the optimum with the 10% cutoff.
# The shaded band is the 500-tick window whose mean fitness decides the class.
#
# Input: figure/r_data/steady_tracking{,_long}/data_summary_*.rds,
#        analysis/steady_tracking/server/plots_fit500/tracking_class/tracking_class_per_sim.csv
# Output: analysis/steady_tracking/server/plots_fit500/tracking_class/
#   example_fail_to_success.pdf, close_calls_trajectories.pdf
#
# Run from the repo root (~1 min):
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_tracking_class_examples.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")

out_dir    <- "analysis/steady_tracking/server/plots_fit500/tracking_class"
example_id <- 4184     # 2.3 ave psi 0.5 rep 39: 421 ticks below 10%, ends at ~18%, steady lag
close_tol  <- 0.01     # close calls: |fit_500 - cutoff| <= close_tol
cutoff     <- 0.1
peak       <- dnorm(0, 0, sqrt(2000))
col_z      <- "#2a6fdb"
col_opt    <- "grey45"

cls <- keep_speeds(fread(file.path(out_dir, "tracking_class_per_sim.csv")))   # plot_speeds
unres <- cls[class != "aborted" & abs(fit_500 - cutoff) <= close_tol][order(dataset, reg, psi, shift_speed, rep)]
want  <- rbind(cls[dataset == "steady_tracking" & ID == example_id], unres)

cols <- c("ID", "tick", "optimum", "z_mean", "z_opt", "fitness_mean")
load_ids <- function(f, ids) as.data.table(readRDS(f))[ID %in% ids & tick >= 9900, ..cols]
ids_st <- want[dataset == "steady_tracking", ID]
d_st <- rbindlist(lapply(st_files, function(f) load_ids(file.path(st_data_dir, f[1]), ids_st)))
d_lg <- load_ids("figure/r_data/steady_tracking_long/data_summary_all.rds", want[dataset == "long", ID])
traj <- rbind(d_st[, dataset := "steady_tracking"], d_lg[, dataset := "long"])

plot_sim <- function(row) {
  x <- traj[dataset == row$dataset & ID == row$ID]
  end_tick <- max(x$tick)
  long <- rbind(
    x[, .(tick, panel = "phenotype", series = "z_mean", value = z_mean)],
    x[, .(tick, panel = "phenotype", series = "optimum", value = optimum)],
    x[, .(tick, panel = "lag (optimum - z_mean)", series = "z_mean", value = z_opt)],
    x[, .(tick, panel = "log10 fitness / peak", series = "z_mean", value = log10(fitness_mean / peak))])
  long[, panel := factor(panel, levels = c("phenotype", "lag (optimum - z_mean)", "log10 fitness / peak"))]
  cut_line <- data.table(panel = factor("log10 fitness / peak", levels = levels(long$panel)), y = log10(cutoff))
  meta_txt <- if (isTRUE(row$meta_failed)) sprintf("metadata: failed (%d ticks < 10%%)", row$n_fail_ticks) else "metadata: not failed"
  ggplot(long, aes(tick, value, colour = series, linetype = series)) +
    annotate("rect", xmin = end_tick - 500, xmax = end_tick, ymin = -Inf, ymax = Inf, fill = "grey92") +
    geom_vline(xintercept = 10000, colour = "grey70", linewidth = 0.3) +
    geom_hline(data = cut_line, aes(yintercept = y), colour = "#c0392b", linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 0.4) +
    facet_grid(panel ~ ., scales = "free_y", switch = "y") +
    scale_colour_manual(values = c(z_mean = col_z, optimum = col_opt), name = NULL) +
    scale_linetype_manual(values = c(z_mean = "solid", optimum = "dashed"), name = NULL) +
    labs(x = "tick", y = NULL,
         title = sprintf("%s, speed %.1f, psi %.1f, rep %d (ID %d, %s)",
                         reg_names[row$reg], row$shift_speed, row$psi, row$rep, row$ID, row$dataset),
         subtitle = sprintf("class: %s (%s)   |   %s\nlast 500 ticks: mean fitness %.3f of peak, lag %.1f, lag slope %.1f%% of speed",
                            row$class, row$outcome, meta_txt, row$fit_500, row$lag_500,
                            100 * row$lag_slope_500 / row$shift_speed),
         caption = "shaded: last 500 ticks (classification window); red dashed: 10% cutoff; grey line: shift at tick 10000") +
    theme_bw(base_size = 10) +
    theme(strip.placement = "outside", strip.background = element_blank(), legend.position = "top",
          panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
}

pdf(file.path(out_dir, "example_fail_to_success.pdf"), width = 8, height = 8)
print(plot_sim(want[dataset == "steady_tracking" & ID == example_id]))
invisible(dev.off())

pdf(file.path(out_dir, "close_calls_trajectories.pdf"), width = 8, height = 8)
for (i in seq_len(nrow(unres))) print(plot_sim(unres[i]))
invisible(dev.off())
cat("wrote", nrow(unres), "close-call pages and the example to", out_dir, "\n")
