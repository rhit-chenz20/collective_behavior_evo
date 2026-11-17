#!/usr/bin/env Rscript
# Heterozygosity H against the lag (optimum - mean phenotype) for every parameter combination
# of the steady_tracking run (speeds 0.5-5.0 x 3 reg x 4 psi; see steady_data.R). Tests whether the loss of H
# is set by how far the population falls behind, whatever the speed / reg / psi.
#
# Input: figure/r_data/steady_tracking/data_summary_all.rds (per tick: H, z_opt = lag)
#        plots_fit500/tracking_class/tracking_class_per_sim.csv (failed: mean fitness over the
#        last 500 ticks < 10% of optimum, classify_tracking.R)
# Output: analysis/steady_tracking/server/plots_fit500/H_vs_lag/
#   H_vs_lag_end.pdf      one point per sim: mean H and lag over ticks end_window
#                         page 1: all sims in one panel, colour = speed, shape = reg
#                         page 2: panel per reg x psi, colour = speed
#                         page 3: panel per reg x psi, colour = failed or not
#   H_vs_lag_windows.pdf  one point per sim x 100-tick window from window_from on, panel per
#                         reg x psi, colour = window (time); shows the relation while it changes
#   H_vs_lag_paths.pdf    mean over reps per combination and 50-tick window, drawn as a path
#                         from the shift to tick 12000 (dot = end), panel per reg x psi,
#                         colour = speed
# The lag axis is log10: lags run from ~10 to ~1500.
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_H_vs_lag.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")   # read_summary(), read_meta(): both runs

# ---- settings ----
data_dir    <- "figure/r_data/steady_tracking"
out_dir     <- "analysis/steady_tracking/server/plots_fit500/H_vs_lag"
shift_tick  <- 10000
end_window  <- c(11800, 12000)  # end state: mean over these ticks
window_from <- 10200            # time windows start here (skips the first rise of the lag)
window_len  <- 100              # ticks per window in H_vs_lag_windows.pdf
path_len    <- 50               # ticks per window in H_vs_lag_paths.pdf

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_summary()
d <- d[tick >= shift_tick, .(shift_speed, reg, psi, rep, ID, tick, H, lag = z_opt)]
m <- read_class()[, .(ID, failed)]

speeds <- sort(unique(d$shift_speed))
fmt <- function(x) {
  x[, `:=`(speed_f = factor(shift_speed, levels = speeds),
           reg = factor(reg, levels = c("ave", "ext", "fit")),
           psi_f = factor(psi, labels = paste("psi", sort(unique(psi)))))]
  x
}
speed_cols <- scale_colour_viridis_d("speed", option = "plasma", end = 0.9)
base_theme <- list(
  scale_x_log10(),
  annotation_logticks(sides = "b", linewidth = 0.2),
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"))
)
lab <- labs(x = "lag (optimum - mean phenotype), log scale", y = "heterozygosity H")

# ---- 1. end state, one point per sim ----
e <- d[tick >= end_window[1] & tick <= end_window[2],
       .(H = mean(H), lag = mean(lag)), by = .(shift_speed, reg, psi, ID)]
e <- fmt(merge(e, m, by = "ID"))
sub_e <- sprintf("one point per sim (%d sims): mean over ticks %d-%d", nrow(e), end_window[1], end_window[2])
pdf(file.path(out_dir, "H_vs_lag_end.pdf"), width = 11, height = 8)
print(ggplot(e, aes(lag, H, colour = speed_f, shape = reg)) +
        geom_point(size = 1, alpha = 0.6) +
        speed_cols + scale_shape_manual("reg", values = c(ave = 16, ext = 17, fit = 15)) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
               shape = guide_legend(override.aes = list(size = 3))) +
        lab + labs(title = "H vs lag at the end of the run, all combinations", subtitle = sub_e) +
        base_theme)
print(ggplot(e, aes(lag, H, colour = speed_f)) +
        geom_point(size = 0.8, alpha = 0.6) +
        facet_grid(reg ~ psi_f) +
        speed_cols + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        lab + labs(title = "H vs lag at the end of the run, by reg and psi", subtitle = sub_e) +
        base_theme)
print(ggplot(e, aes(lag, H, colour = ifelse(failed, "failed", "tracked"))) +
        geom_point(size = 0.8, alpha = 0.6) +
        facet_grid(reg ~ psi_f) +
        scale_colour_manual(NULL, values = c(tracked = "#2166ac", failed = "#b2182b")) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        lab + labs(title = "H vs lag at the end of the run, coloured by whether the sim failed", subtitle = sub_e) +
        base_theme)
invisible(dev.off())

# ---- 2. time windows per sim ----
w <- d[tick >= window_from][, win := window_from + ((tick - window_from) %/% window_len) * window_len + window_len / 2]
w <- w[, .(H = mean(H), lag = mean(lag)), by = .(shift_speed, reg, psi, ID, win)]
w <- fmt(w)
p <- ggplot(w, aes(lag, H, colour = win)) +
  geom_point(size = 0.3, alpha = 0.3, shape = 16) +
  facet_grid(reg ~ psi_f) +
  scale_colour_viridis_c("tick (window\nmidpoint)", option = "viridis") +
  lab + labs(title = "H vs lag over time, all speeds",
             subtitle = sprintf("one point per sim and %d-tick window from tick %d", window_len, window_from)) +
  base_theme
ggsave(file.path(out_dir, "H_vs_lag_windows.pdf"), p, width = 13, height = 8)

# ---- 3. mean path per combination ----
pa <- d[, win := shift_tick + (tick - shift_tick) %/% path_len * path_len + path_len / 2]
pa <- pa[, .(H = mean(H), lag = mean(lag)), by = .(shift_speed, reg, psi, win)]
pa <- fmt(pa[lag > 0])
setorder(pa, shift_speed, reg, psi, win)
p <- ggplot(pa, aes(lag, H, colour = speed_f, group = speed_f)) +
  geom_path(linewidth = 0.5) +
  geom_point(data = pa[, .SD[.N], by = .(speed_f, reg, psi_f)], size = 1.5) +
  facet_grid(reg ~ psi_f) +
  speed_cols +
  lab + labs(title = "Mean path of each combination through H-lag space",
             subtitle = sprintf("mean of 60 reps per %d-tick window, from the shift to tick 12000 (dot = end)", path_len)) +
  base_theme
ggsave(file.path(out_dir, "H_vs_lag_paths.pdf"), p, width = 13, height = 8)

fwrite(e, file.path(out_dir, "H_vs_lag_end_per_sim.csv"))
cat("plots written to", out_dir, "\n")
