#!/usr/bin/env Rscript
# fitness_mean over time for every parameter combination of the steady_tracking run
# (moving optimum, n = 8; shift_speed x reg x psi; speeds 0.5-2.9 with 60 reps and 3.2-5.0 with 20 reps,
# the fast_speeds run; combinations with aborted reps are dropped, see steady_data.R).
#
# Input (made by analysis/load_data/server/load_data.R -d steady_tracking -t data,metadata --tag all):
#   figure/r_data/steady_tracking/data_summary_all.rds   per tick, 9950-12000
#   figure/r_data/steady_tracking/metadata_all.rds       per sim (fit_ratio, for the fitness at the optimum)
#   plots_fit500/tracking_class/tracking_class_per_sim.csv  failed = mean fitness over the last
#                                                         500 ticks < 10% of optimum (classify_tracking.R)
# Output: analysis/steady_tracking/server/plots_fit500/
#   fitness_mean_summary.pdf   one panel per reg x speed, a line per psi: mean over reps,
#                              ribbon = 10-90% quantiles across reps
#   fitness_mean_reps.pdf      one panel per combination (reg x psi rows, speed columns):
#                              each rep as a thin line, coloured by whether it failed
#   fitness_mean_reps_<reg>.pdf  the same, one file per reg (easier to read)
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_fitness_mean.R
# Reading the 1.8 GB data file takes ~30 s.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")   # read_summary(), read_meta(): both runs

# ---- settings ----
data_dir   <- "figure/r_data/steady_tracking"
out_dir    <- "analysis/steady_tracking/server/plots_fit500"
shift_tick <- 10000   # optimum starts moving here
tick_step  <- 5       # plot every 5th tick in the per-rep figures (2051 -> ~410 points per sim)
fail_fit_ratio <- 0.1 # a tick fails if fitness_mean < this x fitness at the optimum

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_summary()
W <- speed_width(uniqueN(d$shift_speed))   # page width: one facet column per speed
d <- d[, .(shift_speed, reg, psi, rep, ID, tick, fitness_mean)]
m <- read_meta()
m <- merge(m[, .(ID, fit_ratio)], read_class()[, .(ID, failed)], by = "ID")

# fitness at the optimum (the same in every sim), from fit_ratio = fitness_mean(12000) / fit_opt
end <- merge(d[tick == max(tick), .(ID, fitness_mean)], m, by = "ID")
fit_opt <- median(end$fitness_mean / end$fit_ratio)
fail_line <- fail_fit_ratio * fit_opt
cat(sprintf("fitness at optimum %.6g, fail threshold %.6g\n", fit_opt, fail_line))

d <- merge(d, m[, .(ID, failed)], by = "ID")
d[, psi_f := factor(psi)]
d[, speed_f := factor(shift_speed, levels = sort(unique(shift_speed)),
                      labels = paste("speed", sort(unique(shift_speed))))]
d[, reg := reg_factor(reg)]

common <- list(
  geom_vline(xintercept = shift_tick, linetype = "dotted", colour = "grey40"),
  geom_hline(yintercept = fail_line, linetype = "dashed", colour = "grey40"),
  labs(x = "tick", y = "mean fitness"),
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_text(angle = 45, hjust = 1))
)

# ---- 1. summary: mean over reps per combination ----
s <- d[, .(mean = mean(fitness_mean),
           lo = quantile(fitness_mean, 0.1), hi = quantile(fitness_mean, 0.9),
           frac_failed = mean(failed)),
       by = .(reg, speed_f, psi_f, tick)]
p1 <- ggplot(s, aes(tick, mean, colour = psi_num(psi_f), fill = psi_num(psi_f), group = psi_f)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.5) +
  facet_grid(reg ~ speed_f) +
  psi_scale() +
  common +
  labs(title = "Mean fitness over time (line: mean of 60 reps, ribbon: 10-90% quantiles)",
       subtitle = sprintf("dotted: optimum starts moving (tick %d); dashed: %.1f x fitness at optimum; failed = mean fitness over the last 500 ticks below it",
                          shift_tick, fail_fit_ratio))
ggsave(file.path(out_dir, "fitness_mean_summary.pdf"), p1, width = W, height = 7)

# ---- 2. every rep, one panel per combination ----
r <- d[(tick - min(tick)) %% tick_step == 0 | tick == max(tick)]
r[, outcome := factor(ifelse(failed, "failed", "tracked"), levels = c("tracked", "failed"))]
# panel labels with how many reps failed
nf <- unique(d[, .(reg, psi_f, speed_f, ID, failed)])[, .(n_failed = sum(failed), n = .N),
                                                         by = .(reg, psi_f, speed_f)]
rep_plot <- function(dat, lab, rows) {
  ggplot(dat, aes(tick, fitness_mean, group = ID, colour = outcome)) +
    geom_line(linewidth = 0.2, alpha = 0.35) +
    geom_text(data = lab, aes(x = Inf, y = Inf, label = sprintf("%d/%d failed", n_failed, n)),
              inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 2.8) +
    facet_grid(rows, labeller = labeller(psi_f = function(x) paste("psi", x),
                                         reg = label_value)) +
    scale_colour_manual(NULL, values = c(tracked = "#2166ac", failed = "#b2182b")) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1, alpha = 1))) +
    common
}
p2 <- rep_plot(r, nf, reg + psi_f ~ speed_f) +
  labs(title = "Mean fitness of every rep, one panel per combination (rows: rule x psi)")
ggsave(file.path(out_dir, "fitness_mean_reps.pdf"), p2, width = W, height = 22)

for (g in levels(d$reg)) {
  p3 <- rep_plot(r[reg == g], nf[reg == g], psi_f ~ speed_f) +
    labs(title = sprintf("Mean fitness of every rep, %s (rows: psi)", g))
  ggsave(file.path(out_dir, sprintf("fitness_mean_reps_%s.pdf", reg_code(g))), p3, width = W, height = 9)
}

fwrite(s[tick %in% c(shift_tick, shift_tick + 500, shift_tick + 1000, max(tick))],
       file.path(out_dir, "fitness_mean_summary_selected_ticks.csv"))
cat("plots written to", out_dir, "\n")
