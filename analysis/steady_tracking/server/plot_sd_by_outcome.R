#!/usr/bin/env Rscript
# SD of genotypic value (g_sd) and SD of phenotype (z_sd) over time, split by tracking outcome
# (classify_tracking.R): one file with only the reps that tracked (success), one with only the reps
# that failed. Same layout as plot_trajectories.R: panel per rule x speed, line per psi = mean over
# the reps of that outcome, ribbon = 10-90% quantiles. Each panel lists how many reps per psi
# went in.
#
# Failed = too_unfit only: sims aborted by SLiM stop early, and averaging a set of reps that shrinks
# over time makes jumps in the mean, so they are left out (their count is printed and written).
# Combinations with aborted reps are kept (unlike read_summary()'s default).
#
# Input: figure/r_data/steady_tracking/data_summary_*.rds,
#        plots_fit500/tracking_class/tracking_class_per_sim.csv
# Output: analysis/steady_tracking/server/plots_fit500/trajectories/
#   sd_success.pdf, sd_failed.pdf   page 1 g_sd, page 2 z_sd
#   sd_by_outcome_counts.csv        reps per combination and outcome (incl. aborted)
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_sd_by_outcome.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")

out_dir    <- file.path(plot_root, "trajectories")
shift_tick <- 10000
vars <- c(g_sd = "SD of genotypic value", z_sd = "SD of phenotype")

d <- read_summary(cols = names(vars), keep_partial = TRUE)
cls <- keep_speeds(fread(class_file)[dataset == "steady_tracking", .(ID, class, outcome)])
d <- merge(d, cls[class != "aborted"], by = "ID")

counts <- merge(unique(d[, .(ID, shift_speed, reg, psi)]), cls, by = "ID")[
  , .N, by = .(shift_speed, reg, psi, class)]
counts <- rbind(counts, keep_speeds(read_crashed())[, .N, by = .(shift_speed, reg, psi)][, class := "aborted"])
fwrite(dcast(counts, shift_speed + reg + psi ~ class, value.var = "N", fill = 0)[order(reg, psi, shift_speed)],
       file.path(out_dir, "sd_by_outcome_counts.csv"))

speeds <- sort(unique(d$shift_speed))
W <- speed_width(length(speeds))
d[, `:=`(psi_f = factor(psi), speed_f = factor(shift_speed, levels = speeds, labels = paste("speed", speeds)),
         reg = reg_factor(reg))]

base_theme <- list(
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_text(angle = 45, hjust = 1))
)

sd_plot <- function(x, v, label, what) {
  s <- x[!is.na(get(v)), .(mean = mean(get(v)), lo = quantile(get(v), 0.1), hi = quantile(get(v), 0.9)),
         by = .(reg, speed_f, psi_f, tick)]
  # reps per psi in each panel (0 where none of that outcome)
  n <- CJ(reg = levels(x$reg), speed_f = levels(x$speed_f), psi_f = levels(x$psi_f))
  n <- merge(n, x[, .(n = uniqueN(ID)), by = .(reg = as.character(reg), speed_f = as.character(speed_f),
                                                  psi_f = as.character(psi_f))],
             all.x = TRUE)[is.na(n), n := 0L]
  n <- n[, .(lab = paste0("n = ", paste(n[order(as.numeric(psi_f))], collapse = "/"))), by = .(reg, speed_f)]
  n[, `:=`(reg = factor(reg, levels = levels(x$reg)), speed_f = factor(speed_f, levels = levels(x$speed_f)))]
  ggplot(s, aes(tick, mean, colour = psi_num(psi_f), fill = psi_num(psi_f), group = psi_f)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.5) +
    geom_vline(xintercept = shift_tick, linetype = "dotted", colour = "grey40") +
    geom_text(data = n, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = -0.05, vjust = 1.4, size = 2.2) +
    facet_grid(reg ~ speed_f) +
    psi_scale() +
    labs(x = "tick", y = label, title = sprintf("%s  [%s], %s", label, v, what),
         subtitle = sprintf("line: mean over the %s, ribbon: 10-90%% quantiles; n = reps per psi (%s); dotted: optimum starts moving",
                            what, paste(levels(x$psi_f), collapse = "/"))) +
    base_theme
}

sets <- list(success = list(d[outcome == "success"], "reps that tracked (success)"),
             failed  = list(d[outcome == "failed"],  "reps that failed (too unfit; aborted sims left out)"))
for (k in names(sets)) {
  x <- sets[[k]][[1]]
  pdf(file.path(out_dir, sprintf("sd_%s.pdf", k)), width = W, height = 7)
  for (v in names(vars)) print(sd_plot(x, v, vars[[v]], sets[[k]][[2]]))
  invisible(dev.off())
  cat(k, ":", uniqueN(x$ID), "sims\n")
}
cat("aborted sims left out:", counts[class == "aborted", sum(N)], "\n")
cat("plots written to", out_dir, "\n")
