#!/usr/bin/env Rscript
# Copy-the-average at speed 2.6, psi 0.2 and 0.5: every rep's mean phenotype against the optimum
# (top row, zoomed to phenotype 0-500, i.e. the first ~200 ticks after the shift) and its lag,
# optimum - z_mean (bottom row); one column per psi. Reps are coloured by their outcome from
# classify_tracking.R (failed = mean fitness over the last 500 ticks < 10% of the fitness at the
# optimum, or aborted by SLiM). Dashed red line in the lag row: the lag at the 10% cutoff, the
# median lag_500 of the sims with fit_500 within 0.1 +- 0.01 (as in plot_lag_zmean.R).
#
# Input: figure/r_data/steady_tracking/data_summary_all.rds,
#        plots_fit500/tracking_class/tracking_class_per_sim.csv
# Output: analysis/steady_tracking/server/plots_fit500/trajectories/zmean_lag_ave_speed2.6.pdf
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_zmean_lag_example.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("analysis/steady_tracking/server/steady_data.R")

rule   <- "ave"
speed  <- 2.6
psis   <- c(0.2, 0.5)
z_lim  <- c(0, 500)              # top row: phenotype range shown
out    <- file.path(plot_root, "trajectories", sprintf("zmean_lag_%s_speed%.1f.pdf", rule, speed))
out_cols <- c(tracked = "#2166ac", failed = "#b2182b")

d <- read_summary(cols = c("optimum", "z_mean", "z_opt"), runs = "all", keep_partial = TRUE,
                  speeds = NULL)
d <- d[reg == rule & abs(shift_speed - speed) < 1e-9 & round(psi, 2) %in% psis & tick >= 9900]
cls <- fread(class_file)[dataset == "steady_tracking", .(ID, outcome)]
d <- merge(d, cls, by = "ID")
d[, outcome := factor(fifelse(outcome == "success", "tracked", "failed"), levels = names(out_cols))]
d[, psi_f := factor(psi, levels = psis, labels = sprintf("psi = %.1f", psis))]
n_lab <- d[, .(lab = sprintf("tracked %d / failed %d", uniqueN(ID[outcome == "tracked"]),
                             uniqueN(ID[outcome == "failed"]))), by = psi_f]

# lag at the 10% cutoff
fail_lag <- read_class()[abs(fit_500 - 0.1) <= 0.01, median(lag_500)]
cat(sprintf("lag at the 10%% cutoff: %.1f\n", fail_lag))

# top row: ticks until the optimum passes the top of the zoom
t_max <- ceiling(10000 + z_lim[2] / speed) + 20
top <- d[tick <= t_max]
opt <- unique(top[, .(psi_f, tick, optimum)])

theme_row <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"), panel.grid.minor = element_blank(),
        legend.position = "top")
shift_line <- geom_vline(xintercept = 10000, colour = "grey70", linewidth = 0.3)
rep_cols <- scale_colour_manual(NULL, values = out_cols,
                                labels = c(tracked = "rep that tracked", failed = "rep that failed"),
                                drop = FALSE)

p1 <- ggplot(top, aes(tick, z_mean)) +
  shift_line +
  geom_line(aes(group = ID, colour = outcome), linewidth = 0.25, alpha = 0.4) +
  geom_line(data = opt, aes(y = optimum, linetype = "optimum"), colour = "grey20", linewidth = 0.6) +
  geom_text(data = n_lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.5, size = 3) +
  facet_grid(. ~ psi_f) +
  rep_cols + scale_linetype_manual(NULL, values = c(optimum = "dashed")) +
  coord_cartesian(ylim = z_lim, expand = FALSE) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1, alpha = 1))) +
  labs(x = "tick", y = "mean phenotype (z_mean)") +
  theme_row

p2 <- ggplot(d, aes(tick, z_opt)) +
  shift_line +
  geom_hline(aes(yintercept = fail_lag, linetype = "lag at the 10% fitness cutoff"),
             colour = "#c0392b", linewidth = 0.5) +
  geom_line(aes(group = ID, colour = outcome), linewidth = 0.2, alpha = 0.35, show.legend = FALSE) +
  facet_grid(. ~ psi_f) +
  rep_cols + scale_linetype_manual(NULL, values = c(`lag at the 10% fitness cutoff` = "dashed")) +
  labs(x = "tick", y = "lag (optimum - z_mean)") +
  theme_row

p <- (p1 / p2) +
  plot_annotation(title = sprintf("%s, speed %.1f", reg_names[[rule]], speed),
                  subtitle = sprintf("thin lines: individual reps; top: zoomed to phenotype %d-%d; red dashed: lag %.0f, where mean fitness is 10%% of peak; grey line: shift at tick 10000",
                                     z_lim[1], z_lim[2], fail_lag),
                  theme = theme(plot.title = element_text(face = "bold"),
                                plot.subtitle = element_text(size = 9)))

ggsave(out, p, width = 9, height = 7.5)
cat("wrote", out, "\n")
