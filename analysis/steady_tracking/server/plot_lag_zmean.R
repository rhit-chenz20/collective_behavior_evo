#!/usr/bin/env Rscript
# Lag (z_opt = optimum - z_mean) and mean phenotype vs. the moving optimum, for every
# parameter combination of the steady_tracking run (shift_speed x reg x psi, 60 reps each).
#
# Input (from analysis/load_data/server/load_data.R -d steady_tracking -t data,metadata --tag all):
#   figure/r_data/steady_tracking/data_summary_all.rds, metadata_all.rds
# Output: analysis/steady_tracking/server/plots_fit500/
#   lag_summary.pdf          panel per reg x speed, line per psi: mean lag over reps,
#                            ribbon = 10-90% quantiles
#   lag_reps_<reg>.pdf       panel per psi x speed, every rep (blue tracked / red failed: mean
#                            fitness over the last 500 ticks < 10% of optimum, classify_tracking.R)
#   zmean_opt_<reg>.pdf      panel per psi x speed: every rep's z_mean, optimum in black
#                            (y axis free per panel: the optimum ends at 10000 + 2000 x speed)
#   zmean_opt_zoom_<reg>.pdf the same for ticks in zoom_ticks, where the lag builds up
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_lag_zmean.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")   # read_summary(), read_meta(): both runs

# ---- settings ----
data_dir   <- "figure/r_data/steady_tracking"
out_dir    <- "analysis/steady_tracking/server/plots_fit500"
shift_tick <- 10000
tick_step  <- 5                  # per-rep figures: every 5th tick
zoom_ticks <- c(9950, 10300)     # window for the zoomed z_mean figures (every tick)
lag_zoom   <- 250                # y-axis cut for lag_summary_zoom.pdf

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_summary()
W <- speed_width(uniqueN(d$shift_speed))   # page width: one facet column per speed
d <- d[, .(shift_speed, reg, psi, rep, ID, tick, optimum, z_mean, lag = z_opt)]
m <- read_meta()
m <- read_class()
# lag at which the end fitness is at the cutoff: median lag_500 of the sims with fit_500 within 0.1 +- 0.01
fail_lag <- m[abs(fit_500 - 0.1) <= 0.01, median(lag_500)]
cat(sprintf("lag at the 10%% cutoff: %.1f\n", fail_lag))

d <- merge(d, m[, .(ID, failed)], by = "ID")
speeds <- sort(unique(d$shift_speed))
d[, `:=`(psi_f = factor(psi),
         speed_f = factor(shift_speed, levels = speeds, labels = paste("speed", speeds)),
         reg = reg_factor(reg),
         outcome = factor(ifelse(failed, "failed", "tracked"), levels = c("tracked", "failed")))]

out_cols <- c(tracked = "#2166ac", failed = "#b2182b")
base_theme <- list(
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_text(angle = 45, hjust = 1))
)
shift_line <- geom_vline(xintercept = shift_tick, linetype = "dotted", colour = "grey40")
lag_line   <- geom_hline(yintercept = fail_lag, linetype = "dashed", colour = "grey40")
row_lab    <- labeller(psi_f = function(x) paste("psi", x))
nf <- unique(d[, .(reg, psi_f, speed_f, ID, failed)])[, .(n_failed = sum(failed), n = .N),
                                                         by = .(reg, psi_f, speed_f)]
fail_text <- function(lab) geom_text(data = lab, aes(x = Inf, y = Inf, label = sprintf("%d/%d failed", n_failed, n)),
                                     inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 2.8)
rep_colours <- list(scale_colour_manual(NULL, values = out_cols),
                    guides(colour = guide_legend(override.aes = list(linewidth = 1, alpha = 1))))

# ---- 1. lag summary ----
s <- d[, .(mean = mean(lag), lo = quantile(lag, 0.1), hi = quantile(lag, 0.9)),
       by = .(reg, speed_f, psi_f, tick)]
p <- ggplot(s, aes(tick, mean, colour = psi_num(psi_f), fill = psi_num(psi_f), group = psi_f)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.5) +
  shift_line + lag_line +
  facet_grid(reg ~ speed_f, scales = "free_y") +
  psi_scale() +
  labs(x = "tick", y = "lag (optimum - mean phenotype)",
       title = "Lag over time (line: mean of 60 reps, ribbon: 10-90% quantiles; y axis free per row)",
       subtitle = sprintf("dotted: optimum starts moving; dashed: lag at which mean fitness is 10%% of optimum (%.0f)", fail_lag)) +
  base_theme
ggsave(file.path(out_dir, "lag_summary.pdf"), p, width = W, height = 7)
# same, y axis cut at lag_zoom so the lags of populations that keep up can be compared
p_zoom <- p + facet_grid(reg ~ speed_f) + coord_cartesian(ylim = c(NA, lag_zoom)) +
  labs(title = sprintf("Lag over time, y axis cut at %d (lines above it run off the top)", lag_zoom))
ggsave(file.path(out_dir, "lag_summary_zoom.pdf"), p_zoom, width = W, height = 7)

# ---- per-rep figures, one file per reg ----
r <- d[(tick - min(tick)) %% tick_step == 0 | tick == max(tick)]
z <- d[tick >= zoom_ticks[1] & tick <= zoom_ticks[2]]
opt_of <- function(dat) unique(dat[, .(reg, psi_f, speed_f, tick, optimum)])

for (g in levels(d$reg)) {
  rg <- r[reg == g]; lab <- nf[reg == g]

  # 2. lag of every rep
  p <- ggplot(rg, aes(tick, lag, group = ID, colour = outcome)) +
    geom_line(linewidth = 0.2, alpha = 0.35) +
    shift_line + lag_line + fail_text(lab) +
    facet_grid(psi_f ~ speed_f, labeller = row_lab, scales = "free_y") +
    rep_colours +
    labs(x = "tick", y = "lag (optimum - mean phenotype)",
         title = sprintf("Lag of every rep, %s (rows: psi; y axis free per row)", g),
         subtitle = sprintf("dashed: lag at which mean fitness is 10%% of optimum (%.0f)", fail_lag)) +
    base_theme
  ggsave(file.path(out_dir, sprintf("lag_reps_%s.pdf", reg_code(g))), p, width = W, height = 9)

  # 3. z_mean of every rep with the optimum
  zmean_plot <- function(dat, title) {
    ggplot(dat, aes(tick, z_mean)) +
      geom_line(aes(group = ID, colour = outcome), linewidth = 0.2, alpha = 0.35) +
      geom_line(data = opt_of(dat), aes(tick, optimum), colour = "black", linewidth = 0.5) +
      shift_line + fail_text(lab) +
      facet_wrap(psi_f ~ speed_f, ncol = length(speeds), scales = "free_y",
                 labeller = labeller(psi_f = function(x) paste("psi", x), .multi_line = FALSE)) +
      rep_colours +
      labs(x = "tick", y = "mean phenotype", title = title,
           subtitle = "black: optimum; each panel has its own y axis") +
      base_theme
  }
  ggsave(file.path(out_dir, sprintf("zmean_opt_%s.pdf", reg_code(g))),
         zmean_plot(rg, sprintf("Mean phenotype of every rep and the optimum, %s (rows: psi)", g)),
         width = W, height = 10)
  ggsave(file.path(out_dir, sprintf("zmean_opt_zoom_%s.pdf", reg_code(g))),
         zmean_plot(z[reg == g], sprintf("Mean phenotype and the optimum, ticks %d-%d, %s (rows: psi)",
                                          zoom_ticks[1], zoom_ticks[2], g)),
         width = W, height = 10)
}
cat("plots written to", out_dir, "\n")
