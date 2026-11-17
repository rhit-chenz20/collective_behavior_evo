#!/usr/bin/env Rscript
# The long_borderline run: 11 combinations that were not settled at tick 12000 in
# steady_tracking, rerun (20 new reps each) to tick 15000. Do they settle at a steady lag, or
# run away?
#
# Input:
#   figure/r_data/steady_tracking_long/data_summary_all.rds   (load_data.R -d steady_tracking_long -t data,metadata --tag all)
#   figure/r_data/steady_tracking/data_summary_all.rds        same combinations to 12000, 60 reps (for comparison)
# Output: analysis/steady_tracking/server/plots_fit500/long_borderline/
#   long_borderline.pdf   one page per quantity (lag, response ratio, fitness_mean, z_sd);
#                         one panel per combination: each rep thin, mean thick; dashed grey =
#                         mean of the 60 steady_tracking reps (to 12000)
#   long_borderline_trend.pdf  per 1000-tick window: lag gain per tick / speed (0 = steady), per
#                         combination
#   long_borderline_windows.csv
#
# Run from the repo root (~1 min):
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_long_borderline.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir    <- "analysis/steady_tracking/server/plots_fit500/long_borderline"
shift_tick <- 10000
tick_step  <- 10
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols <- c("shift_speed", "reg", "psi", "rep", "ID", "tick", "z_opt", "fitness_mean", "z_sd", "dg_obs")
source("analysis/steady_tracking/server/steady_data.R")   # keep_speeds(): plot_speeds
L <- keep_speeds(as.data.table(readRDS("figure/r_data/steady_tracking_long/data_summary_all.rds"))[, ..cols])
combos <- unique(L[, .(shift_speed, reg, psi)])
O <- as.data.table(readRDS("figure/r_data/steady_tracking/data_summary_all.rds"))[, ..cols]
O <- O[combos, on = .(shift_speed, reg, psi), nomatch = 0]

lab_of <- function(x) sprintf("%s psi %s, speed %s", reg_names[x$reg], x$psi, x$shift_speed)
combos[, combo := lab_of(combos)]
setorder(combos, reg, psi, shift_speed)
lev <- combos$combo
add_lab <- function(x) { x[, combo := factor(lab_of(x), levels = lev)]; x }
L <- add_lab(L); O <- add_lab(O)

# response ratio: dg_obs / speed, smoothed over 50 ticks within each sim
for (x in list(L, O)) {
  setorder(x, ID, tick)
  x[, response := frollmean(dg_obs, 50, align = "center") / shift_speed, by = ID]
}

qs <- list(
  list("z_opt",        "lag (optimum - mean phenotype)"),
  list("response",     "response ratio (dg_obs / speed, 50-tick mean)"),
  list("fitness_mean", "mean fitness"),
  list("z_sd",         "SD of phenotype"))

base_theme <- list(theme_bw(base_size = 10),
                   theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92")))
reg_cols <- setNames(rule_cols, names(reg_names))   # keyed by the data codes

pdf(file.path(out_dir, "long_borderline.pdf"), width = 15, height = 9)
for (q in qs) {
  v <- q[[1]]
  reps <- L[tick %% tick_step == 0, .(combo, reg, ID, tick, y = get(v))]
  mL <- L[, .(y = mean(get(v), na.rm = TRUE)), by = .(combo, reg, tick)]
  mO <- O[, .(y = mean(get(v), na.rm = TRUE)), by = .(combo, tick)]
  p <- ggplot() +
    geom_vline(xintercept = c(shift_tick, 12000), linetype = "dotted", colour = "grey40") +
    geom_line(data = reps, aes(tick, y, group = ID, colour = reg), linewidth = 0.15, alpha = 0.35) +
    geom_line(data = mL, aes(tick, y, colour = reg), linewidth = 0.8) +
    geom_line(data = mO, aes(tick, y), colour = "grey20", linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~ combo, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = reg_cols, guide = "none") +
    labs(x = "tick", y = q[[2]], title = sprintf("long_borderline to tick 15000: %s", q[[2]]),
         subtitle = "thin: each of 20 reps; thick: their mean; dashed: mean of the 60 steady_tracking reps (to 12000)") +
    base_theme
  if (v == "response") p <- p + geom_hline(yintercept = 1, colour = "black", linewidth = 0.3) + coord_cartesian(ylim = c(0, 1.5))
  print(p)
}
invisible(dev.off())

# ---- trend per 1000-tick window ----
L[, win := (tick - shift_tick - 1) %/% 1000]
w <- L[tick > shift_tick + 1000, .(lag_first = mean(z_opt[tick <= min(tick) + 99]), lag_last = mean(z_opt[tick >= max(tick) - 99]),
                                    lag = mean(z_opt)),
       by = .(combo, reg, psi, shift_speed, ID, win)]
w[, gain_rel := (lag_last - lag_first) / 900 / shift_speed]
ws <- w[, .(lag = mean(lag), gain_rel = mean(gain_rel), gain_se = sd(gain_rel) / sqrt(.N)),
        by = .(combo, reg, psi, shift_speed, win)]
ws[, window := sprintf("%d-%d", shift_tick + 1000 * win, shift_tick + 1000 * (win + 1))]
fwrite(ws[order(combo, win)], file.path(out_dir, "long_borderline_windows.csv"))

pdf(file.path(out_dir, "long_borderline_trend.pdf"), width = 12, height = 6)
for (pp in list(list("gain_rel", "gain_se", "lag gain per tick / speed (0 = steady lag, 1 = mean not moving)"))) {
  print(ggplot(ws, aes(window, get(pp[[1]]), colour = reg, group = combo)) +
          geom_hline(yintercept = 0, colour = "grey40") +
          geom_pointrange(aes(ymin = get(pp[[1]]) - 2 * get(pp[[2]]), ymax = get(pp[[1]]) + 2 * get(pp[[2]])), size = 0.25) +
          geom_line(linewidth = 0.4) +
          facet_wrap(~ combo, ncol = 4) +
          scale_colour_manual(values = reg_cols, guide = "none") +
          labs(x = "ticks (1000-tick windows)", y = pp[[3]], title = pp[[3]],
               subtitle = "mean of 20 reps, +/- 2 SE; first and last 100 ticks of each window compared") +
          base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}
invisible(dev.off())

print(dcast(ws, combo ~ window, value.var = "gain_rel", fun.aggregate = function(x) round(x, 3)))
print(dcast(ws, combo ~ window, value.var = "lag", fun.aggregate = function(x) round(x)))
cat("plots written to", out_dir, "\n")
