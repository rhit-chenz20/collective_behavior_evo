#!/usr/bin/env Rscript
# Fluctuating optimum (fluc_env_adjusted_shift_size, n = 8): genetic sd of psi = x relative to
# psi = 0, on the amplitude x period grid. Same layout, windows and label versions as
# plot_fitness_heatmap.R, but the per-sim value is the MEAN of g_sd over the window (not a sum).
#   Tile value = mean over reps of mean g_sd(psi = x) - mean over reps of mean g_sd(psi = 0), same
#   (reg, fluc_opt, fluc_interval). Rows = psi_show, columns = copying rule.
# Windows (one figure each), both starting at the first optimum flip at or after win_from:
#   cycles : n_cycles full cycles (4 * fluc_interval ticks: 80, 200, 320, 440)
#   fixed  : fixed_len ticks for every interval
#   x labels: value = fluc_opt (9.4 ... 103.4); sd = fluc_opt in units of sd_unit (1, 3, ..., 11)
#   y labels: period = 2 * fluc_interval (40 ... 220); freq = 1 / period (1/40 ... 1/220)
#   files: gsd_heatmap_n8_<window>_<value|sd>_<period|freq>_8.7cm.pdf and a .csv per window
#
# Input:  figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_track_gen9950-11000.rds
#         (load_data.R -d fluc_env_adjusted_shift_size -t data --n 8 --psi 0,0.2,0.5,0.8
#          --cols optimum,z_mean,g_mean,z_sd,g_sd,H,seg_site_pop --gens 9950:11000
#          --tag n8_track_gen9950-11000)
# Output: analysis/fluctuating_env/server/plots/gsd_heatmap/
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/plot_gsd_heatmap.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- settings ----
data_file <- "figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_track_gen9950-11000.rds"
out_dir   <- "analysis/fluctuating_env/server/plots/gsd_heatmap"
n_tag     <- "n8"
psi_show  <- c(0.2, 0.5, 0.8)
win_from  <- 10000          # window starts at the first flip at or after this tick
n_cycles  <- 2
fixed_len <- 300
reg_names <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")
low_col   <- "#2B4BE0"
high_col  <- "#E0402B"
sd_unit   <- 9.4            # standing genetic sd at psi = 0 under stabilizing selection

# ---- figure size ----
fig_w    <- 8.7 / 2.54
fig_h    <- 2.85
tick_pt  <- 4.5
txt_pt   <- 5
strip_pt <- 5.5
sym_pt   <- 6

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
s <- as.data.table(readRDS(data_file))

# first flip at or after win_from, per (fluc_opt, fluc_interval), from one sim's logged optimum
flip0 <- s[, .(ID1 = ID[1]), by = .(fluc_opt, fluc_interval)][
  s, on = .(fluc_opt, fluc_interval)][ID == ID1][order(tick)][
  , .(t0 = { f <- tick[which(diff(optimum) != 0) + 1]; min(f[f >= win_from]) }),
  by = .(fluc_opt, fluc_interval)]

make_tiles <- function(win) {
  w <- flip0[, .(fluc_opt, fluc_interval, t0,
                 len = if (win == "cycles") 2 * n_cycles * fluc_interval else fixed_len)]
  w[, t1 := t0 + len - 1]
  stopifnot(all(w$t1 <= max(s$tick)))
  d <- s[w, on = .(fluc_opt, fluc_interval)][tick >= t0 & tick <= t1]
  per_sim <- d[, .(W = mean(g_sd), ticks = .N, len = len[1]),
               by = .(reg, fluc_opt, fluc_interval, psi, rep)]
  short <- per_sim[ticks < len]
  if (nrow(short)) cat(sprintf("  %s: %d sims with missing ticks in the window (dropped)\n", win, nrow(short)))
  per_sim <- per_sim[ticks == len]
  m <- per_sim[, .(W = mean(W), reps = .N, len = len[1]), by = .(reg, fluc_opt, fluc_interval, psi)]
  base <- m[abs(psi) < 1e-9, .(reg, fluc_opt, fluc_interval, W0 = W)]
  m <- base[m, on = .(reg, fluc_opt, fluc_interval)][psi > 1e-9]
  m[, diff := W - W0]
  m[, reps_min := min(reps), by = .(reg, fluc_opt, fluc_interval)]
  cat(sprintf("  %s: window lengths %s ticks; %d-%d reps per tile; diff range %.2f to %.2f\n", win,
              paste(sort(unique(m$len)), collapse = "/"), min(m$reps), max(m$reps),
              min(m$diff), max(m$diff)))
  m
}

# x tick labels: the float value (one decimal) or the multiple of sd_unit
x_lab <- function(v, style) {
  if (style == "value") return(sprintf("%.1f", v))
  k <- v / sd_unit
  stopifnot(all(abs(k - round(k)) < 1e-9))
  sprintf("%d", as.integer(round(k)))
}

# y: period (full cycle = 2 flips) or frequency; frequency increases upward
y_factor <- function(interval, style) {
  per <- 2 * interval
  lv <- sort(unique(per))
  if (style == "period") return(factor(per, levels = lv))
  factor(sprintf("1/%d", per), levels = sprintf("1/%d", rev(lv)))
}

plot_tiles <- function(m, win, xlab_style, ylab_style) {
  m <- copy(m)
  m[, `:=`(reg = factor(reg_names[reg], levels = reg_names),
           psi_lab = factor(sprintf("psi == %s", format(psi)), levels = sprintf("psi == %s", format(psi_show))),
           x = factor(x_lab(fluc_opt, xlab_style), levels = x_lab(sort(unique(fluc_opt)), xlab_style)),
           y = y_factor(fluc_interval, ylab_style))]
  lim <- max(abs(m$diff))
  sub <- if (win == "cycles") sprintf("%d full cycles = %s gen.", n_cycles,
                                      paste(sort(unique(m$len)), collapse = "/"))
         else sprintf("%d generations", fixed_len)
  ggplot(m, aes(x, y, fill = diff)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_grid(psi_lab ~ reg, switch = "y",
               labeller = labeller(psi_lab = label_parsed, reg = label_value)) +
    scale_fill_gradient2(bquote("Mean genetic sd," ~ psi == x ~ "minus" ~ psi == 0 ~ "(" * .(sub) * ")"),
                         low = low_col, mid = "white", high = high_col, midpoint = 0,
                         limits = c(-lim, lim), n.breaks = 5,
                         guide = guide_colourbar(direction = "horizontal", title.position = "top",
                                                 title.hjust = 0.5, barwidth = unit(4, "cm"),
                                                 barheight = unit(0.13, "cm"))) +
    scale_x_discrete(expand = expansion(0)) +
    scale_y_discrete(expand = expansion(0)) +
    labs(x = if (xlab_style == "value") "Amplitude (distance between optima)"
             else bquote("Amplitude (distance between optima, in" ~ sigma[G] == .(sd_unit) * ")"),
         y = if (ylab_style == "period") "Period of fluctuation (generations)"
             else "Frequency of fluctuation (cycles per generation)") +
    theme_bw(base_size = txt_pt) +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(3, "pt"),
          axis.text = element_text(size = tick_pt),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks.length = unit(1.2, "pt"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(size = strip_pt, margin = margin(1, 0, 1.5, 0)),
          strip.text.y.left = element_text(size = sym_pt, angle = 90, margin = margin(0, 1, 0, 0)),
          legend.position = "bottom",
          legend.title = element_text(size = txt_pt),
          legend.text = element_text(size = tick_pt),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "pt"),
          plot.margin = margin(2, 4, 2, 2))
}

for (win in c("cycles", "fixed")) {
  m <- make_tiles(win)
  tag <- if (win == "cycles") sprintf("%dcycles", n_cycles) else sprintf("%dgen", fixed_len)
  fwrite(m[order(reg, psi, fluc_interval, fluc_opt)],
         file.path(out_dir, sprintf("gsd_heatmap_%s_%s.csv", n_tag, tag)))
  for (xs in c("value", "sd")) for (ys in c("period", "freq")) {
    f <- file.path(out_dir, sprintf("gsd_heatmap_%s_%s_%s_%s_8.7cm.pdf", n_tag, tag, xs, ys))
    ggsave(f, plot_tiles(m, win, xs, ys), width = fig_w, height = fig_h)
    cat(sprintf("  written %s\n", f))
  }
}
