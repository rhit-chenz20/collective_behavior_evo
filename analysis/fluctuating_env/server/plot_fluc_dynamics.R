#!/usr/bin/env Rscript
# Fluctuating optimum (fluc_env_adjusted_shift_size, n = 8): the optimum alternates between
# +fluc_opt/2 and -fluc_opt/2, flipping every fluc_interval ticks (at ticks 50 + k * interval).
# One figure per corner of the amplitude x frequency grid:
#   small / large amplitude  = fluc_opt 9.4 / 103.4
#   frequent / infrequent    = fluc_interval 20 / 110
# Each figure shows a window of two full cycles (4 x fluc_interval ticks) starting at the first
# flip at or after tick win_from. Columns = copying rule.
#   row 1: phenotype distribution over time at psi_heat, all reps pooled (heatmap: fraction of
#          individuals per bin); red bar = mean of each recorded generation, spanning its column
#   row 2: same for the genotype
#   rows 1-2: dashed = the optimum; density bar and mean key in a strip under row 2
#   row 3: population mean fitness for psi_show (dots = mean over reps, ribbon = 95% CI),
#          relative to the maximum dnorm(0, 0, sqrt(Vs)); dotted verticals = optimum flips
#
# Input (analysis/load_data/server/load_data.R -d fluc_env_adjusted_shift_size --n 8
#        --fluc-opt 9.4,103.4 --fluc-interval 20,110 --gens 10000:10510 ...):
#   data_summary_n8_corners_gen10000-10510.rds   (-t data --psi 0,0.2,0.5,0.8
#                                                 --cols optimum,z_mean,g_mean,fitness_mean)
#   {phenotype,genotype}_n8_corners_psi0.8_gen10000-10510.rds   (-t phenotype,genotype --psi 0.8)
#   in figure/r_data/fluc_env_adjusted_shift_size/. Per-individual data are recorded every
#   tick at 9995-10009, then every 5 ticks.
# Output: analysis/fluctuating_env/server/plots/fluc_dynamics/
#         fluc_dynamics_n8_<small|large>_<frequent|infrequent>_opt<..>_int<..>_8.7cm.pdf
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/plot_fluc_dynamics.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(MetBrewer)
  library(patchwork)
})

# ---- settings ----
data_dir    <- "figure/r_data/fluc_env_adjusted_shift_size"
data_tag    <- "n8_corners_gen10000-10510"
ind_tag     <- "n8_corners_psi0.8_gen10000-10510"
out_dir     <- "analysis/fluctuating_env/server/plots/fluc_dynamics"
corners <- data.table(fluc_opt      = c(9.4, 9.4, 103.4, 103.4),
                      fluc_interval = c(20, 110, 20, 110),
                      name = c("small_frequent", "small_infrequent", "large_frequent", "large_infrequent"))
win_from    <- 10010          # window starts at the first flip at or after this tick
n_cycles    <- 2              # full cycles shown (one cycle = two flips)
psi_show    <- c(0, 0.2, 0.5, 0.8)
psi_grid    <- seq(0, 1, by = 0.1)
psi_heat    <- 0.8
heat_margin <- 35             # heatmap y range: +-(fluc_opt / 2 + heat_margin)
bin_width   <- 1
Vs          <- 2000           # fitness = dnorm(z, optimum, sqrt(Vs)) (fluc_*.slim)
opt_colour  <- "#2F5BEA"
mean_colour <- "#E41A1C"
reg_names   <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")

# ---- figure size (same as analysis/abrupt_shift/server/plot_shift_response.R) ----
fig_w       <- 8.7 / 2.54
fig_h       <- 2.72
tick_pt     <- 5
txt_pt      <- 5
strip_pt    <- 5.5
sym_pt      <- 6.5
gap         <- unit(4, "pt")
panel_ratio <- 11 / 16        # every panel is 16:11 (width:height)
mean_lwd    <- 0.55
point_size  <- 0.7
ref_lwd     <- 0.3
fig_theme <- theme(plot.margin = margin(2, 4, 2, 0), legend.background = element_blank(),
                   legend.key = element_blank(), legend.margin = margin(0, 1, 1, 0),
                   legend.text = element_text(size = txt_pt - 0.5, margin = margin(t = 0.5)),
                   axis.text = element_text(size = tick_pt),
                   axis.title.y = element_text(margin = margin(r = 1)),
                   axis.ticks.length = unit(1.2, "pt"),
                   aspect.ratio = panel_ratio,
                   panel.spacing.x = gap,
                   strip.text = element_text(size = strip_pt, margin = margin(1, 0, 1.5, 0)),
                   strip.background = element_blank())
base_theme <- list(
  theme_bw(base_size = txt_pt),
  fig_theme,
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size = sym_pt, angle = 0, vjust = 0.5, margin = margin(r = 1)))
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
reg_factor <- function(r) factor(reg_names[as.character(r)], levels = reg_names)
pal <- as.character(met.brewer("Hokusai2", length(psi_grid)))
fmax <- dnorm(0, 0, sqrt(Vs))

s_all <- as.data.table(readRDS(file.path(data_dir, sprintf("data_summary_%s.rds", data_tag))))
read_ind <- function(type) {
  x <- as.data.table(readRDS(file.path(data_dir, sprintf("%s_%s.rds", type, ind_tag))))
  setnames(x, type, "value")
  x[abs(psi - psi_heat) < 1e-9, .(fluc_opt, fluc_interval, reg, rep, generation, value)]
}
ind_all <- list(phenotype = read_ind("phenotype"), genotype = read_ind("genotype"))

for (k in seq_len(nrow(corners))) {
  cn <- corners[k]
  sel <- function(x) x[abs(fluc_opt - cn$fluc_opt) < 1e-9 & fluc_interval == cn$fluc_interval]
  s <- sel(s_all)

  # window: two full cycles from the first flip at or after win_from
  opt1 <- s[ID == ID[1]][order(tick)]
  flips <- opt1$tick[which(diff(opt1$optimum) != 0) + 1]
  t0 <- min(flips[flips >= win_from])
  t1 <- t0 + 2 * n_cycles * cn$fluc_interval - 1   # last tick before the next flip
  stopifnot(t1 <= max(s$tick))
  flips_in <- flips[flips > t0 & flips < t1] - t0
  # optimum as separate horizontal segments (no vertical connectors)
  opt <- opt1[tick >= t0 & tick <= t1, .(x = tick - t0, y = optimum)]
  opt[, seg := rleid(y)]
  opt <- opt[, .(x = c(min(x) - 0.5, max(x) + 0.5), y = y[1]), by = seg]
  x_lim <- c(-0.5, t1 - t0 + 0.5)
  y_lim_heat <- c(-1, 1) * (ceiling(cn$fluc_opt / 2 + heat_margin) + 0.5)
  breaks <- seq(y_lim_heat[1], y_lim_heat[2], by = bin_width)

  # ---- rows 1-2: heatmaps ----
  heat_data <- function(type) {
    d <- sel(ind_all[[type]])[generation >= t0 & generation <= t1]
    d[, x := generation - t0]
    d[, bin := findInterval(value, breaks, rightmost.closed = TRUE)]
    d[, total := .N, by = .(reg, x)]
    h <- d[bin >= 1 & bin < length(breaks), .(frac = .N / total[1]), by = .(reg, x, bin)]
    # fill each tick with the nearest recorded generation (regular grid for geom_raster)
    gens <- sort(unique(d$x))
    ticks <- data.table(tick = seq(0, t1 - t0))   # whole window, incl. ticks after the last record
    ticks[, x := gens[findInterval(tick, (head(gens, -1) + tail(gens, -1)) / 2) + 1]]
    h <- h[ticks, on = "x", allow.cartesian = TRUE, nomatch = NULL]
    h[, `:=`(y = (breaks[bin] + breaks[bin + 1]) / 2, reg = reg_factor(reg))]
    m <- d[, .(mean = mean(value)), by = .(reg, x)]
    span <- ticks[, .(xmin = min(tick) - 0.5, xmax = max(tick) + 0.5), by = x]
    m <- span[m, on = "x"]
    m[, reg := reg_factor(reg)]
    cat(sprintf("  %s: %d reps pooled, %d recorded generations, %.3f%% outside the y range\n", type,
                uniqueN(d$rep), length(gens), 100 * d[bin < 1 | bin >= length(breaks), .N] / nrow(d)))
    list(h = h, m = m)
  }
  cat(sprintf("%s: fluc_opt %s, interval %d, window %d-%d\n", cn$name, format(cn$fluc_opt),
              cn$fluc_interval, t0, t1))
  hz <- heat_data("phenotype")
  hg <- heat_data("genotype")
  fill_lim <- c(0, max(hz$h$frac, hg$h$frac))

  heat_plot <- function(hd, ylab, show_fill, show_mean, top) {
    g <- ggplot() +
      geom_raster(data = hd$h, aes(tick, y, fill = frac)) +
      geom_line(data = opt, aes(x, y, group = seg), colour = opt_colour,
                linetype = "dashed", linewidth = ref_lwd) +
      geom_segment(data = hd$m, aes(x = xmin, xend = xmax, y = mean, yend = mean,
                                   colour = "Distribution mean"), linewidth = mean_lwd) +
      annotate("text", x = -Inf, y = Inf, label = sprintf("psi == %s", psi_heat), parse = TRUE,
               hjust = -0.15, vjust = 1.4, size = txt_pt / .pt) +
      facet_wrap(~ reg, nrow = 1) +
      scale_fill_gradient("Density", low = "white", high = "black", limits = fill_lim,
                          breaks = c(0, floor(fill_lim[2] / 0.05) * 0.05), labels = function(x) format(x, drop0trailing = TRUE),
                          guide = if (show_fill) guide_colourbar(order = 2, direction = "horizontal", title.position = "left",
                                                                 theme = theme(legend.title = element_text(size = txt_pt, vjust = 0.5,
                                                                                                           margin = margin(r = 3))),
                                                                 barwidth = unit(1.25, "cm"),
                                                                 barheight = unit(0.11, "cm")) else "none") +
      scale_colour_manual(NULL, values = c("Distribution mean" = mean_colour),
                          guide = if (show_mean) guide_legend(order = 1, keywidth = unit(0.25, "cm"),
                                                              keyheight = unit(0.2, "cm")) else "none") +
      scale_x_continuous(limits = x_lim, expand = expansion(0), n.breaks = 3) +
      scale_y_continuous(limits = y_lim_heat, expand = expansion(0), n.breaks = 4) +
      labs(x = NULL, y = ylab) +
      base_theme +
      theme(axis.text.x = element_blank(),
            legend.position = "bottom", legend.justification = "right",
            legend.box = "horizontal", legend.spacing.x = unit(6, "pt"),
            legend.box.spacing = unit(1, "pt"), legend.margin = margin(0, 0, 0, 0))
    if (!top) g <- g + theme(strip.text = element_blank())
    g
  }
  p1 <- heat_plot(hz, "z", show_fill = FALSE, show_mean = FALSE, top = TRUE)
  p2 <- heat_plot(hg, "g", show_fill = TRUE, show_mean = TRUE, top = FALSE)

  # ---- row 3: mean fitness ----
  w <- s[psi %in% psi_show & tick >= t0 & tick <= t1]
  w[, `:=`(x = tick - t0, rel = fitness_mean / fmax)]
  wm <- w[, .(mean = mean(rel), se = sd(rel) / sqrt(.N), reps = .N), by = .(reg, psi, x)]
  wm[, reg := reg_factor(reg)]
  p3 <- ggplot(wm[order(psi)], aes(x, mean, colour = psi, fill = psi, group = psi)) +
    geom_vline(xintercept = flips_in, linetype = "dotted", colour = "grey60", linewidth = ref_lwd) +
    geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), colour = NA, alpha = 0.3) +
    geom_point(size = point_size, shape = 16) +
    facet_wrap(~ reg, nrow = 1) +
    scale_colour_gradientn(expression(psi), colours = pal, values = psi_grid, limits = range(psi_grid),
                           breaks = psi_show, labels = format(psi_show, drop0trailing = TRUE),
                           aesthetics = c("colour", "fill"),
                           guide = guide_colourbar(nbin = length(psi_grid), display = "rectangles",
                                                   direction = "horizontal", title.position = "left",
                                                   theme = theme(legend.title = element_text(size = sym_pt, vjust = 0.75,
                                                                                             margin = margin(r = 2))),
                                                   barwidth = unit(1.25, "cm"), barheight = unit(0.11, "cm"))) +
    scale_x_continuous(limits = x_lim, expand = expansion(0), n.breaks = 3) +
    scale_y_continuous(n.breaks = 4, expand = expansion(mult = c(0.45, 0.05))) +   # room for the psi bar
    labs(x = "Generation since window start", y = expression(bar(w) / w[max])) +
    base_theme +
    theme(strip.text = element_blank(),
          legend.position = "inside", legend.position.inside = c(0.995, 0.02),
          legend.justification.inside = c(1, 0))

  p_all <- (p1 / p2 / p3) + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
  f <- file.path(out_dir, sprintf("fluc_dynamics_n8_%s_opt%s_int%d_8.7cm.pdf", cn$name,
                                  format(cn$fluc_opt), cn$fluc_interval))
  ggsave(f, p_all, width = fig_w, height = fig_h)
  cat(sprintf("  row 3: %d reps per point; written %s\n", max(wm$reps), f))
}
