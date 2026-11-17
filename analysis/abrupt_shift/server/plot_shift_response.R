#!/usr/bin/env Rscript
# Response to an abrupt shift of the optimum (cons_env_data_with_adjusted_shift_size, n = 8,
# optimum jumps from 0 to shift_size at tick 10000). Columns = copying rule.
#   row 1: mean phenotype z_mean over time for psi_show (dots = mean over reps, ribbon = 95% CI);
#          dashed: the optimum
#   row 2: phenotype distribution over time at psi_heat, all reps pooled (heatmap: fraction
#          of individuals per bin); red bar = mean of each recorded generation's distribution,
#          spanning that generation's column of the heatmap (not connected)
#   row 3: same for the genotype
#
# Input (analysis/load_data/server/load_data.R -d cons_env_data_with_adjusted_shift_size --n 8
#        --shift-size 28.2 ...):
#   data_summary_n8_sz28.2_gen9950-10300.rds          (-t data --psi 0,0.2,0.5,0.8
#                                                       --gens 9950:10300 --cols z_mean,z_opt,g_mean)
#   {phenotype,genotype}_n8_sz28.2_psi0.8_gen9950-10300.rds   (-t phenotype,genotype --psi 0.8
#                                                       --gens 9950:10300)
#   in figure/r_data/cons_env_data_with_adjusted_shift_size/. Per-individual data are recorded
#   at 9950-9995 every 5, 9996-10010 every tick, then every 5.
# Output: analysis/abrupt_shift/server/plots/shift_response/shift_response_n8_sz28.2.pdf
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/abrupt_shift/server/plot_shift_response.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(MetBrewer)
  library(patchwork)
})

# ---- settings ----
data_dir    <- "figure/r_data/cons_env_data_with_adjusted_shift_size"
data_tag    <- "n8_sz28.2_gen9950-10300"
ind_tag     <- "n8_sz28.2_psi0.8_gen9950-10300"
out_dir     <- "analysis/abrupt_shift/server/plots/shift_response"
shift_sz    <- 28.2   # optimum jump (the data column is shift_size)
shift_tick  <- 10000          # x = tick - shift_tick; the optimum is shift_size from x = 0
psi_show    <- c(0, 0.2, 0.5, 0.8)
psi_grid    <- seq(0, 1, by = 0.1)   # psi colour bar, as in plot_equilibrium_dist.R
psi_heat    <- 0.8
x_lim_mean  <- c(-10, 100)    # row 1 (pre-shift mostly left out)
y_lim_mean  <- c(-4, 30)
x_lim_heat  <- c(-50, 300)    # rows 2-3
bin_width   <- 1
y_lim_heat  <- c(-25.5, 55.5) # bin edges of the heatmaps
opt_colour  <- "#2F5BEA"
mean_colour <- "#E41A1C"
# ---- figure size (same as plot_equilibrium_dist.R): single journal column, 8.7 cm wide ----
fig_w     <- 8.7 / 2.54   # inches (3.42 in); reference width for the heights below
# every figure is also saved at these widths (cm): single column, 1.5 column, full page.
# Height stays that of the 8.7 cm version (same row heights): wider figures only stretch the
# panels horizontally. Text and line sizes stay fixed.
fig_widths_cm <- c(8.7)
save_widths <- function(plot, stem, height_at_ref) {
  for (wcm in fig_widths_cm) {
    w <- wcm / 2.54
    f <- file.path(out_dir, sprintf("%s_%scm.pdf", stem, format(wcm)))
    ggsave(f, plot, width = w, height = height_at_ref)
    cat("written", f, "\n")
  }
}
tick_pt   <- 5            # axis tick labels
txt_pt    <- 5            # axis titles, legend labels, in-panel text
strip_pt  <- 5.5          # copying-rule column titles (one line)
sym_pt    <- 6.5          # single-symbol titles: psi / n colour bars, z-bar / z / g axes
gap       <- unit(4, "pt")  # horizontal space between panels
panel_ratio <- 11 / 16    # every panel is 16:11 (width:height)
curve_lwd <- 0.6          # distribution curves
mean_lwd  <- 0.55         # mean lines / mean bars
point_size <- 0.7         # dots
ref_lwd   <- 0.3          # dashed reference lines
err_lwd   <- 0.25         # error bars
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
reg_names   <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
reg_factor <- function(r) factor(reg_names[as.character(r)], levels = reg_names)
optimum_line <- function(xl) data.table(x = c(xl[1], 0, 0, xl[2]),
                                        y = c(0, 0, shift_sz, shift_sz), seg = c(1, 1, 2, 2))

base_theme <- list(
  theme_bw(base_size = txt_pt),
  fig_theme,
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size = sym_pt, angle = 0, vjust = 0.5, margin = margin(r = 1)))
)

# ---- row 1: mean phenotype ----
s <- as.data.table(readRDS(file.path(data_dir, sprintf("data_summary_%s.rds", data_tag))))
s <- s[abs(shift_size - shift_sz) < 1e-9 & psi %in% psi_show]
s[, x := tick - shift_tick]
s <- s[x >= x_lim_mean[1] & x <= x_lim_mean[2]]
zm <- s[, .(mean = mean(z_mean), se = sd(z_mean) / sqrt(.N), reps = .N), by = .(reg, psi, x)]
zm[, reg := reg_factor(reg)]
cat(sprintf("row 1: %d-%d reps per point\n", min(zm$reps), max(zm$reps)))

pal <- as.character(met.brewer("Hokusai2", length(psi_grid)))
p1 <- ggplot(zm[order(psi)], aes(x, mean, colour = psi, fill = psi, group = psi)) +
  geom_line(data = optimum_line(x_lim_mean), aes(x, y, group = seg),
            colour = opt_colour, linetype = "dashed", linewidth = ref_lwd, inherit.aes = FALSE) +
  # label the dashed line in the first panel instead of a legend key
  geom_text(data = data.table(reg = factor(reg_names[1], levels = reg_names)),
            aes(x = -9, y = shift_sz - 2.5, label = "optimum"), colour = opt_colour, hjust = 0,
            size = (txt_pt - 1) / .pt, inherit.aes = FALSE) +
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
  scale_x_continuous(limits = x_lim_mean, expand = expansion(0.01), breaks = c(0, 50, 100),
                     labels = c("0", "50", "")) +
  # y axis labels the new optimum
  scale_y_continuous(breaks = c(0, 10, 20, shift_sz), labels = as.character(c(0, 10, 20, shift_sz)),
                     expand = expansion(0)) +
  coord_cartesian(ylim = y_lim_mean) +
  labs(x = NULL, y = expression(bar(z))) +
  base_theme +
  # psi colour bar in the bottom-right corner of the Copy-the-extreme panel
  theme(legend.position = "inside", legend.position.inside = c(0.995, 0.02),
        legend.justification.inside = c(1, 0))

# ---- rows 2-3: distributions over time at psi_heat ----
read_type <- function(type) {
  x <- as.data.table(readRDS(file.path(data_dir, sprintf("%s_%s.rds", type, ind_tag))))
  x <- x[abs(psi - psi_heat) < 1e-9]
  setnames(x, type, "value")
  x[, .(reg, rep, x = generation - shift_tick, value)][x >= x_lim_heat[1] & x <= x_lim_heat[2]]
}
breaks <- seq(y_lim_heat[1], y_lim_heat[2], by = bin_width)
heat_data <- function(type) {
  d <- read_type(type)
  d[, bin := findInterval(value, breaks, rightmost.closed = TRUE)]
  d[, total := .N, by = .(reg, x)]
  h <- d[bin >= 1 & bin < length(breaks), .(frac = .N / total[1]), by = .(reg, x, bin)]
  # generations are recorded every 1 or 5 ticks: fill each tick with the nearest recorded
  # generation so the heatmap is a regular grid (drawn with geom_raster, no seams)
  gens <- sort(unique(d$x))
  ticks <- data.table(tick = seq(min(gens), max(gens)))
  ticks[, x := gens[findInterval(tick, (head(gens, -1) + tail(gens, -1)) / 2) + 1]]
  h <- h[ticks, on = "x", allow.cartesian = TRUE, nomatch = NULL]
  h[, `:=`(y = (breaks[bin] + breaks[bin + 1]) / 2, reg = reg_factor(reg))]
  m <- d[, .(mean = mean(value)), by = .(reg, x)]
  # the ticks each recorded generation covers in the raster -> extent of its mean bar
  span <- ticks[, .(xmin = min(tick) - 0.5, xmax = max(tick) + 0.5), by = x]
  m <- span[m, on = "x"]
  m[, reg := reg_factor(reg)]
  cat(sprintf("%s: %d reps pooled per generation; %.4f%% outside the y range\n", type,
              uniqueN(d$rep), 100 * d[bin < 1 | bin >= length(breaks), .N] / nrow(d)))
  list(h = h, m = m, reps = uniqueN(d$rep))
}
hz <- heat_data("phenotype")
hg <- heat_data("genotype")
fill_lim <- c(0, max(hz$h$frac, hg$h$frac))

heat_plot <- function(hd, ylab, show_fill, show_mean) {
  ggplot() +
    geom_raster(data = hd$h, aes(tick, y, fill = frac)) +
    geom_line(data = optimum_line(x_lim_heat), aes(x, y, group = seg), colour = opt_colour,
              linetype = "dashed", linewidth = ref_lwd) +
    geom_segment(data = hd$m, aes(x = xmin, xend = xmax, y = mean, yend = mean,
                                 colour = "Distribution\nmean"), linewidth = mean_lwd) +
    annotate("text", x = -Inf, y = Inf, label = sprintf("psi == %s", psi_heat), parse = TRUE,
             hjust = -0.15, vjust = 1.4, size = txt_pt / .pt) +
    facet_wrap(~ reg, nrow = 1) +
    scale_fill_gradient("Density", low = "white", high = "black", limits = fill_lim,
                        guide = if (show_fill) guide_colourbar(order = 2, direction = "vertical",
                                                                 title.position = "left",
                                                                 theme = theme(legend.title = element_text(angle = 90, vjust = 0.5, hjust = 0.5,
                                                                                                           size = txt_pt - 1)),
                                                                 barheight = unit(0.55, "cm"),
                                                                 barwidth = unit(0.11, "cm")) else "none") +
    scale_colour_manual(NULL, values = c("Distribution\nmean" = mean_colour),
                        guide = if (show_mean) guide_legend(order = 1, keywidth = unit(0.25, "cm"),
                                                              keyheight = unit(0.2, "cm")) else "none") +
    scale_x_continuous(limits = x_lim_heat + c(-0.5, 0.5), expand = expansion(0), breaks = c(0, 100, 200, 300),
                       labels = c("0", "", "200", "")) +
    scale_y_continuous(limits = y_lim_heat, expand = expansion(0)) +
    labs(x = NULL, y = ylab) +
    base_theme +
    theme(strip.text = element_blank(),
          legend.position = "inside", legend.position.inside = c(0.995, 0.02),
          legend.justification.inside = c(1, 0),
          legend.box = "horizontal", legend.box.just = "bottom", legend.spacing.x = unit(1, "pt"))
}
p2 <- heat_plot(hz, "z", show_fill = FALSE, show_mean = TRUE) +
  theme(axis.text.x = element_blank())   # same x axis as row 3
p3 <- heat_plot(hg, "g", show_fill = TRUE, show_mean = FALSE) + labs(x = "Generation after optimum shift")

p_all <- (p1 / p2 / p3) + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
save_widths(p_all, sprintf("shift_response_n8_sz%s", format(shift_sz)), 2.61)
cat(sprintf("n = 8, shift size %s. Row 1: mean over %d reps, ribbon = 95%% CI. Rows 2-3: psi = %s, %d reps pooled per generation\n",
            format(shift_sz), max(zm$reps), psi_heat, hz$reps))
