# Shared helpers for the fluctuating-environment heatmaps (layout of plot_fitness_heatmap.R):
# rows = psi, columns = copying rule, x = amplitude (fluc_opt in multiples of sd_unit),
# y = period of the fluctuation (2 * fluc_interval). Sourced by the plot_*.R scripts here.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

sd_unit   <- 9.4            # standing genetic sd at psi = 0 under stabilizing selection
reg_names <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")
reg_factor <- function(r) factor(reg_names[as.character(r)], levels = reg_names)
fig_w <- 8.7 / 2.54
fig_h <- 3.5
hm_tick_pt <- 4.5
hm_txt_pt  <- 5

# whole-number-of-periods window: largest multiple of period that fits in max_len ticks
whole_cycles <- function(period, max_len) floor(max_len / period) * period

# fundamental (complex amplitude) of x(t) at period: x ~ mean + Mod(c) cos(2 pi t / period + Arg(c))
fourier1 <- function(x, t, period) {
  om <- 2 * pi / period
  complex(real = 2 / length(x) * sum(x * cos(om * t)), imaginary = -2 / length(x) * sum(x * sin(om * t)))
}

# tile data: needs reg, psi, fluc_opt, fluc_interval and a column val
hm_prepare <- function(d) {
  d <- copy(d)
  k <- d$fluc_opt / sd_unit
  stopifnot(all(abs(k - round(k)) < 1e-9))
  psis <- sort(unique(d$psi))
  d[, `:=`(reg = reg_factor(reg),
           psi_lab = factor(sprintf("psi == %s", format(psi)), levels = sprintf("psi == %s", format(psis))),
           x = factor(as.integer(round(fluc_opt / sd_unit))),
           y = factor(2L * fluc_interval, levels = sort(unique(2L * fluc_interval))))]
  d
}

# fill: "seq" (magma), "seq_log", or "div" (blue-white-red centred on div_mid)
hm_plot <- function(d, title, fill = "seq", div_mid = 0, legend_title = NULL, limits = NULL,
                    breaks = waiver(), labels = waiver()) {
  d <- hm_prepare(d)
  guide <- guide_colourbar(direction = "horizontal", barwidth = unit(4, "cm"), barheight = unit(0.13, "cm"),
                           title.position = "top", title.hjust = 0.5)
  fs <- switch(fill,
    seq = scale_fill_viridis_c(legend_title, option = "magma", direction = -1, limits = limits,
                               oob = scales::squish,
                               breaks = breaks, labels = labels, guide = guide),
    seq_log = scale_fill_viridis_c(legend_title, option = "magma", direction = -1, transform = "log10",
                                   limits = limits, breaks = breaks, labels = labels, guide = guide),
    div = {
      lim <- if (is.null(limits)) div_mid + c(-1, 1) * max(abs(d$val - div_mid), na.rm = TRUE) else limits
      scale_fill_gradient2(legend_title, low = "#2B4BE0", mid = "white", high = "#E0402B",
                           midpoint = div_mid, limits = lim, breaks = breaks, labels = labels, guide = guide)
    })
  ggplot(d, aes(x, y, fill = val)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_grid(psi_lab ~ reg, switch = "y", labeller = labeller(psi_lab = label_parsed, reg = label_value)) +
    fs +
    scale_x_discrete(expand = expansion(0)) +
    scale_y_discrete(expand = expansion(0)) +
    labs(x = bquote("Amplitude (distance between optima, in" ~ sigma[G] == .(sd_unit) * ")"),
         y = "Period of fluctuation (generations)", title = title) +
    theme_bw(base_size = hm_txt_pt) +
    theme(panel.grid = element_blank(), panel.spacing = unit(3, "pt"),
          plot.title = element_text(size = hm_txt_pt, hjust = 0.5),
          axis.text = element_text(size = hm_tick_pt), axis.ticks.length = unit(1.2, "pt"),
          strip.background = element_blank(), strip.placement = "outside",
          strip.text.x = element_text(size = 5.5, margin = margin(1, 0, 1.5, 0)),
          strip.text.y.left = element_text(size = 6, angle = 90, margin = margin(0, 1, 0, 0)),
          legend.position = "bottom", legend.title = element_text(size = hm_txt_pt),
          legend.text = element_text(size = hm_tick_pt), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "pt"), plot.margin = margin(2, 4, 2, 2))
}
