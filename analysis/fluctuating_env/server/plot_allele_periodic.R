#!/usr/bin/env Rscript
# Heatmaps of the per-sim allele summaries from allele_periodic_summary.R (mean over reps), in the
# layout of plot_fitness_heatmap.R: rows = psi, columns = copying rule, x = amplitude (fluc_opt,
# in multiples of sd_unit), y = period of the fluctuation (2 * fluc_interval).
# One PDF per measure (definitions in allele_periodic_summary.R):
#   swing         mean periodic swing of an allele's frequency (alleles reaching 0.05), log colour
#   coherence     share of allele-frequency variance at the forcing period, log colour
#   gain          periodic swing of the genetic mean / that of the optimum
#   gain_inphase  same, only the part in phase with the optimum
#   lag           delay of the genetic mean's swing behind the optimum, in generations
#   lag_frac      the same as a fraction of the period
#
# Input:  figure/r_data/fluc_env_adjusted_shift_size/allele_periodic_n8_rep1_10.rds
# Output: analysis/fluctuating_env/server/plots/allele_periodic/allele_<measure>_n8_8.7cm.pdf
#         and allele_periodic_n8.csv (tile means)
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/plot_allele_periodic.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- settings ----
in_file   <- "figure/r_data/fluc_env_adjusted_shift_size/allele_periodic_n8_rep1_10.rds"
out_dir   <- "analysis/fluctuating_env/server/plots/allele_periodic"
n_tag     <- "n8"
sd_unit   <- 9.4            # standing genetic sd at psi = 0 under stabilizing selection
reg_names <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")
measures <- data.table(
  name  = c("swing", "coherence", "gain", "gain_inphase", "lag", "lag_frac"),
  title = c("Periodic swing of allele frequency",
            "Coherence: share of allele-frequency variance at the forcing period",
            "Genetic tracking gain: swing of genetic mean / swing of optimum",
            "In-phase genetic tracking gain",
            "Lag of the genetic mean behind the optimum (generations)",
            "Lag of the genetic mean behind the optimum (fraction of period)"),
  log   = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))

# ---- figure size (as plot_fitness_heatmap.R) ----
fig_w    <- 8.7 / 2.54
fig_h    <- 3.5
tick_pt  <- 4.5
txt_pt   <- 5
strip_pt <- 5.5
sym_pt   <- 6

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
s <- as.data.table(readRDS(in_file))
s[, lag_frac := lag / period]
m <- s[, c(lapply(.SD, mean), list(reps = .N)), by = .(reg, fluc_opt, period, psi), .SDcols = measures$name]
fwrite(m[order(reg, psi, period, fluc_opt)], file.path(out_dir, sprintf("allele_periodic_%s.csv", n_tag)))
cat(sprintf("%d sims, %d-%d reps per tile\n", nrow(s), min(m$reps), max(m$reps)))

psi_show <- sort(unique(m$psi))
k <- m$fluc_opt / sd_unit
stopifnot(all(abs(k - round(k)) < 1e-9))
m[, `:=`(reg = factor(reg_names[reg], levels = reg_names),
         psi_lab = factor(sprintf("psi == %s", format(psi)), levels = sprintf("psi == %s", format(psi_show))),
         x = factor(as.integer(round(fluc_opt / sd_unit))),
         y = factor(period, levels = sort(unique(period))))]

for (i in seq_len(nrow(measures))) {
  ms <- measures[i]
  d <- copy(m)[, val := get(ms$name)]
  fill_scale <- scale_fill_viridis_c(NULL, option = "magma", direction = -1,
                                     transform = if (ms$log) "log10" else "identity",
                                     n.breaks = 5,
                                     guide = guide_colourbar(direction = "horizontal",
                                                             barwidth = unit(4, "cm"),
                                                             barheight = unit(0.13, "cm")))
  p <- ggplot(d, aes(x, y, fill = val)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_grid(psi_lab ~ reg, switch = "y",
               labeller = labeller(psi_lab = label_parsed, reg = label_value)) +
    fill_scale +
    scale_x_discrete(expand = expansion(0)) +
    scale_y_discrete(expand = expansion(0)) +
    labs(x = bquote("Amplitude (distance between optima, in" ~ sigma[G] == .(sd_unit) * ")"),
         y = "Period of fluctuation (generations)", title = ms$title) +
    theme_bw(base_size = txt_pt) +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(3, "pt"),
          plot.title = element_text(size = txt_pt, hjust = 0.5),
          axis.text = element_text(size = tick_pt),
          axis.ticks.length = unit(1.2, "pt"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(size = strip_pt, margin = margin(1, 0, 1.5, 0)),
          strip.text.y.left = element_text(size = sym_pt, angle = 90, margin = margin(0, 1, 0, 0)),
          legend.position = "bottom",
          legend.text = element_text(size = tick_pt),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "pt"),
          plot.margin = margin(2, 4, 2, 2))
  f <- file.path(out_dir, sprintf("allele_%s_%s_8.7cm.pdf", ms$name, n_tag))
  ggsave(f, p, width = fig_w, height = fig_h)
  cat(sprintf("  %-12s range %.4g - %.4g; written %s\n", ms$name, min(d$val), max(d$val), f))
}
