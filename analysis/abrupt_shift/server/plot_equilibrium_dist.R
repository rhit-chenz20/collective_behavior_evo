#!/usr/bin/env Rscript
# Phenotype (z) and genotype (g) distributions at equilibrium, just before the abrupt shift
# of the optimum (cons_env_data_with_adjusted_shift_size, n = 8): individuals of every recorded
# generation in gen_window and every rep are pooled, so each curve is the distribution
# averaged over those generations. Columns = copying rule, rows = phenotype / genotype,
# colour = psi (MetBrewer Hokusai2 over psi 0-1; only psi_show is drawn).
# Before the shift the optimum is 0 and shift_size has not acted yet, so one shift size suffices.
#
# Input (load_data.R -d cons_env_data_with_adjusted_shift_size -t phenotype,genotype --n 8
#        --shift-size 18.8 --psi 0,0.2,0.5,0.8 --gens 9950:9999 --tag n8_sz18.8_gen9950-9999):
#   figure/r_data/cons_env_data_with_adjusted_shift_size/{phenotype,genotype}_n8_sz18.8_gen9950-9999.rds
#   (recorded generations: 9950-9995 every 5, then 9996-9999; 60 reps)
# Output: analysis/abrupt_shift/server/plots/distributions/equilibrium_dist_n8.pdf
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/abrupt_shift/server/plot_equilibrium_dist.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(MetBrewer)
  library(patchwork)
})

# ---- settings ----
data_dir   <- "figure/r_data/cons_env_data_with_adjusted_shift_size"
tag        <- "n8_sz18.8_gen9950-9999"
out_dir    <- "analysis/abrupt_shift/server/plots/distributions"
gen_window <- c(9950, 9999)
psi_show   <- c(0, 0.2, 0.5, 0.8)
psi_grid   <- seq(0, 1, by = 0.1)       # colour bar: one Hokusai2 colour per psi
bin_width  <- 1
bin_range  <- c(-40.5, 40.5)            # bins centred on the integers -40..40
reg_names  <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")
# genetic variance vs mutation-drift balance (third row of the combined figure)
vg_tag      <- "n2_8_20_50_gen9950-9999"   # data_summary_<vg_tag>.rds, metadata_<vg_n_tag>.rds
vg_meta_tag <- "n2_8_20_50"
vg_n        <- c(2, 8, 20, 50)
vg_shift    <- 18.8   # pre-shift data of the three shift sizes share burn-ins; use one
Ea2         <- 1      # E[a^2] of the effect size a = +-sqrt(Exp(1)) (cons_*.slim)
type_names <- c(phenotype = "Phenotype (z)\nproportion", genotype = "Genotype (g)\nproportion")  # y titles
# ---- figure size: single journal column, 8.7 cm wide ----
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

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_type <- function(type) {
  x <- as.data.table(readRDS(file.path(data_dir, sprintf("%s_%s.rds", type, tag))))
  x <- x[generation >= gen_window[1] & generation <= gen_window[2] & psi %in% psi_show]
  setnames(x, type, "value")
  x[, .(reg, psi, rep, generation, value, type = type)]
}
d <- rbind(read_type("phenotype"), read_type("genotype"))
cat(sprintf("%s individuals; generations %s; %d reps\n", format(nrow(d), big.mark = ","),
            paste(sort(unique(d$generation)), collapse = ","), uniqueN(d$rep)))

# proportion of individuals per bin; every generation and rep has the same population size,
# so pooling equals averaging the per-generation distributions
# (values outside bin_range are left out of the bins but stay in the denominator)
breaks <- seq(bin_range[1], bin_range[2], by = bin_width)
d[, bin := findInterval(value, breaks, rightmost.closed = TRUE)]
d[, total := .N, by = .(type, reg, psi)]
cat(sprintf("%.4f%% of individuals fall outside %s\n",
            100 * d[bin < 1 | bin >= length(breaks), .N] / nrow(d), paste(bin_range, collapse = " to ")))
h <- d[bin >= 1 & bin < length(breaks), .(prop = .N / total[1]), by = .(type, reg, psi, bin)]
# fill empty bins with 0 so the curves reach the baseline
h <- h[CJ(type = unique(h$type), reg = unique(h$reg), psi = unique(h$psi),
          bin = seq_len(length(breaks) - 1)), on = .(type, reg, psi, bin)]
h[is.na(prop), prop := 0]
h[, `:=`(x = breaks[bin] + bin_width / 2,
         reg = factor(reg_names[reg], levels = reg_names))]

pal <- as.character(met.brewer("Hokusai2", length(psi_grid)))

# one plot per row so each has a real y title; patchwork lines the titles up at the far left
dist_row <- function(tp, top) {
  g <- ggplot(h[type == tp][order(psi)], aes(x, prop, colour = psi, group = psi)) +
    geom_line(linewidth = curve_lwd) +
    facet_wrap(~ reg, nrow = 1) +
    scale_colour_gradientn(expression(psi), colours = pal, values = psi_grid, limits = range(psi_grid),
                           breaks = psi_show, labels = format(psi_show, drop0trailing = TRUE),
                           guide = if (!top) "none" else
                             guide_colourbar(nbin = length(psi_grid), display = "rectangles",
                                             direction = "horizontal", title.position = "left",
                                             theme = theme(legend.title = element_text(size = sym_pt, vjust = 0.75,
                                                                                       margin = margin(r = 2))),
                                             barwidth = unit(1.25, "cm"), barheight = unit(0.11, "cm"))) +
    # label -30 / 0 / 30 only, so labels of neighbouring panels cannot touch
    scale_x_continuous(limits = bin_range, expand = expansion(0), breaks = c(-30, 0, 30)) +
    scale_y_continuous(n.breaks = 4) +
    labs(x = if (top) NULL else "Trait value", y = type_names[[tp]]) +
    theme_bw(base_size = txt_pt) +
    fig_theme +
    theme(panel.grid = element_blank(),
          # colour bar in the top-right panel (Copy-the-extreme, phenotype)
          legend.position = "inside", legend.position.inside = c(0.995, 0.99),
          legend.justification.inside = c(1, 1))
  # top row shares the x axis of the row below: no tick labels; bottom row: no column titles
  g <- g + if (top) theme(axis.text.x = element_blank()) else theme(strip.text = element_blank())
  g
}
p <- (dist_row("phenotype", TRUE) / dist_row("genotype", FALSE)) + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
save_widths(p, "equilibrium_dist_n8", 1.74)
cat(sprintf("rows 1-2: n = 8, distributions averaged over generations %d-%d, %d reps pooled\n",
            gen_window[1], gen_window[2], uniqueN(d$rep)))

# ---- equilibrium genetic variance relative to mutation-drift balance ----
# Neutral mutation-drift balance for an additive trait in a diploid Wright-Fisher population:
#   Vg = 2 Ne Vm,  Vm = 2 U E[a^2]  (new genetic variance per generation; U per haplosome)
# with Ne = N_pop. Per sim, Vg = mean of g_sd^2 over gen_window; then mean +- 95% CI over reps.
meta <- as.data.table(readRDS(file.path(data_dir, sprintf("metadata_%s.rds", vg_meta_tag))))
stopifnot(uniqueN(meta[, .(N_pop, U)]) == 1)
Ne <- meta$N_pop[1]; U <- meta$U[1]
Vm <- 2 * U * Ea2
Vg_md <- 2 * Ne * Vm
cat(sprintf("mutation-drift balance: Ne = %d, U = %g, Vm = %g, 2 Ne Vm = %g\n", Ne, U, Vm, Vg_md))

s <- as.data.table(readRDS(file.path(data_dir, sprintf("data_summary_%s.rds", vg_tag))))
s <- s[abs(shift_size - vg_shift) < 1e-9 & n %in% vg_n & tick >= gen_window[1] & tick <= gen_window[2]]
v <- s[, .(ratio = mean(g_sd^2) / Vg_md), by = .(reg, n, psi, rep)]
v <- v[, .(mean = mean(ratio), se = sd(ratio) / sqrt(.N), reps = .N), by = .(reg, n, psi)]
v[, `:=`(reg = factor(reg_names[reg], levels = reg_names), n_f = factor(n, levels = vg_n))]
cat(sprintf("Vg ratio: %d-%d reps per point\n", min(v$reps), max(v$reps)))

# colour bar: one reversed-Greek colour per group size tested in the dataset (light = small n);
# lines take the colour of their n, and only the plotted n are labelled
n_tested <- c(2, 4, 5, 8, 10, 20, 25, 40, 50, 100)
stopifnot(all(vg_n %in% n_tested))
n_pal    <- rev(as.character(met.brewer("Greek", length(n_tested))))
v[, n_idx := match(n, n_tested)]
dodge <- position_dodge(width = 0.03)
p_vg <- ggplot(v, aes(psi, mean, colour = n_idx, group = n_f)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = ref_lwd) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                width = 0.03, linewidth = err_lwd, position = dodge) +
  geom_line(linewidth = mean_lwd, position = dodge) +
  geom_point(size = point_size, position = dodge) +
  facet_wrap(~ reg, nrow = 1) +
  scale_colour_gradientn("n", colours = n_pal, values = scales::rescale(seq_along(n_tested)),
                         limits = c(1, length(n_tested)),
                         breaks = match(vg_n, n_tested), labels = vg_n,
                         guide = guide_colourbar(nbin = length(n_tested), display = "rectangles",
                                                 direction = "horizontal", title.position = "left",
                                                 theme = theme(legend.title = element_text(size = sym_pt, vjust = 0.75,
                                                                                           margin = margin(r = 2))),
                                                 barwidth = unit(1.0, "cm"), barheight = unit(0.11, "cm"))) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = expansion(0.06)) +
  labs(x = expression(psi), y = expression(bar(V)[g] / (2 * N[e] * V[m]))) +
  theme_bw(base_size = txt_pt) +
  fig_theme +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_text(size = sym_pt),
        # n colour bar in the empty upper-left of the Copy-the-extreme panel (white ticks = plotted n)
        legend.position = "inside", legend.position.inside = c(0.685, 0.9),
        legend.justification.inside = c(0, 1))

p_all <- (p / p_vg) + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
save_widths(p_all, "equilibrium_dist_vg_n8", 2.66)
cat(sprintf("row 3: Vg = mean g_sd^2 over generations %d-%d, mean +- 95%% CI over %d reps; 2NeVm = %g\n",
            gen_window[1], gen_window[2], max(v$reps), Vg_md))
