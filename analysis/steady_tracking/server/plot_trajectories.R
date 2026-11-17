#!/usr/bin/env Rscript
# Trajectories of genetic-variation, allele-effect and selection quantities for every parameter
# combination of the steady_tracking run (shift_speed x reg x psi, 60 reps each).
#
# Input (from analysis/load_data/server/load_data.R -d steady_tracking -t data,metadata --tag all):
#   figure/r_data/steady_tracking/data_summary_all.rds
# Output: analysis/steady_tracking/server/plots_fit500/trajectories/
#   <var>_summary.pdf     panel per reg x speed, line per psi: mean over reps, ribbon = 10-90%
#                         quantiles across reps
#   trajectories_all.pdf  all of the above, one page per quantity
#   <var>_reps_<reg>.pdf  every rep (only for the quantities in rep_vars)
#
# Not plotted, and why: H, PI and Ne_heterozygosity (every QTL is under selection and there are no
# neutral sites to compare against, so H and an Ne from it say nothing about drift),
# g_opt (= z_opt, since g_mean == z_mean in this model),
# z_opt_closest, parent_* (parent_z_sd / parent_g_sd are empty), abs_effect_* (mostly the
# growing sum over fixed mutations), fitness_sd, sfs (a string), *_sd of per-ind means.
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_trajectories.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")   # read_summary(), read_meta(): both runs

# ---- settings ----
data_dir   <- "figure/r_data/steady_tracking"
out_dir    <- "analysis/steady_tracking/server/plots_fit500/trajectories"
shift_tick <- 10000
tick_step  <- 5            # per-rep figures: every 5th tick
rep_vars   <- c("seg_site_pop")        # quantities that also get per-rep figures

# quantities to plot: column (in data, or derived below), axis label, rolling-mean window in
# ticks applied to each sim before averaging (1 = none; for quantities that are noisy per tick),
# and an optional reference line
vars <- rbindlist(list(
  # genetic variation
  list("seg_site_pop",       "segregating sites (population)",                     1,  NA),
  list("num_seg_mean",       "unique mutations per individual (mean)",             1,  NA),
  list("fixed_since_shift",  "mutations fixed since tick 9999",                    1,  NA),
  list("g_sd",               "SD of genotypic value",                              1,  NA),
  list("z_sd",               "SD of phenotype",                                    1,  NA),
  list("group_vars_mean",    "within-group variance (mean over individuals)",      1,  NA),
  # allele effects
  list("mean_effect_mean",   "mean effect of an individual's mutations",           1,  0),
  list("mean_abs_effect_mean", "mean |effect| of an individual's mutations",       1,  NA),
  # selection and response
  list("s",                  "selection differential on phenotype (20-tick mean)", 20, 0),
  list("b_gz",               "regression of genotype on phenotype b_gz (20-tick mean)", 20, 1),
  list("dg_obs",             "change in mean genotype per tick (20-tick mean)",    20, NA),
  list("response_ratio",     "dg_obs / shift_speed (1 = keeping pace; 20-tick mean)", 20, 1)
))
setnames(vars, c("var", "label", "smooth", "ref"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_summary()
W <- speed_width(uniqueN(d$shift_speed))   # page width: one facet column per speed
base_cols <- intersect(vars$var, names(d))
d <- d[, c("shift_speed", "reg", "psi", "rep", "ID", "tick", base_cols, "num_fixed_mut"), with = FALSE]
setkey(d, ID, tick)
d[, fixed_since_shift := num_fixed_mut - num_fixed_mut[tick == shift_tick - 1], by = ID]
d[, response_ratio := dg_obs / shift_speed]
d[, num_fixed_mut := NULL]

# rolling means (centred) per sim for the noisy quantities
for (i in which(vars$smooth > 1)) {
  v <- vars$var[i]; k <- vars$smooth[i]
  d[, (v) := frollmean(get(v), k, align = "center"), by = ID]
}

speeds <- sort(unique(d$shift_speed))
d[, `:=`(psi_f = factor(psi),
         speed_f = factor(shift_speed, levels = speeds, labels = paste("speed", speeds)),
         reg = reg_factor(reg))]

base_theme <- list(
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_text(angle = 45, hjust = 1))
)
shift_line <- geom_vline(xintercept = shift_tick, linetype = "dotted", colour = "grey40")

summary_plot <- function(v, label, ref) {
  s <- d[!is.na(get(v)), .(mean = mean(get(v)), lo = quantile(get(v), 0.1), hi = quantile(get(v), 0.9)),
         by = .(reg, speed_f, psi_f, tick)]
  ggplot(s, aes(tick, mean, colour = psi_num(psi_f), fill = psi_num(psi_f), group = psi_f)) +
    { if (!is.na(ref)) geom_hline(yintercept = ref, linetype = "dashed", colour = "grey40") } +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.5) +
    shift_line +
    facet_grid(reg ~ speed_f) +
    psi_scale() +
    labs(x = "tick", y = label, title = sprintf("%s  [%s]", label, v),
         subtitle = "line: mean of 60 reps, ribbon: 10-90% quantiles; dotted: optimum starts moving") +
    base_theme
}

plots <- list()
for (i in seq_len(nrow(vars))) {
  v <- vars$var[i]
  cat("plotting", v, "\n")
  p <- summary_plot(v, vars$label[i], vars$ref[i])
  ggsave(file.path(out_dir, sprintf("%s_summary.pdf", v)), p, width = W, height = 7)
  plots[[v]] <- p
}
pdf(file.path(out_dir, "trajectories_all.pdf"), width = W, height = 7)
for (p in plots) print(p)
invisible(dev.off())

# ---- per-rep figures for rep_vars ----
r <- d[(tick - min(tick)) %% tick_step == 0 | tick == max(tick)]
for (v in rep_vars) {
  lab <- vars[var == v, label]
  for (g in levels(d$reg)) {
    p <- ggplot(r[reg == g & !is.na(get(v))], aes(tick, .data[[v]], group = ID, colour = psi_num(psi_f))) +
      geom_line(linewidth = 0.2, alpha = 0.3) +
      shift_line +
      facet_grid(psi_f ~ speed_f, labeller = labeller(psi_f = function(x) paste("psi", x))) +
      psi_scale(aesthetics = "colour", guide = FALSE) +
      labs(x = "tick", y = lab, title = sprintf("%s of every rep, %s (rows: psi)", lab, g)) +
      base_theme
    ggsave(file.path(out_dir, sprintf("%s_reps_%s.pdf", v, reg_code(g))), p, width = W, height = 9)
  }
}
cat("plots written to", out_dir, "\n")
