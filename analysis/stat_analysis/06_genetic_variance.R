# 6. genetic variance (always uses g_sd, independent of `trait`)
# Section 6 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/06_genetic_variance.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
source("analysis/stat_analysis/00_load_data.R")

# Settings for this section. Results go to their own folder inside base_plot_folder.
gv_folder <- file.path(base_plot_folder, "genetic_variance")
dir.create(gv_folder, recursive = TRUE, showWarnings = FALSE)
gv_path <- function(name) file.path(gv_folder, name)

# Neutral expectation V_g = 2 * N_e * V_m (diploid, additive), from the SLiM scripts:
# N = 1000; U = 0.03 mutations per haploid genome per generation; effect a = +/- sqrt(Exp(1)),
# so E[a^2] = 1 and V_m = 2 * U * E[a^2] per diploid per generation.
# N_e = N assumes an ideal Wright-Fisher population; selection on the group can lower it.
N_e <- 1000
V_m <- 2 * 0.03 * 1
V_g_neutral <- 2 * N_e * V_m

# Per-run genetic values: runs are psi/n/reg/rep; one shift_size kept per run (see above).
gv_df <- cons_data_df[cons_data_df[[tick_col]] < burnin_end & cons_data_df$reg %in% models,
                      c(tick_col, "psi", "n", "reg", "rep", "shift_size", "g_sd")]
names(gv_df)[1] <- "tick"
gv_key  <- interaction(gv_df$psi, gv_df$n, gv_df$reg, gv_df$rep, drop = TRUE)
gv_df   <- gv_df[gv_df$shift_size == ave(gv_df$shift_size, gv_key, FUN = function(x) min(x, na.rm = TRUE)), ]
gv_df   <- gv_df[!duplicated(gv_df[c("psi", "n", "reg", "rep", "tick")]), ]
gv_df   <- gv_df[order(gv_df$reg, gv_df$psi, gv_df$n, gv_df$rep, gv_df$tick), ]

gv_run <- do.call(rbind, lapply(split(gv_df, list(gv_df$reg, gv_df$psi, gv_df$n, gv_df$rep), drop = TRUE),
  function(d) {
    data.frame(model = d$reg[1], psi = d$psi[1], n = d$n[1], rep = d$rep[1],
               g_sd  = mean(d$g_sd,   na.rm = TRUE),   # mean genetic SD over the window
               g_var = mean(d$g_sd^2, na.rm = TRUE))   # mean genetic variance over the window
  }))
gv_run$model <- factor(gv_run$model, levels = models)
rownames(gv_run) <- NULL
write.csv(gv_run, gv_path("gv_per_run.csv"), row.names = FALSE)

# Save one PNG per model plus a multi-page PDF with the same plots.
save_by_model <- function(plots, stem, width = 7, height = 6) {
  for (mdl in names(plots)) ggsave(gv_path(paste0(stem, "_", mdl, ".pdf")), plots[[mdl]],
                                   width = width, height = height, dpi = 300)
  pdf(gv_path(paste0(stem, "_by_model.pdf")), width = width, height = height)
  for (p in plots) print(p)
  dev.off()
}

# Overview: mean g_sd (+/- 95% CI across runs) vs psi, one line per n, one panel per model.
# The dashed line is the neutral ceiling sqrt(2 N_e V_m).
gv_sum <- aggregate(g_sd ~ model + psi + n, data = gv_run,
                    FUN = function(x) c(m = mean(x), se = sd(x) / sqrt(length(x))))
gv_sum <- do.call(data.frame, gv_sum)
p_overview <- ggplot(gv_sum, aes(psi, g_sd.m, color = factor(n), group = factor(n))) +
  geom_hline(yintercept = sqrt(V_g_neutral), linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = g_sd.m - 1.96 * g_sd.se, ymax = g_sd.m + 1.96 * g_sd.se), width = 0.01) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
  facet_wrap(~ model) + scale_color_viridis_d(end = 0.9) +
  labs(x = expression(psi), y = "Genetic SD", color = "n",
       caption = "Mean over runs, 95% CI; dashed = neutral ceiling sqrt(2 Ne Vm)") +
  theme_bw()
ggsave(gv_path("g_sd_vs_psi.pdf"), p_overview, width = 11, height = 4.5, dpi = 300)
p_overview


# ---- 6.0.1 genetic SD vs psi, selected n ----
# Same plot as g_sd_vs_psi.pdf, restricted to the n values listed here.
n_sel <- c(2, 8, 25, 100)
n_sel <- sort(intersect(n_sel, unique(gv_sum$n)))   # drop values not in the data
gv_sel <- gv_sum[gv_sum$n %in% n_sel, ]

p_sel_base <- ggplot(gv_sel, aes(psi, g_sd.m, group = factor(n))) +
  geom_hline(yintercept = sqrt(V_g_neutral), linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = g_sd.m - 1.96 * g_sd.se, ymax = g_sd.m + 1.96 * g_sd.se), width = 0.01) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
  facet_wrap(~ model) +
  labs(x = expression(psi), y = "Genetic SD", color = "n",
       caption = "Mean over runs, 95% CI; dashed = neutral ceiling sqrt(2 Ne Vm)") +
  theme_bw()

# Option A: continuous log-scale color bar, labelled at the selected n.
p_sel_bar <- p_sel_base + aes(color = n) + labs(color = "n") +
  scale_color_viridis_c(end = 0.9, trans = "log10", breaks = n_sel, labels = n_sel) +
  guides(color = guide_colorbar(barwidth = unit(12, "lines"), barheight = unit(0.6, "lines"),
                                title.vjust = 0.8)) +
  theme(legend.position = "bottom")
ggsave(gv_path("g_sd_vs_psi_sel_colorbar.pdf"), p_sel_bar, width = 11, height = 5)
p_sel_bar

# Option B: discrete legend as a single row below the plot.
p_sel_row <- p_sel_base + aes(color = factor(n)) + labs(color = "n") +
  scale_color_viridis_d(end = 0.9) +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")
ggsave(gv_path("g_sd_vs_psi_sel_legendrow.pdf"), p_sel_row, width = 11, height = 5)
p_sel_row

# Options C/D: the two legends above, placed inside the top of the ext panel
# (ext is the rightmost panel; position is in plot-area fractions, top-right corner).
legend_in_ext <- theme(legend.position = c(0.995, 0.02), legend.justification = c(1, 0),
                       legend.direction = "horizontal",
                       legend.background = element_rect(fill = alpha("white", 0.8), color = NA))
p_sel_bar_in <- p_sel_bar + legend_in_ext +
  guides(color = guide_colorbar(barwidth = unit(8, "lines"), barheight = unit(0.5, "lines"),
                                title.vjust = 0.8))
ggsave(gv_path("g_sd_vs_psi_sel_colorbar_inext.pdf"), p_sel_bar_in, width = 11, height = 4.5)
p_sel_bar_in

p_sel_row_in <- p_sel_row + legend_in_ext
ggsave(gv_path("g_sd_vs_psi_sel_legendrow_inext.pdf"), p_sel_row_in, width = 11, height = 4.5)
p_sel_row_in

# Options E/F: stepped color bar with one tile per n tested in the dataset (colors
# spread over all n, so a plotted line matches its tile), labelled only at n_sel.
# Only the lines feed the legend, drawn as rect keys so adjacent tiles form a bar.
n_all <- sort(unique(gv_sum$n))
p_sel_tiles <- ggplot(gv_sel, aes(psi, g_sd.m, color = factor(n, levels = n_all), group = factor(n))) +
  geom_hline(yintercept = sqrt(V_g_neutral), linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = g_sd.m - 1.96 * g_sd.se, ymax = g_sd.m + 1.96 * g_sd.se),
                width = 0.01, show.legend = FALSE) +
  geom_line(linewidth = 0.6, key_glyph = "rect") +
  geom_point(size = 1.5, show.legend = FALSE) +
  facet_wrap(~ model) +
  scale_color_viridis_d(end = 0.9, drop = FALSE,
                        labels = function(x) ifelse(x %in% as.character(n_sel), x, "")) +
  guides(color = guide_legend(nrow = 1, label.position = "bottom", label.hjust = 0.5,
                              title.vjust = 0.8)) +
  labs(x = expression(psi), y = "Genetic SD", color = "n",
       caption = "Mean over runs, 95% CI; dashed = neutral ceiling sqrt(2 Ne Vm)") +
  theme_bw() +
  theme(legend.spacing.x = unit(0, "pt"),
        legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.6, "lines"))

p_sel_tiles_bottom <- p_sel_tiles + theme(legend.position = "bottom")
ggsave(gv_path("g_sd_vs_psi_sel_tiles.pdf"), p_sel_tiles_bottom, width = 11, height = 5)
p_sel_tiles_bottom

p_sel_tiles_in <- p_sel_tiles + legend_in_ext +
  theme(legend.key.width = unit(0.7, "lines"), legend.key.height = unit(0.5, "lines"))
ggsave(gv_path("g_sd_vs_psi_sel_tiles_inext.pdf"), p_sel_tiles_in, width = 11, height = 4.5)
p_sel_tiles_in


# ---- 6.0.2 genetic SD relative to the neutral ceiling, selected n ----
# Same selected-n plot as 6.0.1, with g_sd divided by the neutral ceiling sqrt(2 N_e V_m),
# so 1 = neutral expectation. The ceiling is a constant, so the mean and its CI are
# scaled by the same factor; the CI can exceed 1. Linear y axis so the ceiling reads directly.
# Summarise per-run genetic variance across replicate runs; t-based 95% CI.
gv_sum_v <- aggregate(g_var ~ model + psi + n, data = gv_run,
                      FUN = function(x) c(m = mean(x), se = sd(x) / sqrt(length(x)), k = length(x)))
gv_sum_v <- do.call(data.frame, gv_sum_v)
tq <- qt(0.975, df = gv_sum_v$g_var.k - 1)

# Ratio to the neutral expectation on the variance scale (1 = neutral),
# then sqrt of estimate and endpoints to put it back on the SD scale.
gv_sum_v <- transform(gv_sum_v,
  rel_var    = g_var.m / V_g_neutral,
  rel_var_lo = (g_var.m - tq * g_var.se) / V_g_neutral,
  rel_var_hi = (g_var.m + tq * g_var.se) / V_g_neutral)
gv_sum_v <- transform(gv_sum_v,
  rel    = sqrt(rel_var),
  rel_lo = sqrt(pmax(rel_var_lo, 0)),
  rel_hi = sqrt(rel_var_hi))

gv_sel_rel <- gv_sum_v[gv_sum_v$n %in% n_sel, ]

p_rel_tiles <- ggplot(gv_sel_rel, aes(psi, rel_var, color = factor(n, levels = n_all), group = factor(n))) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = rel_var_lo, ymax = rel_var_hi), width = 0.01, show.legend = FALSE) +
  geom_line(linewidth = 0.6, key_glyph = "rect") +
  geom_point(size = 1.5, show.legend = FALSE) +
  facet_wrap(~ model) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  coord_cartesian(ylim = c(0.6, NA)) +   # zooms only; drop this line to show the full range
  scale_color_viridis_d(end = 0.9, drop = FALSE,
                        labels = function(x) ifelse(x %in% as.character(n_sel), x, "")) +
  guides(color = guide_legend(nrow = 1, label.position = "bottom", label.hjust = 0.5,
                              title.vjust = 0.8)) +
  labs(x = expression(psi), y = expression(bar(V)[g] / (2 * N[e] * V[m])), color = "n",
       caption = "Mean genetic variance over runs relative to neutral expectation 2NeVm; 95% t-interval; dashed = 1") +
  theme_bw() +
  theme(legend.spacing.x = unit(0, "pt"),
        legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.6, "lines"))
p_rel_tiles_bottom <- p_rel_tiles + theme(legend.position = "bottom")
ggsave(gv_path("g_sd_rel_neutral_vs_psi_sel_tiles.pdf"), p_rel_tiles_bottom, width = 11, height = 5)
p_rel_tiles_bottom

p_rel_tiles_in <- p_rel_tiles + legend_in_ext +
  theme(legend.key.width = unit(0.7, "lines"), legend.key.height = unit(0.5, "lines"))
ggsave(gv_path("g_sd_rel_neutral_vs_psi_sel_tiles_inext.pdf"), p_rel_tiles_in, width = 11, height = 4.5)
p_rel_tiles_in


# ---- 6.1 genetic SD relative to psi = 0 ----
# Ratio mean g_sd(psi) / mean g_sd(psi = 0) within each model and n; bootstrap
# resamples runs within each psi group. Classified E/D/? like the model ratios.
psi0_list <- list()
for (mdl in models) for (nn in sort(unique(gv_run$n))) {
  sub <- gv_run[gv_run$model == mdl & gv_run$n == nn, ]
  ref <- sub$g_sd[sub$psi == 0]
  if (length(ref) < 2) next
  for (ps in setdiff(sort(unique(sub$psi)), 0)) {
    x <- sub$g_sd[sub$psi == ps]
    if (length(x) < 2) next
    ci <- boot_diff_ratio(ref, x, n_boot)
    psi0_list[[length(psi0_list) + 1]] <- data.frame(
      model = mdl, n = nn, psi = ps, ratio = mean(x) / mean(ref),
      ratio_lo90 = ci[["ratio_lo90"]], ratio_hi90 = ci[["ratio_hi90"]],
      ratio_lo   = ci[["ratio_lo"]],   ratio_hi   = ci[["ratio_hi"]])
  }
}
psi0_df <- do.call(rbind, psi0_list)
psi0_df$class <- classify_equiv(psi0_df$ratio_lo90, psi0_df$ratio_hi90, equiv_margin)
write.csv(psi0_df, gv_path("g_sd_ratio_vs_psi0.csv"), row.names = FALSE)

# One heatmap per model, each with its own color bar (psi = 0 is the reference row, not shown).
psi0_plots <- setNames(lapply(models, function(mdl)
  ratio_heatmap(psi0_df[psi0_df$model == mdl, ],
                paste0("Genetic SD relative to psi = 0 (", mdl, ")"),
                fill_name = "g_sd ratio")), models)
save_by_model(psi0_plots, "g_sd_ratio_vs_psi0_heatmap")


# ---- 6.2 trend test: does g_sd increase with psi? ----
# For each model and n (runs as units):
#   Spearman rho between psi and per-run g_sd;
#   Jonckheere-Terpstra test for an increasing trend across ordered psi groups
#   (one-sided, normal approximation; exact ties are negligible for continuous g_sd);
#   effect size = slope of log(g_sd) on psi (lm, 95% CI), also as % change per 0.1 psi.
jt_test <- function(x, g) {
  lv <- sort(unique(g)); ni <- as.vector(table(factor(g, levels = lv))); N <- length(x)
  J <- 0
  for (i in seq_along(lv)[-length(lv)]) for (j in (i + 1):length(lv)) {
    xi <- x[g == lv[i]]; xj <- x[g == lv[j]]
    J <- J + sum(outer(xi, xj, "<")) + 0.5 * sum(outer(xi, xj, "=="))
  }
  EJ <- (N^2 - sum(ni^2)) / 4
  VJ <- (N^2 * (2 * N + 3) - sum(ni^2 * (2 * ni + 3))) / 72
  z <- (J - EJ) / sqrt(VJ)
  c(JT = J, z = z, p = pnorm(z, lower.tail = FALSE))
}

trend_list <- list()
for (mdl in models) for (nn in sort(unique(gv_run$n))) {
  sub <- gv_run[gv_run$model == mdl & gv_run$n == nn & gv_run$g_sd > 0, ]
  if (length(unique(sub$psi)) < 3) next
  sp  <- suppressWarnings(cor.test(sub$psi, sub$g_sd, method = "spearman", exact = FALSE))
  jt  <- jt_test(sub$g_sd, sub$psi)
  fit <- lm(log(g_sd) ~ psi, data = sub)
  ci  <- confint(fit)["psi", ]
  trend_list[[length(trend_list) + 1]] <- data.frame(
    model = mdl, n = nn, n_runs = nrow(sub),
    slope_log_gsd = coef(fit)[["psi"]], slope_lo = ci[[1]], slope_hi = ci[[2]],
    pct_per_0.1psi = 100 * (exp(0.1 * coef(fit)[["psi"]]) - 1),
    spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
    jt_z = jt[["z"]], jt_p_increasing = jt[["p"]])
}
trend_df <- do.call(rbind, trend_list)
# Holm across n within each model
trend_df$spearman_p_holm <- ave(trend_df$spearman_p,      trend_df$model, FUN = function(p) p.adjust(p, "holm"))
trend_df$jt_p_holm       <- ave(trend_df$jt_p_increasing, trend_df$model, FUN = function(p) p.adjust(p, "holm"))
write.csv(trend_df, gv_path("g_sd_trend_tests.csv"), row.names = FALSE)


# ---- 6.3 mechanism check: relaxed selection on own genotype ----
# If a fraction psi of the phenotype no longer comes from own genotype, selection on
# own genotype is weakened by (1 - psi)^2, i.e. effective Vs -> Vs / (1 - psi)^2.
# Plot log(V_g) against x = -2 log(1 - psi) for psi < 1:
#   slope ~ 1   -> V_g proportional to Vs (house-of-cards regime)
#   slope ~ 0.5 -> V_g proportional to sqrt(Vs) (Gaussian / Kimura regime)
mech <- gv_run[gv_run$psi < 1 & gv_run$g_var > 0, ]
mech$x <- -2 * log(1 - mech$psi)
mech$y <- log(mech$g_var)

mech_list <- list()
for (mdl in models) for (nn in sort(unique(mech$n))) {
  sub <- mech[mech$model == mdl & mech$n == nn, ]
  if (length(unique(sub$x)) < 2) next
  fit <- lm(y ~ x, data = sub); ci <- confint(fit)["x", ]
  mech_list[[length(mech_list) + 1]] <- data.frame(
    model = mdl, n = nn, slope = coef(fit)[["x"]], slope_lo = ci[[1]], slope_hi = ci[[2]],
    r_squared = summary(fit)$r.squared)
}
mech_df <- do.call(rbind, mech_list)
write.csv(mech_df, gv_path("mechanism_slopes.csv"), row.names = FALSE)

# Per-psi means over runs, with the fitted line per n and a slope-1 reference line
# anchored at the psi = 0 mean (dashed).
mech_mean <- aggregate(y ~ model + n + x, data = mech, FUN = mean)
ref1 <- merge(mech_mean[mech_mean$x == 0, c("model", "n", "y")], data.frame(slope = 1))
p_mech <- ggplot(mech, aes(x, y, color = factor(n))) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_point(data = mech_mean, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
  geom_abline(data = ref1, aes(intercept = y, slope = slope, color = factor(n)),
              linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  facet_wrap(~ model) + scale_color_viridis_d(end = 0.9) +
  labs(x = expression(-2 * log(1 - psi)), y = expression(log(g[sd]^2)), color = "n",
       caption = "Small points = runs, large = mean over runs; solid = fit, dashed = slope 1 from psi = 0") +
  theme_bw()
ggsave(gv_path("mechanism_logVg_vs_log1mpsi.pdf"), p_mech, width = 11, height = 4.5, dpi = 300)
p_mech


# ---- 6.4 neutral ceiling at psi = 1 ----
# Observed mean V_g at psi = 1 vs the neutral expectation 2 N_e V_m, with a bootstrap
# CI over runs. Note: at psi = 1 own genotype has no direct effect on own phenotype,
# but in fit/ext a genotype still matters when it is the one copied, so full
# neutrality is not guaranteed.
neut_list <- list()
at1 <- gv_run[gv_run$psi == 1, ]
for (mdl in models) for (nn in sort(unique(at1$n))) {
  v <- at1$g_var[at1$model == mdl & at1$n == nn]
  if (length(v) < 2) next
  bm <- colMeans(matrix(v[sample.int(length(v), length(v) * n_boot, replace = TRUE)], nrow = length(v)))
  neut_list[[length(neut_list) + 1]] <- data.frame(
    model = mdl, n = nn, n_runs = length(v), mean_Vg = mean(v),
    expected_Vg = V_g_neutral, ratio = mean(v) / V_g_neutral,
    ratio_lo = quantile(bm, 0.025, names = FALSE) / V_g_neutral,
    ratio_hi = quantile(bm, 0.975, names = FALSE) / V_g_neutral)
}
neut_df <- if (length(neut_list)) do.call(rbind, neut_list) else NULL
if (!is.null(neut_df)) write.csv(neut_df, gv_path("neutral_ceiling_psi1.csv"), row.names = FALSE)


# ---- 6.5 effect of group size on genetic SD ----
# Within each model and psi (runs as units):
#   Kruskal-Wallis test of g_sd across n;
#   effect size = slope of log(g_sd) on log(n) (lm, 95% CI): the % change in g_sd per doubling of n;
#   ratio g_sd(n) / g_sd(smallest n) with bootstrap CI and E/D/? label, shown as a heatmap.
n_ref <- min(gv_run$n)
gs_list <- list(); nratio_list <- list()
for (mdl in models) for (ps in sort(unique(gv_run$psi))) {
  sub <- gv_run[gv_run$model == mdl & gv_run$psi == ps & gv_run$g_sd > 0, ]
  if (length(unique(sub$n)) < 2) next
  kw  <- kruskal.test(g_sd ~ factor(n), data = sub)
  fit <- lm(log(g_sd) ~ log2(n), data = sub); ci <- confint(fit)["log2(n)", ]
  gs_list[[length(gs_list) + 1]] <- data.frame(
    model = mdl, psi = ps, n_runs = nrow(sub),
    slope_log_gsd_per_log2n = coef(fit)[["log2(n)"]], slope_lo = ci[[1]], slope_hi = ci[[2]],
    pct_per_doubling_n = 100 * (exp(coef(fit)[["log2(n)"]]) - 1),
    kw_chisq = unname(kw$statistic), kw_df = unname(kw$parameter), kw_p = kw$p.value)

  ref <- sub$g_sd[sub$n == n_ref]
  if (length(ref) < 2) next
  for (nn in setdiff(sort(unique(sub$n)), n_ref)) {
    x <- sub$g_sd[sub$n == nn]
    if (length(x) < 2) next
    ci <- boot_diff_ratio(ref, x, n_boot)
    nratio_list[[length(nratio_list) + 1]] <- data.frame(
      model = mdl, psi = ps, n = nn, ratio = mean(x) / mean(ref),
      ratio_lo90 = ci[["ratio_lo90"]], ratio_hi90 = ci[["ratio_hi90"]],
      ratio_lo   = ci[["ratio_lo"]],   ratio_hi   = ci[["ratio_hi"]])
  }
}
gs_df <- do.call(rbind, gs_list)
gs_df$kw_p_holm <- ave(gs_df$kw_p, gs_df$model, FUN = function(p) p.adjust(p, "holm"))
write.csv(gs_df, gv_path("group_size_tests.csv"), row.names = FALSE)

nratio_df <- do.call(rbind, nratio_list)
nratio_df$class <- classify_equiv(nratio_df$ratio_lo90, nratio_df$ratio_hi90, equiv_margin)
write.csv(nratio_df, gv_path("g_sd_ratio_vs_smallest_n.csv"), row.names = FALSE)
nratio_plots <- setNames(lapply(models, function(mdl)
  ratio_heatmap(nratio_df[nratio_df$model == mdl, ],
                paste0("Genetic SD relative to n = ", n_ref, " (", mdl, ")"),
                fill_name = "g_sd ratio")), models)
save_by_model(nratio_plots, "g_sd_ratio_vs_smallest_n_heatmap")


# ---- 6.6 genetic variance report ----
gv_rep <- c(
  "==================================================================",
  " GENETIC VARIANCE (burn-in, generation < 10000)",
  "==================================================================",
  sprintf(" Per-run value: mean g_sd (and mean g_sd^2) over all time points."),
  " Unit of analysis: one run. Bootstrap CIs resample runs within groups",
  sprintf(" (%d resamples). E/D/? uses the 90%% CI vs margin [%.2f, %.2f].",
          n_boot, equiv_margin[1], equiv_margin[2]),
  ""
)

# 1. ratio to psi = 0
gv_rep <- c(gv_rep, rule, " 1. Genetic SD relative to psi = 0: ratio [90% CI] [95% CI]", rule)
for (mdl in models) {
  d <- psi0_df[psi0_df$model == mdl, ]
  if (!nrow(d)) next
  cnt <- table(factor(d$class, levels = names(class_code)))
  gv_rep <- c(gv_rep, "", sprintf(" %s: equivalent %d, different %d, inconclusive %d (of %d cells)",
                                  mdl, cnt[["equivalent"]], cnt[["different"]], cnt[["inconclusive"]], nrow(d)))
  grid <- tapply(sprintf("%.2f%s", d$ratio, class_code[d$class]), list(d$psi, d$n), identity)
  grid[is.na(grid)] <- "-"
  gv_rep <- c(gv_rep, "   Grid (rows psi, columns n):", paste0("     ", capture.output(print(noquote(grid)))))
}

# 2. trend
gv_rep <- c(gv_rep, "", rule, " 2. Trend of g_sd with psi (per n)",
            "    effect size: slope of log(g_sd) on psi [95% CI], % change per 0.1 psi",
            "    tests: Spearman rho, Jonckheere-Terpstra (one-sided, increasing); Holm across n", rule)
for (i in seq_len(nrow(trend_df))) {
  t <- trend_df[i, ]
  gv_rep <- c(gv_rep, sprintf("   %-4s n = %-4s slope = %s (%+.1f%% per 0.1 psi); rho = %s, p_holm %s; JT z = %s, p_holm %s",
                              t$model, t$n, ci3(t$slope_log_gsd, t$slope_lo, t$slope_hi), t$pct_per_0.1psi,
                              f3(t$spearman_rho), fmt_p(t$spearman_p_holm), f3(t$jt_z), fmt_p(t$jt_p_holm)))
}

# 3. mechanism
gv_rep <- c(gv_rep, "", rule, " 3. Mechanism: slope of log(g_sd^2) on -2 log(1 - psi), psi < 1",
            "    ~1 = V_g proportional to Vs (house of cards); ~0.5 = proportional to sqrt(Vs) (Gaussian)", rule)
for (i in seq_len(nrow(mech_df))) {
  t <- mech_df[i, ]
  gv_rep <- c(gv_rep, sprintf("   %-4s n = %-4s slope = %s, R^2 = %s",
                              t$model, t$n, ci3(t$slope, t$slope_lo, t$slope_hi), f3(t$r_squared)))
}

# 4. neutral ceiling
gv_rep <- c(gv_rep, "", rule,
            sprintf(" 4. Neutral ceiling at psi = 1: V_g expected = 2 Ne Vm = 2 x %g x %g = %g", N_e, V_m, V_g_neutral),
            "    observed mean V_g / expected [95% CI]", rule)
if (is.null(neut_df)) {
  gv_rep <- c(gv_rep, "   no runs at psi = 1")
} else for (i in seq_len(nrow(neut_df))) {
  t <- neut_df[i, ]
  gv_rep <- c(gv_rep, sprintf("   %-4s n = %-4s V_g = %-9s ratio = %s",
                              t$model, t$n, fmt(t$mean_Vg), ci3(t$ratio, t$ratio_lo, t$ratio_hi)))
}

# 5. group size
gv_rep <- c(gv_rep, "", rule, " 5. Effect of group size n on g_sd (within each psi)",
            "    effect size: % change in g_sd per doubling of n [95% CI of log slope]",
            "    test: Kruskal-Wallis across n, Holm across psi within model", rule)
for (i in seq_len(nrow(gs_df))) {
  t <- gs_df[i, ]
  gv_rep <- c(gv_rep, sprintf("   %-4s psi = %-5s %+6.1f%% per doubling (log slope %s); KW chi2(%d) = %s, p_holm %s",
                              t$model, t$psi, t$pct_per_doubling_n,
                              ci3(t$slope_log_gsd_per_log2n, t$slope_lo, t$slope_hi),
                              t$kw_df, fmt(t$kw_chisq), fmt_p(t$kw_p_holm)))
}
for (mdl in models) {
  d <- nratio_df[nratio_df$model == mdl, ]
  if (!nrow(d)) next
  cnt <- table(factor(d$class, levels = names(class_code)))
  gv_rep <- c(gv_rep, "", sprintf(" %s, g_sd(n) / g_sd(n = %s): equivalent %d, different %d, inconclusive %d (of %d cells)",
                                  mdl, n_ref, cnt[["equivalent"]], cnt[["different"]], cnt[["inconclusive"]], nrow(d)))
  grid <- tapply(sprintf("%.2f%s", d$ratio, class_code[d$class]), list(d$psi, d$n), identity)
  grid[is.na(grid)] <- "-"
  gv_rep <- c(gv_rep, "   Grid (rows psi, columns n):", paste0("     ", capture.output(print(noquote(grid)))))
}

cat(gv_rep, sep = "\n")
writeLines(gv_rep, gv_path("genetic_variance_report.txt"))
