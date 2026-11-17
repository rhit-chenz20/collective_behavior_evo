# 7. phenotypic variance (always uses z_sd and g_sd, independent of `trait`)
# Section 7 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/07_phenotypic_variance.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
source("analysis/stat_analysis/00_load_data.R")

# Settings for this section. Results go to their own folder inside base_plot_folder.
pv_folder <- file.path(base_plot_folder, "phenotypic_variance")
dir.create(pv_folder, recursive = TRUE, showWarnings = FALSE)
pv_path <- function(name) file.path(pv_folder, name)

n_mc_inds      <- 20000   # simulated individuals per prediction for the fit/ext order statistics

# Per-run values over all time points, as in section 6; one shift_size per run.
pv_df <- cons_data_df[cons_data_df[[tick_col]] < burnin_end & cons_data_df$reg %in% models,
                      c(tick_col, "psi", "n", "reg", "rep", "shift_size", "z_sd", "g_sd", "g_mean")]
names(pv_df)[1] <- "tick"
pv_key <- interaction(pv_df$psi, pv_df$n, pv_df$reg, pv_df$rep, drop = TRUE)
pv_df  <- pv_df[pv_df$shift_size == ave(pv_df$shift_size, pv_key, FUN = function(x) min(x, na.rm = TRUE)), ]
pv_df  <- pv_df[!duplicated(pv_df[c("psi", "n", "reg", "rep", "tick")]), ]
pv_df  <- pv_df[order(pv_df$reg, pv_df$psi, pv_df$n, pv_df$rep, pv_df$tick), ]

pv_run <- do.call(rbind, lapply(split(pv_df, list(pv_df$reg, pv_df$psi, pv_df$n, pv_df$rep), drop = TRUE),
  function(d) {
    data.frame(model = d$reg[1], psi = d$psi[1], n = d$n[1], rep = d$rep[1],
               z_sd  = mean(d$z_sd,   na.rm = TRUE), z_var = mean(d$z_sd^2, na.rm = TRUE),
               g_sd  = mean(d$g_sd,   na.rm = TRUE), g_var = mean(d$g_sd^2, na.rm = TRUE),
               g_mean = mean(d$g_mean, na.rm = TRUE))
  }))
pv_run$model <- factor(pv_run$model, levels = models)
pv_run$zg_ratio <- pv_run$z_sd / pv_run$g_sd   # per-run shrinkage (<1) / amplification (>1)
rownames(pv_run) <- NULL

# Save one PNG per model plus a multi-page PDF with the same plots.
save_by_model_pv <- function(plots, stem, width = 7, height = 6) {
  for (mdl in names(plots)) ggsave(pv_path(paste0(stem, "_", mdl, ".pdf")), plots[[mdl]],
                                   width = width, height = height, dpi = 300)
  pdf(pv_path(paste0(stem, "_by_model.pdf")), width = width, height = height)
  for (p in plots) print(p)
  dev.off()
}


# ---- 7.0 phenotypic variance factor Var(P) / V_g under each social rule ----
# Factor c such that Var(phenotype) = c * V_g, for genotypes a ~ Normal within the
# population and random groups of size n (exact SLiM rules):
#   ave: P_i = [(n-1)(1-psi) a_i + psi * sum(group a)] / (n-1+psi)   -> closed form
#   fit: P_i = (1-psi) a_i + psi * a_best, best = member closest to the optimum
#   ext: P_i = (1-psi) a_i + psi * a_ext,  ext  = member farthest from the group mean
# fit and ext are order statistics, computed numerically with common random numbers
# (the same standard-normal draws for every call, so predictions vary smoothly).
# d = (g_mean - optimum) / g_sd matters only for fit.
old_seed <- .Random.seed
set.seed(1)
Z_cache <- lapply(setNames(nm = sort(unique(pv_run$n))), function(nn)
  matrix(rnorm(ceiling(n_mc_inds / nn) * nn), ncol = nn))   # rows = groups, cols = members
.Random.seed <- old_seed   # keep the bootstrap RNG stream unchanged

pheno_var_factor <- function(model, psi, n, d = 0) {
  if (model == "ave") return((((n - 1) * (1 - psi) + psi)^2 + (n - 1) * psi^2) / (n - 1 + psi)^2)
  A <- Z_cache[[as.character(n)]] + d
  pick <- if (model == "fit") max.col(-abs(A), ties.method = "first") else
                              max.col(abs(A - rowMeans(A)), ties.method = "first")
  copied <- A[cbind(seq_len(nrow(A)), pick)]
  var(as.vector((1 - psi) * A + psi * copied)) / var(as.vector(A))
}

# Theory grid at d = 0 (population mean on the optimum), used as reference curves.
theory_df <- expand.grid(model = models, psi = sort(unique(pv_run$psi)), n = sort(unique(pv_run$n)),
                         stringsAsFactors = FALSE)
theory_df$factor <- mapply(pheno_var_factor, theory_df$model, theory_df$psi, theory_df$n)
theory_df$sd_ratio_theory <- sqrt(theory_df$factor)
theory_df$model <- factor(theory_df$model, levels = models)
write.csv(theory_df, pv_path("theory_variance_factor.csv"), row.names = FALSE)


# ---- 7.1 phenotypic SD relative to psi = 0 ----
# Same as 6.1 but for z_sd: mean z_sd(psi) / mean z_sd(psi = 0) within each model and n.
pz0_list <- list()
for (mdl in models) for (nn in sort(unique(pv_run$n))) {
  sub <- pv_run[pv_run$model == mdl & pv_run$n == nn, ]
  ref <- sub$z_sd[sub$psi == 0]
  if (length(ref) < 2) next
  for (ps in setdiff(sort(unique(sub$psi)), 0)) {
    x <- sub$z_sd[sub$psi == ps]
    if (length(x) < 2) next
    ci <- boot_diff_ratio(ref, x, n_boot)
    pz0_list[[length(pz0_list) + 1]] <- data.frame(
      model = mdl, n = nn, psi = ps, ratio = mean(x) / mean(ref),
      ratio_lo90 = ci[["ratio_lo90"]], ratio_hi90 = ci[["ratio_hi90"]],
      ratio_lo   = ci[["ratio_lo"]],   ratio_hi   = ci[["ratio_hi"]])
  }
}
pz0_df <- do.call(rbind, pz0_list)
pz0_df$class <- classify_equiv(pz0_df$ratio_lo90, pz0_df$ratio_hi90, equiv_margin)
write.csv(pz0_df, pv_path("z_sd_ratio_vs_psi0.csv"), row.names = FALSE)

pz0_plots <- setNames(lapply(models, function(mdl)
  ratio_heatmap(pz0_df[pz0_df$model == mdl, ],
                paste0("Phenotypic SD relative to psi = 0 (", mdl, ")"),
                fill_name = "z_sd ratio")), models)
save_by_model_pv(pz0_plots, "z_sd_ratio_vs_psi0_heatmap")


# ---- 7.2 phenotypic SD against n ----
# Left: observed z_sd vs n (log-log), one line per psi, one panel per model.
# Right: z_sd / g_sd vs n, with the order-statistic theory sqrt(c(psi, n)) as dashed lines;
#   dividing by g_sd removes the change in V_g with n, isolating the social rule.
# Expected scaling (at fixed V_g): ave falls ~ 1/sqrt(n) at psi = 1, ext rises with n,
# fit falls fastest.
cell_mean <- function(v) {
  a <- aggregate(pv_run[[v]], pv_run[c("model", "psi", "n")],
                 function(x) c(m = mean(x), se = sd(x) / sqrt(length(x))))
  a <- do.call(data.frame, a); names(a)[4:5] <- c("m", "se"); a
}
zc <- cell_mean("z_sd"); rc <- cell_mean("zg_ratio")

p_zn <- ggplot(zc, aes(n, m, color = psi, group = psi)) +
  geom_errorbar(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se), width = 0.03) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.3) +
  facet_wrap(~ model) + scale_x_log10(breaks = sort(unique(zc$n))) + scale_y_log10() +
  scale_color_viridis_c(end = 0.9) +
  labs(x = "n (log scale)", y = "Phenotypic SD (log scale)", color = expression(psi),
       caption = "Mean over runs, 95% CI") +
  theme_bw()
ggsave(pv_path("z_sd_vs_n.pdf"), p_zn, width = 11, height = 4.5, dpi = 300)

p_rn <- ggplot(rc, aes(n, m, color = psi, group = psi)) +
  geom_line(data = theory_df, aes(y = sd_ratio_theory), linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se), width = 0.03) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.3) +
  facet_wrap(~ model) + scale_x_log10(breaks = sort(unique(rc$n))) + scale_y_log10() +
  scale_color_viridis_c(end = 0.9) +
  labs(x = "n (log scale)", y = "z_sd / g_sd (log scale)", color = expression(psi),
       caption = "Solid = observed mean over runs (95% CI); dashed = order-statistic theory at fixed V_g") +
  theme_bw()
ggsave(pv_path("zg_ratio_vs_n.pdf"), p_rn, width = 11, height = 4.5, dpi = 300)

# Effect size: % change per doubling of n (slope of log value on log2 n, lm over runs),
# for observed z_sd, observed z_sd/g_sd, and theory sqrt(c).
scale_list <- list()
for (mdl in models) for (ps in sort(unique(pv_run$psi))) {
  sub <- pv_run[pv_run$model == mdl & pv_run$psi == ps & pv_run$z_sd > 0, ]
  if (length(unique(sub$n)) < 2) next
  fz <- lm(log(z_sd) ~ log2(n), data = sub);     cz <- confint(fz)["log2(n)", ]
  fr <- lm(log(zg_ratio) ~ log2(n), data = sub); cr <- confint(fr)["log2(n)", ]
  th <- theory_df[theory_df$model == mdl & theory_df$psi == ps, ]
  ft <- lm(log(sd_ratio_theory) ~ log2(n), data = th)
  pct <- function(b) 100 * (exp(b) - 1)
  scale_list[[length(scale_list) + 1]] <- data.frame(
    model = mdl, psi = ps,
    z_sd_pct_per_doubling = pct(coef(fz)[[2]]), z_sd_lo = pct(cz[[1]]), z_sd_hi = pct(cz[[2]]),
    ratio_pct_per_doubling = pct(coef(fr)[[2]]), ratio_lo = pct(cr[[1]]), ratio_hi = pct(cr[[2]]),
    theory_pct_per_doubling = pct(coef(ft)[[2]]))
}
scale_df <- do.call(rbind, scale_list)
write.csv(scale_df, pv_path("z_sd_scaling_with_n.csv"), row.names = FALSE)
p_zn; p_rn


# ---- 7.3 decomposition: phenotypic SD / genetic SD ----
# Per run, z_sd / g_sd (same run, same window); cell = mean over runs with a
# bootstrap CI. < 1: the social rule shrinks genetic variation; > 1: amplifies it.
# E/D/? here means "no net change" vs "shrinks/amplifies" relative to 1.
dec_list <- list()
for (mdl in models) for (ps in sort(unique(pv_run$psi))) for (nn in sort(unique(pv_run$n))) {
  x <- pv_run$zg_ratio[pv_run$model == mdl & pv_run$psi == ps & pv_run$n == nn]
  if (length(x) < 2) next
  ci <- boot_mean_ci(x, n_boot)
  dec_list[[length(dec_list) + 1]] <- data.frame(
    model = mdl, psi = ps, n = nn, ratio = mean(x),
    ratio_lo90 = ci[["lo90"]], ratio_hi90 = ci[["hi90"]], ratio_lo = ci[["lo"]], ratio_hi = ci[["hi"]])
}
dec_df <- do.call(rbind, dec_list)
dec_df$class <- classify_equiv(dec_df$ratio_lo90, dec_df$ratio_hi90, equiv_margin)
dec_df <- merge(dec_df, theory_df[c("model", "psi", "n", "sd_ratio_theory")])
dec_df <- dec_df[order(dec_df$model, dec_df$psi, dec_df$n), ]
write.csv(dec_df, pv_path("zg_ratio_decomposition.csv"), row.names = FALSE)

dec_plots <- setNames(lapply(models, function(mdl)
  ratio_heatmap(dec_df[dec_df$model == mdl, ],
                paste0("Phenotypic SD / genetic SD (", mdl, ")"),
                fill_name = "z_sd / g_sd")), models)
save_by_model_pv(dec_plots, "zg_ratio_heatmap")


# ---- 7.4 analytical prediction of phenotypic variance ----
# Predict each run's phenotypic variance from its own V_g (mean g_sd^2 over the window):
#   V_z_pred = c(model, psi, n, d) * V_g, with d = (g_mean - optimum) / g_sd (fit only).
# ave also gets the formula (1-psi)^2 V_g + psi^2 V_g / n for comparison; it differs from
# the SLiM rule (different weighting, and it omits the 2 psi (1-psi) V_g / n covariance
# between a_i and the group mean).
# Assumes normally distributed genotypes; fit/ext predictions depend on the tails.
pv_run$d <- (pv_run$g_mean - optimum_burnin) / pv_run$g_sd
pv_run$factor_pred <- mapply(function(m, p, n, d) pheno_var_factor(as.character(m), p, n, d),
                             pv_run$model, pv_run$psi, pv_run$n, pv_run$d)
pv_run$z_var_pred <- pv_run$factor_pred * pv_run$g_var
pv_run$z_var_pred_simple <- ifelse(pv_run$model == "ave",
                                   ((1 - pv_run$psi)^2 + pv_run$psi^2 / pv_run$n) * pv_run$g_var, NA)
pv_run$obs_over_pred <- pv_run$z_var / pv_run$z_var_pred
write.csv(pv_run, pv_path("pv_per_run_with_predictions.csv"), row.names = FALSE)

p_pred <- ggplot(pv_run[pv_run$z_var > 0 & pv_run$z_var_pred > 0, ],
                 aes(z_var_pred, z_var, color = psi, shape = factor(n))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_wrap(~ model) + scale_x_log10() + scale_y_log10() +
  scale_color_viridis_c(end = 0.9) +
  labs(x = "Predicted phenotypic variance c * V_g (log)", y = "Observed phenotypic variance (log)",
       color = expression(psi), shape = "n", caption = "One point per run; dashed = 1:1") +
  theme_bw()
ggsave(pv_path("pred_vs_obs_z_var.pdf"), p_pred, width = 11, height = 4.5, dpi = 300)
p_pred

# Agreement per model and per cell: mean observed / predicted (bootstrap CI over runs),
# plus log-log slope and R^2 of observed on predicted across runs.
agree_df <- do.call(rbind, lapply(models, function(mdl) {
  s <- pv_run[pv_run$model == mdl, ]
  cbind(model = mdl, prediction = "SLiM rule", agree_one(s$z_var, s$z_var_pred))
}))
s_ave <- pv_run[pv_run$model == "ave", ]
agree_df <- rbind(agree_df, cbind(model = "ave", prediction = "(1-psi)^2 Vg + psi^2 Vg/n",
                                  agree_one(s_ave$z_var, s_ave$z_var_pred_simple)))
write.csv(agree_df, pv_path("prediction_agreement.csv"), row.names = FALSE)

cell_agree <- aggregate(obs_over_pred ~ model + psi + n, data = pv_run, FUN = mean)
write.csv(cell_agree, pv_path("prediction_obs_over_pred_by_cell.csv"), row.names = FALSE)


# ---- 7.5 phenotypic variance report ----
pv_rep <- c(
  "==================================================================",
  " PHENOTYPIC VARIANCE (burn-in, generation < 10000)",
  "==================================================================",
  sprintf(" Per-run value: mean over all time points. Unit of analysis: one run."),
  sprintf(" Bootstrap CIs resample runs (%d resamples). E/D/? uses the 90%% CI vs [%.2f, %.2f].",
          n_boot, equiv_margin[1], equiv_margin[2]),
  ""
)

pv_rep <- c(pv_rep, rule, " 1. Phenotypic SD relative to psi = 0 (grid rows psi, columns n)", rule)
for (mdl in models) {
  d <- pz0_df[pz0_df$model == mdl, ]; if (!nrow(d)) next
  cnt <- table(factor(d$class, levels = names(class_code)))
  pv_rep <- c(pv_rep, "", sprintf(" %s: equivalent %d, different %d, inconclusive %d (of %d cells)",
                                  mdl, cnt[["equivalent"]], cnt[["different"]], cnt[["inconclusive"]], nrow(d)),
              grid_txt(d))
}

pv_rep <- c(pv_rep, "", rule, " 2. Scaling with n: % change per doubling of n",
            "    z_sd = observed [95% CI]; z/g = observed z_sd/g_sd [95% CI]; theory = sqrt(c) at fixed V_g",
            "    expected: ave ~ -29% per doubling at psi = 1 (1/sqrt(n)); ext rising; fit falling fastest", rule)
for (i in seq_len(nrow(scale_df))) {
  t <- scale_df[i, ]
  pv_rep <- c(pv_rep, sprintf("   %-4s psi = %-5s z_sd %+6.1f%% [%+.1f, %+.1f]   z/g %+6.1f%% [%+.1f, %+.1f]   theory %+6.1f%%",
                              t$model, t$psi, t$z_sd_pct_per_doubling, t$z_sd_lo, t$z_sd_hi,
                              t$ratio_pct_per_doubling, t$ratio_lo, t$ratio_hi, t$theory_pct_per_doubling))
}

pv_rep <- c(pv_rep, "", rule, " 3. Decomposition: z_sd / g_sd per cell (<1 shrinks, >1 amplifies)",
            "    grid: observed ratio + E/D/? vs 1; second grid: order-statistic theory sqrt(c)", rule)
for (mdl in models) {
  d <- dec_df[dec_df$model == mdl, ]; if (!nrow(d)) next
  cnt <- table(factor(d$class, levels = names(class_code)))
  th <- d; th$ratio <- th$sd_ratio_theory
  pv_rep <- c(pv_rep, "", sprintf(" %s: equivalent to 1: %d, different: %d, inconclusive: %d (of %d cells)",
                                  mdl, cnt[["equivalent"]], cnt[["different"]], cnt[["inconclusive"]], nrow(d)),
              "   observed:", grid_txt(d), "   theory:", grid_txt(th, "%.2f", with_class = FALSE))
}

pv_rep <- c(pv_rep, "", rule, " 4. Prediction of phenotypic variance from each run's V_g",
            "    obs/pred = mean over runs of observed V_z / predicted V_z [95% CI]; 1 = perfect", rule)
for (i in seq_len(nrow(agree_df))) {
  t <- agree_df[i, ]
  pv_rep <- c(pv_rep, sprintf("   %-4s %-28s obs/pred = %s; log-log slope = %s, R^2 = %s (n = %d runs)",
                              t$model, t$prediction, ci3(t$obs_over_pred, t$lo, t$hi),
                              f3(t$loglog_slope), f3(t$r_squared), t$n_runs))
}
for (mdl in models) {
  d <- cell_agree[cell_agree$model == mdl, ]; if (!nrow(d)) next
  d$ratio <- d$obs_over_pred
  pv_rep <- c(pv_rep, "", sprintf("   %s obs/pred by cell (rows psi, columns n):", mdl),
              grid_txt(d, "%.2f", with_class = FALSE))
}

cat(pv_rep, sep = "\n")
writeLines(pv_rep, pv_path("phenotypic_variance_report.txt"))
