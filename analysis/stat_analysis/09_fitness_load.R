# 9. mean fitness and genetic load
# Section 9 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/09_fitness_load.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
source("analysis/stat_analysis/00_load_data.R")

# Mean population fitness across psi x n, and how it follows from the variance results.
# SLiM fitness: w = dnorm(P, optimum, sqrt(Vs)) with Vs = 2N, so relative fitness is
# w / w_max with w_max = 1 / sqrt(2 pi Vs), and load L = 1 - mean relative fitness.
# Gaussian prediction for phenotypes P ~ N(zbar, V_z):
#   W/W_max = sqrt(Vs / (Vs + V_z)) * exp(-(zbar - optimum)^2 / (2 (Vs + V_z)))
#   variance load L_var = 1 - sqrt(Vs / (Vs + V_z)); the rest is load from the mean offset.
# For small V_z / Vs, L ~ (V_z + (zbar - optimum)^2) / (2 Vs) whatever the shape of P.
fl_folder <- file.path(base_plot_folder, "fitness_load")
dir.create(fl_folder, recursive = TRUE, showWarnings = FALSE)
fl_path <- function(name) file.path(fl_folder, name)

Vs    <- 2000                   # SLiM: Vs = 2 * N
w_max <- 1 / sqrt(2 * pi * Vs)
if (!"fitness_mean" %in% names(cons_data_df)) stop("cons_data_df has no fitness_mean column")

fl_df <- cons_data_df[cons_data_df[[tick_col]] < burnin_end & cons_data_df$reg %in% models,
                      c(tick_col, "psi", "n", "reg", "rep", "shift_size", "fitness_mean", "z_mean", "z_sd")]
names(fl_df)[1] <- "tick"
fl_key <- interaction(fl_df$psi, fl_df$n, fl_df$reg, fl_df$rep, drop = TRUE)
fl_df  <- fl_df[fl_df$shift_size == ave(fl_df$shift_size, fl_key, FUN = function(x) min(x, na.rm = TRUE)), ]
fl_df  <- fl_df[!duplicated(fl_df[c("psi", "n", "reg", "rep", "tick")]), ]
# abs_fit is only filled in SLiM from shift_time - 50; earlier values are 0 = not recorded.
n_unset <- sum(!(fl_df$fitness_mean > 0), na.rm = TRUE) + sum(is.na(fl_df$fitness_mean))
if (n_unset > 0) warning(n_unset, " time points with fitness_mean <= 0 or NA were dropped (fitness not recorded)")
fl_df  <- fl_df[!is.na(fl_df$fitness_mean) & fl_df$fitness_mean > 0, ]
fl_df  <- fl_df[order(fl_df$reg, fl_df$psi, fl_df$n, fl_df$rep, fl_df$tick), ]

# Per time point: observed relative fitness, and the Gaussian prediction from z_mean, z_sd.
fl_df$w_rel      <- fl_df$fitness_mean / w_max
fl_df$w_pred     <- sqrt(Vs / (Vs + fl_df$z_sd^2)) *
                    exp(-(fl_df$z_mean - optimum_burnin)^2 / (2 * (Vs + fl_df$z_sd^2)))
fl_df$w_pred_var <- sqrt(Vs / (Vs + fl_df$z_sd^2))   # variance part only

# Per run: average over all time points, as in sections 6-8.
fl_run <- do.call(rbind, lapply(split(fl_df, list(fl_df$reg, fl_df$psi, fl_df$n, fl_df$rep), drop = TRUE),
  function(d) {
    data.frame(model = d$reg[1], psi = d$psi[1], n = d$n[1], rep = d$rep[1], n_timepoints = nrow(d),
               w_rel = mean(d$w_rel), load = 1 - mean(d$w_rel),
               load_pred = 1 - mean(d$w_pred), load_pred_var = 1 - mean(d$w_pred_var))
  }))
fl_run$model <- factor(fl_run$model, levels = models)
fl_run$load_pred_mean <- fl_run$load_pred - fl_run$load_pred_var     # load from the mean offset
rownames(fl_run) <- NULL
write.csv(fl_run, fl_path("fitness_per_run.csv"), row.names = FALSE)

# Cell summaries: mean over runs with a 95% bootstrap CI.
fl_cell_list <- list()
for (mdl in models) for (ps in sort(unique(fl_run$psi))) for (nn in sort(unique(fl_run$n))) {
  s <- fl_run[fl_run$model == mdl & fl_run$psi == ps & fl_run$n == nn, ]
  if (nrow(s) < 2) next
  cw <- boot_mean_ci(s$w_rel, n_boot); cl <- boot_mean_ci(s$load, n_boot)
  fl_cell_list[[length(fl_cell_list) + 1]] <- data.frame(
    model = mdl, psi = ps, n = nn, n_runs = nrow(s),
    w_rel = mean(s$w_rel), w_rel_lo = cw[["lo"]], w_rel_hi = cw[["hi"]],
    load = mean(s$load), load_lo = cl[["lo"]], load_hi = cl[["hi"]],
    load_pred = mean(s$load_pred),
    share_var = mean(s$load_pred_var) / mean(s$load_pred))   # fraction of predicted load from variance
}
fl_cell <- do.call(rbind, fl_cell_list)
fl_cell$model <- factor(fl_cell$model, levels = models)
write.csv(fl_cell, fl_path("fitness_load_by_cell.csv"), row.names = FALSE)


# ---- 9.1 plots: mean fitness and load across psi x n ----
# (a) Heatmap of mean relative fitness, psi x n, one file per model. Fitness is close to 1
#     everywhere, so cells are colored by load (1 - fitness) on a log scale, shared across
#     models (darker = lower fitness); labels give both fitness and load.
plot_w_heatmap <- function(d, title, lims) {
  ggplot(d, aes(x = factor(n), y = factor(psi), fill = load)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.4f\nL = %.2g", w_rel, load)), size = 2.3, lineheight = 0.85,
              color = heat_label_col) +
    scale_fill_gradient(low = "#f0efec", high = "#5aae61", limits = lims, trans = "log10",
                        name = "Load\n(1 - fitness)") +
    labs(x = "n", y = expression(psi), title = title,
         caption = sprintf("Label: mean relative fitness (w / w_max, Vs = %g) and load L; all time points",
                           Vs)) +
    theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
}
w_lims <- range(fl_cell$load[fl_cell$load > 0], na.rm = TRUE)
w_plots <- setNames(lapply(models, function(mdl)
  plot_w_heatmap(fl_cell[fl_cell$model == mdl, ], paste0("Mean relative fitness (", mdl, ")"), w_lims)), models)
for (mdl in models) ggsave(fl_path(paste0("fitness_heatmap_", mdl, ".pdf")), w_plots[[mdl]],
                           width = 7, height = 6, dpi = 300)
pdf(fl_path("fitness_heatmap_by_model.pdf"), width = 7, height = 6)
for (p in w_plots) print(p)
dev.off()

# (b) Load vs psi, one line per n, one panel per model (log y; 95% CI over runs).
p_load <- ggplot(fl_cell, aes(psi, load, color = factor(n), group = factor(n))) +
  geom_errorbar(aes(ymin = load_lo, ymax = load_hi), width = 0.01) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
  facet_wrap(~ model) + scale_y_log10() + scale_color_viridis_d(end = 0.9) +
  labs(x = expression(psi), y = "Load = 1 - mean relative fitness (log)", color = "n",
       caption = "Mean over runs, 95% bootstrap CI") +
  theme_bw()
ggsave(fl_path("load_vs_psi.pdf"), p_load, width = 11, height = 4.5, dpi = 300)

# (c) Load relative to psi = 0 (bootstrap over runs), E/D/? labels, one file per model.
lr_list <- list()
for (mdl in models) for (nn in sort(unique(fl_run$n))) {
  s <- fl_run[fl_run$model == mdl & fl_run$n == nn, ]
  ref <- s$load[s$psi == 0]
  if (length(ref) < 2 || mean(ref) <= 0) next
  for (ps in setdiff(sort(unique(s$psi)), 0)) {
    x <- s$load[s$psi == ps]
    if (length(x) < 2) next
    ci <- boot_diff_ratio(ref, x, n_boot)
    lr_list[[length(lr_list) + 1]] <- data.frame(
      model = mdl, n = nn, psi = ps, ratio = mean(x) / mean(ref),
      ratio_lo90 = ci[["ratio_lo90"]], ratio_hi90 = ci[["ratio_hi90"]],
      ratio_lo   = ci[["ratio_lo"]],   ratio_hi   = ci[["ratio_hi"]])
  }
}
lr_df <- do.call(rbind, lr_list)
lr_df$class <- classify_equiv(lr_df$ratio_lo90, lr_df$ratio_hi90, equiv_margin)
write.csv(lr_df, fl_path("load_ratio_vs_psi0.csv"), row.names = FALSE)
lr_plots <- setNames(lapply(models, function(mdl)
  ratio_heatmap(lr_df[lr_df$model == mdl, ], paste0("Load relative to psi = 0 (", mdl, ")"),
                fill_name = "load ratio")), models)
for (mdl in models) ggsave(fl_path(paste0("load_ratio_vs_psi0_heatmap_", mdl, ".pdf")), lr_plots[[mdl]],
                           width = 7, height = 6, dpi = 300)
pdf(fl_path("load_ratio_vs_psi0_heatmap_by_model.pdf"), width = 7, height = 6)
for (p in lr_plots) print(p)
dev.off()

# (d) Link to variance: observed load vs the Gaussian prediction from z_mean and z_sd
#     (one point per run; dashed = 1:1).
p_link <- ggplot(fl_run[fl_run$load > 0 & fl_run$load_pred > 0, ],
                 aes(load_pred, load, color = psi, shape = factor(n))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_wrap(~ model) + scale_x_log10() + scale_y_log10() + scale_color_viridis_c(end = 0.9) +
  labs(x = "Predicted load from phenotypic mean and variance (log)", y = "Observed load (log)",
       color = expression(psi), shape = "n", caption = "One point per run; dashed = 1:1") +
  theme_bw()
ggsave(fl_path("load_obs_vs_pred.pdf"), p_link, width = 11, height = 4.5, dpi = 300)
p_load; p_link


# ---- 9.2 fitness and load report ----
fl_agree <- do.call(rbind, lapply(models, function(mdl) {
  s <- fl_run[fl_run$model == mdl, ]
  cbind(model = mdl, agree_one(s$load, s$load_pred))
}))
write.csv(fl_agree, fl_path("load_prediction_agreement.csv"), row.names = FALSE)

cell_grid <- function(d, v, f = "%.4f") {
  g <- tapply(sprintf(f, d[[v]]), list(d$psi, d$n), identity); g[is.na(g)] <- "-"
  paste0("     ", capture.output(print(noquote(g))))
}

fl_rep <- c(
  "==================================================================",
  " MEAN FITNESS AND GENETIC LOAD (burn-in, generation < 10000)",
  "==================================================================",
  sprintf(" Relative fitness = mean w / w_max, w_max = 1/sqrt(2 pi Vs), Vs = %g; load = 1 - relative fitness.", Vs),
  sprintf(" Per-run value: mean over all time points; time points with fitness not recorded dropped."),
  " Unit of analysis: one run; CIs are bootstrap over runs.",
  ""
)
for (mdl in models) {
  d <- fl_cell[fl_cell$model == mdl, ]; if (!nrow(d)) next
  fl_rep <- c(fl_rep, rule, paste0(" ", mdl), rule,
              " Mean relative fitness (rows psi, columns n):", cell_grid(d, "w_rel"),
              " Load (rows psi, columns n):",                  cell_grid(d, "load"),
              " Share of predicted load due to phenotypic variance (rest = mean offset):",
              cell_grid(d, "share_var", "%.2f"))
  dl <- lr_df[lr_df$model == mdl, ]
  if (nrow(dl)) {
    cnt <- table(factor(dl$class, levels = names(class_code)))
    fl_rep <- c(fl_rep, sprintf(" Load relative to psi = 0: equivalent %d, different %d, inconclusive %d (of %d cells)",
                                cnt[["equivalent"]], cnt[["different"]], cnt[["inconclusive"]], nrow(dl)),
                grid_txt(dl))
  }
  t <- fl_agree[fl_agree$model == mdl, ]
  fl_rep <- c(fl_rep, sprintf(" Observed / Gaussian-predicted load = %s; log-log slope = %s, R^2 = %s (%d runs)",
                              ci3(t$obs_over_pred, t$lo, t$hi), f3(t$loglog_slope), f3(t$r_squared), t$n_runs),
              "")
}
cat(fl_rep, sep = "\n")
writeLines(fl_rep, fl_path("fitness_load_report.txt"))
