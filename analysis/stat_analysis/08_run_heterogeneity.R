# 8. run-to-run heterogeneity
# Section 8 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/08_run_heterogeneity.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
source("analysis/stat_analysis/00_load_data.R")
if (!exists("pv_run") || !exists("pv_df")) source("analysis/stat_analysis/07_phenotypic_variance.R")   # needs objects from section 7

# CV across runs of the per-run SD (sd / mean over the runs in a psi x n cell), for
# phenotypic (z_sd) and genetic (g_sd) SD, using the per-run values from section 7.
# High CV marks unstable regions where cell averages may hide diverging runs; those
# cells are flagged and their individual trajectories plotted for inspection.
het_folder <- file.path(base_plot_folder, "heterogeneity")
dir.create(het_folder, recursive = TRUE, showWarnings = FALSE)
het_path <- function(name) file.path(het_folder, name)

cv_flag           <- 0.2   # flag cells whose CV of per-run SD exceeds this
max_inspect_cells <- 24    # at most this many flagged cells (highest CV first) get trajectory plots
het_vars <- c(z_sd = "Phenotypic SD", g_sd = "Genetic SD")

cv_list <- list()
for (v in names(het_vars)) for (mdl in models) for (ps in sort(unique(pv_run$psi))) for (nn in sort(unique(pv_run$n))) {
  x <- pv_run[[v]][pv_run$model == mdl & pv_run$psi == ps & pv_run$n == nn]
  if (length(x) < 3) next
  # bootstrap the CV over runs for a 95% CI
  bx <- matrix(x[sample.int(length(x), length(x) * n_boot, replace = TRUE)], nrow = length(x))
  bcv <- apply(bx, 2, sd) / colMeans(bx)
  # the most atypical run: largest |log(run SD / cell median)|
  dev <- abs(log(x / median(x)))
  reps <- pv_run$rep[pv_run$model == mdl & pv_run$psi == ps & pv_run$n == nn]
  cv_list[[length(cv_list) + 1]] <- data.frame(
    variable = v, model = mdl, psi = ps, n = nn, n_runs = length(x),
    mean_sd = mean(x), cv = sd(x) / mean(x),
    cv_lo = quantile(bcv, 0.025, names = FALSE), cv_hi = quantile(bcv, 0.975, names = FALSE),
    min_sd = min(x), max_sd = max(x),
    most_atypical_rep = reps[which.max(dev)], most_atypical_fold = exp(max(dev)))
}
cv_df <- do.call(rbind, cv_list)
cv_df$flag <- cv_df$cv > cv_flag
write.csv(cv_df, het_path("cv_of_per_run_sd.csv"), row.names = FALSE)

# Heatmap of CV, psi x n, one file per model and variable. The color scale is shared
# across models (per variable) so that CVs can be compared between models;
# "!" marks flagged cells.
plot_cv_heatmap <- function(d, title, lims) {
  d$label <- sprintf("%.2f%s", d$cv, ifelse(d$flag, "\n!", ""))
  ggplot(d, aes(x = factor(n), y = factor(psi), fill = cv)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.6, lineheight = 0.85, color = heat_label_col) +
    scale_fill_gradient(low = "#f0efec", high = "#8073ac", limits = lims, trans = "sqrt",
                        name = "CV across runs") +
    labs(x = "n", y = expression(psi), title = title,
         caption = sprintf("CV = sd / mean of per-run SD over runs (all time points); ! = CV > %.2f",
                           cv_flag)) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank())
}
for (v in names(het_vars)) {
  dv <- cv_df[cv_df$variable == v, ]
  lims <- c(0, max(dv$cv, cv_flag, na.rm = TRUE))
  plots <- setNames(lapply(models, function(mdl)
    plot_cv_heatmap(dv[dv$model == mdl, ], paste0(het_vars[[v]], ": run-to-run CV (", mdl, ")"), lims)),
    models)
  for (mdl in models) ggsave(het_path(paste0("cv_heatmap_", v, "_", mdl, ".pdf")), plots[[mdl]],
                             width = 7, height = 6, dpi = 300)
  pdf(het_path(paste0("cv_heatmap_", v, "_by_model.pdf")), width = 7, height = 6)
  for (p in plots) print(p)
  dev.off()
}

# Trajectories of flagged cells: one line per run over the burn-in (log scale),
# so diverging or non-stationary runs can be spotted before trusting cell averages.
for (v in names(het_vars)) {
  fl <- cv_df[cv_df$variable == v & cv_df$flag, ]
  if (!nrow(fl)) next
  fl <- head(fl[order(-fl$cv), ], max_inspect_cells)
  fl$cell <- sprintf("%s psi=%s n=%s (CV %.2f)", fl$model, fl$psi, fl$n, fl$cv)
  tr <- merge(pv_df, data.frame(reg = fl$model, psi = fl$psi, n = fl$n, cell = fl$cell))
  tr$cell <- factor(tr$cell, levels = fl$cell)
  tr <- tr[tr[[v]] > 0, ]
  p <- ggplot(tr, aes(tick, .data[[v]], group = rep, color = factor(rep))) +
    geom_line(linewidth = 0.3, alpha = 0.8) +
    facet_wrap(~ cell, scales = "free_y") + scale_y_log10() +
    labs(x = "Generation", y = paste(het_vars[[v]], "(log scale)"), color = "rep",
         title = paste("Flagged cells:", het_vars[[v]], "trajectories per run")) +
    theme_bw(base_size = 9) + theme(legend.position = "none")
  ncols <- min(4, nrow(fl))
  ggsave(het_path(paste0("flagged_trajectories_", v, ".pdf")), p,
         width = 3 * ncols, height = 2.5 * ceiling(nrow(fl) / ncols) + 0.6, dpi = 200, limitsize = FALSE)
}

# Text report
het_rep <- c(
  "==================================================================",
  " RUN-TO-RUN HETEROGENEITY (burn-in, generation < 10000)",
  "==================================================================",
  sprintf(" CV = sd / mean over runs of the per-run SD (all time points); 95%% bootstrap CI over runs."),
  " The CV CI is rough: when a few runs drive the CV, resamples that miss them give a",
  " very low CV, so the interval is wide and skewed low. Inspect the trajectories instead.",
  sprintf(" Flag: CV > %.2f. Atypical = run with the largest fold difference from the cell median.", cv_flag),
  ""
)
for (v in names(het_vars)) {
  het_rep <- c(het_rep, rule, paste0(" ", het_vars[[v]], " (", v, ")"), rule)
  for (mdl in models) {
    d <- cv_df[cv_df$variable == v & cv_df$model == mdl, ]
    if (!nrow(d)) next
    g <- tapply(sprintf("%.2f%s", d$cv, ifelse(d$flag, "!", "")), list(d$psi, d$n), identity)
    g[is.na(g)] <- "-"
    het_rep <- c(het_rep, "", sprintf(" %s: %d of %d cells flagged; CV grid (rows psi, columns n):",
                                      mdl, sum(d$flag), nrow(d)),
                 paste0("     ", capture.output(print(noquote(g)))))
    fl <- d[d$flag, ]
    fl <- fl[order(-fl$cv), ]
    if (nrow(fl)) {
      het_rep <- c(het_rep, "   Flagged cells (highest CV first):",
                   sprintf("     psi = %-5s n = %-4s CV = %s  SD range [%s, %s]  most atypical: rep %s (%.2fx median)",
                           fl$psi, fl$n, ci3(fl$cv, fl$cv_lo, fl$cv_hi), fmt(fl$min_sd), fmt(fl$max_sd),
                           fl$most_atypical_rep, fl$most_atypical_fold))
    }
  }
  het_rep <- c(het_rep, "")
}
cat(het_rep, sep = "\n")
writeLines(het_rep, het_path("heterogeneity_report.txt"))
