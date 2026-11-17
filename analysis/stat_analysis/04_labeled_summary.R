# 4. clean, labeled summary
# Section 4 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/04_labeled_summary.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
source("analysis/stat_analysis/00_load_data.R")
if (!exists("pairwise_df")) source("analysis/stat_analysis/03_model_comparison.R")   # needs objects from section 3

report <- c(
  "==================================================================",
  paste0(" MODEL COMPARISON (", trait_word, "): fit vs ave vs ext (burn-in, generation < 10000)"),
  "==================================================================",
  " Unit of analysis: one run (psi, n, reg, rep). Each run is reduced to",
  " one value per metric before any test; generations are never treated",
  " as independent observations.",
  sprintf(" Bootstrap: %d resamples of runs within each model, percentile 95%% CI.", n_boot),
  " Pairwise contrasts are B vs A: diff = mean(B) - mean(A), ratio = mean(B)/mean(A).",
  " p_holm = Holm-adjusted across all parameter sets/pairs for that metric.",
  ""
)

for (ps in unique(run_df$param_set)) {
  report <- c(report, "------------------------------------------------------------------",
              paste0(" PARAMETER SET: ", ps),
              "------------------------------------------------------------------")
  for (m in metrics) {
    report <- c(report, "", paste0(" [", metric_labels[[m]], "]"))

    # group stats
    gs <- group_stats_df[group_stats_df$param_set == ps & group_stats_df$metric == m, ]
    for (i in seq_len(nrow(gs))) {
      report <- c(report, sprintf("   %-4s n = %3d   mean = %-10s sd = %-10s%s",
                                  gs$model[i], gs$n_runs[i], fmt(gs$mean[i]), fmt(gs$sd[i]),
                                  if (m %in% sd_metrics) paste0(" CV = ", fmt(gs$cv[i])) else ""))
    }

    # ANOVA
    an <- anova_df[anova_df$param_set == ps & anova_df$metric == m, ]
    if (nrow(an) == 1) {
      report <- c(report, sprintf("   One-way ANOVA: F(%d, %d) = %s, p %s (p_holm %s)",
                                  an$df_model, an$df_resid, fmt(an$F), fmt_p(an$p), fmt_p(an$p_holm)))
    }

    # pairwise
    pw <- pairwise_df[pairwise_df$param_set == ps & pairwise_df$metric == m, ]
    for (i in seq_len(nrow(pw))) {
      line <- sprintf("   %s vs %s: diff = %s", pw$B[i], pw$A[i], fmt(pw$diff[i]))
      if (m %in% sd_metrics) {
        line <- paste0(line, sprintf(" [%s, %s]; ratio = %s [%s, %s]",
                                     fmt(pw$diff_lo[i]), fmt(pw$diff_hi[i]),
                                     fmt(pw$ratio[i]), fmt(pw$ratio_lo[i]), fmt(pw$ratio_hi[i])))
      }
      line <- paste0(line, sprintf("; Welch t(%s) = %s, p %s (p_holm %s)",
                                   fmt(pw$welch_df[i], 3), fmt(pw$welch_t[i]),
                                   fmt_p(pw$welch_p[i]), fmt_p(pw$welch_p_holm[i])))
      report <- c(report, line)
    }
  }
  report <- c(report, "")
}

cat(report, sep = "\n")
writeLines(report, out_path("model_comparison_summary.txt"))
