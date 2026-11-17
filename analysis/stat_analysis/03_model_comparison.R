# 3. model comparison on per-run values
# Section 3 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/03_model_comparison.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
source("analysis/stat_analysis/00_load_data.R")
if (!exists("run_df")) source("analysis/stat_analysis/02_per_run_summary.R")   # needs objects from section 2

# Comparisons are made within each psi/n parameter set, so model differences are
# not confounded with parameter differences.
metrics    <- c("mean_sd_all", "mean_mean_all")
sd_metrics <- c("mean_sd_all")   # bootstrap diff/ratio and CV are for SD only
model_pairs <- list(c("fit", "ave"), c("fit", "ext"), c("ave", "ext"))  # each pair is (A, B)

group_list <- list(); anova_list <- list(); pair_list <- list()

for (ps in unique(run_df$param_set)) {
  sub <- run_df[run_df$param_set == ps, ]
  sub$model <- factor(sub$model, levels = intersect(models, sub$model))

  for (m in metrics) {
    # (a) group means, SDs, n and, for SD metrics, CV of the per-run SD
    for (g in levels(sub$model)) {
      x <- sub[[m]][sub$model == g]
      group_list[[length(group_list) + 1]] <- data.frame(
        param_set = ps, metric = m, model = g, n_runs = length(x),
        mean = mean(x), sd = sd(x),
        cv = if (m %in% sd_metrics) sd(x) / mean(x) else NA
      )
    }

    # (b) one-way ANOVA across models (runs as replicates)
    if (nlevels(sub$model) >= 2) {
      tab <- summary(aov(sub[[m]] ~ sub$model))[[1]]
      anova_list[[length(anova_list) + 1]] <- data.frame(
        param_set = ps, metric = m,
        F = tab[1, "F value"], df_model = tab[1, "Df"], df_resid = tab[2, "Df"],
        p = tab[1, "Pr(>F)"]
      )
    }

    # (c) pairwise Welch t-tests, plus bootstrap diff/ratio for SD metrics
    for (pr in model_pairs) {
      a <- sub[[m]][sub$model == pr[1]]
      b <- sub[[m]][sub$model == pr[2]]
      if (length(a) < 2 || length(b) < 2) next
      tt <- tryCatch(t.test(b, a, var.equal = FALSE), error = function(e) NULL)
      row <- data.frame(
        param_set = ps, metric = m, A = pr[1], B = pr[2],
        mean_A = mean(a), mean_B = mean(b),
        diff = mean(b) - mean(a), diff_lo = NA, diff_hi = NA,
        ratio = NA, ratio_lo = NA, ratio_hi = NA, ratio_lo90 = NA, ratio_hi90 = NA,
        welch_t  = if (is.null(tt)) NA else unname(tt$statistic),
        welch_df = if (is.null(tt)) NA else unname(tt$parameter),
        welch_p  = if (is.null(tt)) NA else tt$p.value
      )
      if (m %in% sd_metrics) {
        ci <- boot_diff_ratio(a, b, n_boot)
        row$ratio <- mean(b) / mean(a)
        row[names(ci)] <- ci
      }
      pair_list[[length(pair_list) + 1]] <- row
    }
  }
}

group_stats_df <- do.call(rbind, group_list)
anova_df       <- do.call(rbind, anova_list)
pairwise_df    <- do.call(rbind, pair_list)

# Many tests are run (parameter sets x metrics x pairs): add Holm-adjusted p-values
# within each metric so the family-wise error rate is controlled.
anova_df$p_holm       <- ave(anova_df$p,          anova_df$metric,    FUN = function(p) p.adjust(p, "holm"))
pairwise_df$welch_p_holm <- ave(pairwise_df$welch_p, pairwise_df$metric, FUN = function(p) p.adjust(p, "holm"))

write.csv(group_stats_df, out_path("model_group_stats.csv"), row.names = FALSE)
write.csv(anova_df,       out_path("model_anova.csv"),       row.names = FALSE)
write.csv(pairwise_df,    out_path("model_pairwise.csv"),    row.names = FALSE)
