# 5. effect-size report: SD ratio heatmaps + equivalence classification
# Section 5 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/05_effect_size_report.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
source("analysis/stat_analysis/00_load_data.R")
if (!exists("pairwise_df")) source("analysis/stat_analysis/03_model_comparison.R")   # needs objects from section 3

heat_pairs   <- c("ave / fit", "ext / fit")   # panels shown in the heatmap and narrative

# One row per parameter set x SD metric x pair, with psi and n as numbers for plotting.
ratio_df <- pairwise_df[pairwise_df$metric %in% sd_metrics, ]
ratio_df <- merge(ratio_df, unique(run_df[c("param_set", "psi", "n")]), by = "param_set")
ratio_df$pair  <- paste(ratio_df$B, "/", ratio_df$A)
ratio_df$class <- classify_equiv(ratio_df$ratio_lo90, ratio_df$ratio_hi90, equiv_margin)
ratio_df <- ratio_df[order(ratio_df$metric, ratio_df$pair, ratio_df$psi, ratio_df$n), ]
write.csv(ratio_df[c("metric", "pair", "psi", "n", "ratio", "ratio_lo90", "ratio_hi90",
                     "ratio_lo", "ratio_hi", "class")],
          out_path("sd_ratio_equivalence.csv"), row.names = FALSE)

# Heatmap: psi x n, one panel per pair. Fill is log2(ratio) on a diverging scale
# centered on 1, so 0.9 and 1.11 get equally strong color. Each cell is labeled with
# the ratio and its class (E = equivalent, D = different, ? = inconclusive).
plot_ratio_heatmap <- function(m, file) {
  d <- ratio_df[ratio_df$metric == m & ratio_df$pair %in% heat_pairs, ]
  d$label <- sprintf("%.2f\n%s", d$ratio, class_code[d$class])
  lim <- max(abs(log2(d$ratio)), log2(1.10), na.rm = TRUE)
  brk <- c(0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.25)
  brk <- brk[abs(log2(brk)) <= lim]
  p <- ggplot(d, aes(x = factor(n), y = factor(psi), fill = log2(ratio))) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.3, lineheight = 0.85, color = heat_label_col) +
    facet_wrap(~ pair) +
    scale_fill_gradient2(low = "#4575b4", mid = "#f0efec", high = "#d73027", midpoint = 0,
                         limits = c(-lim, lim), breaks = log2(brk), labels = brk,
                         name = "SD ratio") +
    labs(x = "n", y = expression(psi),
         title = paste0("Ratio of ", trait_word, " SD between models (",
                        sub(".*SD, averaged", "averaged", metric_labels[[m]]), ")"),
         caption = sprintf(paste("E = equivalent, D = different, ? = inconclusive;",
                                 "90%% bootstrap CI vs margin [%.2f, %.2f]"),
                           equiv_margin[1], equiv_margin[2])) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank())
  ggsave(file, p, width = 11, height = 6, dpi = 300)
  p
}
for (m in sd_metrics) {
  plot_ratio_heatmap(m, out_path(paste0("sd_ratio_heatmap_", m, ".pdf")))
}

# Separate version: one heatmap per pair, each with its own color bar, because the
# pairs span very different ranges (ext/fit up to ~20x, ave/fit under ~2x); a shared
# scale washes out the ave/fit pattern. Each pair gets its own PNG, and all pairs
# for a metric are also written as pages of one PDF.

plot_ratio_heatmap_single <- function(m, pr) {
  d <- ratio_df[ratio_df$metric == m & ratio_df$pair == pr, ]
  ratio_heatmap(d, paste0(Trait_word, " SD ratio ", pr, " (",
                          sub(".*SD, averaged", "averaged", metric_labels[[m]]), ")"))
}

for (m in sd_metrics) {
  plots <- lapply(heat_pairs, function(pr) plot_ratio_heatmap_single(m, pr))
  for (i in seq_along(heat_pairs)) {
    pair_tag <- gsub(" / ", "_over_", heat_pairs[i])
    ggsave(out_path(paste0("sd_ratio_heatmap_", m, "_", pair_tag, ".pdf")),
           plots[[i]], width = 7, height = 6, dpi = 300)
  }
  pdf(out_path(paste0("sd_ratio_heatmap_", m, "_by_pair.pdf")), width = 7, height = 6)
  for (p in plots) print(p)   # one page per pair
  dev.off()
}

# Text report: effect sizes and equivalence first, Holm-corrected tests as a supplement.
fmt_ci <- function(r, lo, hi) sprintf("%.3f [%.3f, %.3f]", r, lo, hi)

eff <- c(
  "==================================================================",
  paste0(" ", toupper(trait_word), " SD RATIOS BETWEEN MODELS (burn-in, generation < 10000)"),
  "==================================================================",
  paste0(" Effect size: ratio of mean per-run ", trait_word, " SD, B / A, within each (psi, n)."),
  " Unit of analysis: one run; CIs from bootstrapping runs within each model",
  sprintf(" (%d resamples).", n_boot),
  sprintf(" Equivalence margin: SD ratio in [%.2f, %.2f].", equiv_margin[1], equiv_margin[2]),
  " Classification uses the 90% CI (TOST, alpha = 0.05):",
  "   equivalent   = 90% CI fully inside the margin",
  "   different    = 90% CI fully outside the margin",
  "   inconclusive = otherwise",
  " Format: ratio [90% CI] [95% CI].",
  ""
)

for (m in sd_metrics) {
  eff <- c(eff, "------------------------------------------------------------------",
           paste0(" ", metric_labels[[m]]),
           "------------------------------------------------------------------")
  for (pr in heat_pairs) {
    d <- ratio_df[ratio_df$metric == m & ratio_df$pair == pr, ]
    cnt <- table(factor(d$class, levels = names(class_code)))
    eff <- c(eff, "", sprintf(" %s: equivalent %d, different %d, inconclusive %d (of %d parameter sets)",
                              pr, cnt[["equivalent"]], cnt[["different"]], cnt[["inconclusive"]], nrow(d)))

    # Parameter sets that are not equivalent, largest departure from 1 first.
    ex <- d[d$class != "equivalent", ]
    ex <- ex[order(-abs(log(ex$ratio))), ]
    if (nrow(ex) > 0) {
      eff <- c(eff, "   Not equivalent (largest effect first):")
      eff <- c(eff, sprintf("     psi = %-5s n = %-4s %s %s  %s",
                            ex$psi, ex$n,
                            fmt_ci(ex$ratio, ex$ratio_lo90, ex$ratio_hi90),
                            sprintf("[%.3f, %.3f]", ex$ratio_lo, ex$ratio_hi),
                            ex$class))
    }

    # Compact psi x n grid of ratio + class code, same layout as the heatmap.
    grid <- tapply(sprintf("%.2f%s", d$ratio, class_code[d$class]), list(d$psi, d$n), identity)
    grid[is.na(grid)] <- "-"
    eff <- c(eff, "   Grid (rows psi, columns n; E/D/? = class):",
             paste0("     ", capture.output(print(noquote(grid)))))
  }
  eff <- c(eff, "")
}

# Supplementary: formal tests with Holm correction (all pairs, all metrics).
supp <- c(
  "==================================================================",
  " SUPPLEMENTARY: FORMAL TESTS (Holm-adjusted within each metric)",
  "==================================================================",
  "",
  " S1. One-way ANOVA across fit / ave / ext"
)
an <- anova_df[order(anova_df$metric, anova_df$param_set), ]
supp <- c(supp, sprintf("   %-16s %-20s F(%d, %d) = %-9s p %-9s p_holm %s",
                        an$metric, an$param_set, an$df_model, an$df_resid,
                        fmt(an$F), fmt_p(an$p), fmt_p(an$p_holm)))
supp <- c(supp, "", " S2. Pairwise Welch t-tests (B vs A)")
pw <- pairwise_df[order(pairwise_df$metric, pairwise_df$param_set, pairwise_df$B, pairwise_df$A), ]
supp <- c(supp, sprintf("   %-16s %-20s %s vs %s: t(%s) = %-9s p %-9s p_holm %s",
                        pw$metric, pw$param_set, pw$B, pw$A,
                        fmt(pw$welch_df, 3), fmt(pw$welch_t),
                        fmt_p(pw$welch_p), fmt_p(pw$welch_p_holm)))

cat(eff, supp, sep = "\n")
writeLines(c(eff, supp), out_path("sd_ratio_report.txt"))
