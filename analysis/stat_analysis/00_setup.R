# Parameters, packages and shared helper functions for the stat_analysis scripts
# (split from analysis/stat_analysis.Rmd). Every section script sources this file, so edit
# the parameters here (trait, date, n_boot, ...). Paths are relative to the project root.

# parameters for plotting
date <- "10_27_2026"

# set these once for your machine / project
analysis_root <- "analysis"
base_plot_folder     <- file.path("analysis/figure", date)

# ---- data files ----
# load_data.R (analysis/load_data_from_storage) saves each load into its own folder,
#   analysis/rds/<dataset>/<name>/<name>.rds   (<name> = <type>_<tag>; type data is data_summary)
# next to <name>_summary.txt, _params.csv, _sims.csv and _columns.csv describing it.
# The two environments are different datasets, so each has its own <dataset> folder.
rds_root        <- file.path("analysis", "rds")
cons_dataset    <- "cons_env_data_with_adjusted_shift_size"   # constant environment
fluc_dataset    <- "fluc_env_constant_psi"                    # fluctuating environment
cons_rds_folder <- file.path(rds_root, cons_dataset)
fluc_rds_folder <- file.path(rds_root, fluc_dataset)

# path of a loaded .rds by name; files saved before the per-load folders (directly in
# <dataset>/) are still found
rds_file <- function(folder, name) {
  new <- file.path(folder, name, paste0(name, ".rds"))
  old <- file.path(folder, paste0(name, ".rds"))
  if (!file.exists(new) && file.exists(old)) old else new
}

# which file each script reads (names as given by --tag in load_data.R)
burnin_file  <- rds_file(cons_rds_folder, "data_summary_all_nosfs")  # sections 1-9: cons_data_df
px_cons_file <- rds_file(cons_rds_folder, "data_summary_fit_zg_all") # section 10: px_data_df
# section 10 fluc: fluc_data_df; needs fitness_mean at every tick (no --gen-values thinning)
px_fluc_file <- rds_file(fluc_rds_folder, "data_summary_fitness_all")

# readRDS with a clear error naming the files that are in that folder
read_rds_data <- function(file) {
  if (!file.exists(file)) {
    d <- dirname(file)
    if (!dir.exists(d)) d <- dirname(d)   # per-load folder missing: list the dataset folder
    have <- basename(list.files(d, pattern = "\\.rds$", recursive = TRUE))
    stop("Data file not found: ", file, "\n  .rds files in ", d, ": ",
         if (length(have)) paste(have, collapse = ", ") else "(none, or folder missing)",
         "\n  Set the file in 00_setup.R, or create it with load_data.R",
         " (see analysis/load_data_from_storage/load_examples.sh).", call. = FALSE)
  }
  message("Reading ", file)
  readRDS(file)
}

library(ggplot2)
library(MetBrewer)

if (!dir.exists(base_plot_folder)) {
  dir.create(base_plot_folder, recursive = TRUE, showWarnings = FALSE)
}

met_brewer_pallete <- "Hokusai2"  # color pattern for lines

burnin_end <- 10000               # only ticks < burnin_end are used (shift_size irrelevant here)
n_boot     <- 10000               # bootstrap resamples
models     <- c("ave", "fit", "ext")
model_cols <- setNames(met.brewer(name = met_brewer_pallete, n = 3, type = "discrete"), models)
set.seed(20260928)                # makes bootstrap CIs reproducible

# which trait to analyse: "z" = phenotype (z_sd, z_mean), "g" = genotype (g_sd, g_mean).
# Set this and re-run the whole file; outputs are prefixed with the trait so the
# phenotype and genotype results do not overwrite each other.
trait      <- "g"
sd_col     <- paste0(trait, "_sd")
mean_col   <- paste0(trait, "_mean")
trait_word <- c(z = "phenotypic", g = "genotypic")[[trait]]
Trait_word <- c(z = "Phenotypic", g = "Genotypic")[[trait]]
out_path   <- function(name) file.path(base_plot_folder, paste0(trait, "_", name))

# section 10 (psi = x vs psi = 0) output folder, in subfolders:
#   constant_only/  10_psi_vs_psi0.R and 10_psi_vs_psi0_median.R: table csv + pdf
#   results/        10_psi_vs_psi0_fluc.R: the full comparison csv (every setting, metric
#                   and n), its settings csv and the peak / worst phase csv
#   overview/       10_psi_vs_psi0_fluc_plot.R: combined table and per-metric fluctuating plots
#   paper/          10_psi_vs_psi0_figure.R and 10_psi_vs_psi0_figure_n8.R
px_folder <- file.path(base_plot_folder, "psi_vs_psi0")
px_out <- function(sub, name) {
  d <- file.path(px_folder, sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, name)
}

# ---- shared helpers (defined in sections 3-7 of the Rmd, used by later sections) ----
# Bootstrap mean(B) - mean(A) and mean(B) / mean(A) by resampling runs within each model.
boot_diff_ratio <- function(a, b, n_boot) {
  # index-based resampling (sample(x) misbehaves when length(x) == 1)
  ma <- colMeans(matrix(a[sample.int(length(a), length(a) * n_boot, replace = TRUE)], nrow = length(a)))
  mb <- colMeans(matrix(b[sample.int(length(b), length(b) * n_boot, replace = TRUE)], nrow = length(b)))
  c(diff_lo  = unname(quantile(mb - ma, 0.025)), diff_hi  = unname(quantile(mb - ma, 0.975)),
    ratio_lo = unname(quantile(mb / ma, 0.025)), ratio_hi = unname(quantile(mb / ma, 0.975)),
    # 90% CI: used for the equivalence classification (TOST at alpha = 0.05)
    ratio_lo90 = unname(quantile(mb / ma, 0.05)), ratio_hi90 = unname(quantile(mb / ma, 0.95)))
}

fmt   <- function(x, d = 4) ifelse(is.na(x), "NA", formatC(signif(x, d), format = "g", digits = d))
# returns "= 0.0123" or "< 1e-4" so it reads correctly after "p"
fmt_p <- function(p) ifelse(is.na(p), "= NA", ifelse(p < 1e-4, "< 1e-4", paste("=", formatC(p, format = "f", digits = 4))))
metric_labels <- c(
  mean_sd_all      = paste(Trait_word, "SD, averaged over all 50 time points"),
  mean_mean_all    = paste(Trait_word, "mean, averaged over all 50 time points")
)

heat_label_col <- "grey15"

# The narrative is built on the SD ratio (B/A) and its bootstrap CI, not on p-values.
# Equivalence margin: an SD ratio within +/-5% is treated as biologically negligible.
# Classification uses the 90% CI, which is the formal TOST at alpha = 0.05
# (equivalent = both one-sided tests reject). The 95% CI is reported alongside.
equiv_margin <- c(0.95, 1.05)

classify_equiv <- function(lo, hi, margin) {
  ifelse(is.na(lo) | is.na(hi), NA,
  ifelse(lo >= margin[1] & hi <= margin[2], "equivalent",   # CI fully inside margin
  ifelse(hi <  margin[1] | lo >  margin[2], "different",    # CI fully outside margin
                                            "inconclusive")))
}

class_code <- c(equivalent = "E", different = "D", inconclusive = "?")

# Generic psi x n ratio heatmap (d needs columns psi, n, ratio, class); also used in section 6.
ratio_heatmap <- function(d, title, fill_name = "SD ratio") {
  d$label <- sprintf("%.2f\n%s", d$ratio, class_code[d$class])
  # color limits fitted to this plot's range (log scale), always spanning the
  # equivalence margin; 1 stays the neutral midpoint
  lims <- c(min(log2(d$ratio), log2(equiv_margin[1]), na.rm = TRUE),
            max(log2(d$ratio), log2(equiv_margin[2]), na.rm = TRUE))
  lim  <- max(abs(lims))
  # round-number ticks: coarse for large ranges (e.g. ext/fit ~20x), fine for small ones
  brk <- if (lim > log2(6)) c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50, 100) else
         if (lim > log2(1.3)) c(0.25, 0.5, 0.67, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 2, 3, 4) else
           c(0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.25)
  brk <- brk[log2(brk) >= lims[1] & log2(brk) <= lims[2]]
  ggplot(d, aes(x = factor(n), y = factor(psi), fill = log2(ratio))) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.6, lineheight = 0.85, color = heat_label_col) +
    scale_fill_gradient2(low = "#4575b4", mid = "#f0efec", high = "#d73027", midpoint = 0,
                         limits = lims, breaks = log2(brk),
                         labels = brk, name = fill_name) +
    labs(x = "n", y = expression(psi), title = title,
         caption = sprintf(paste("E = equivalent, D = different, ? = inconclusive;",
                                 "90%% bootstrap CI vs margin [%.2f, %.2f]"),
                           equiv_margin[1], equiv_margin[2])) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank())
}

f3 <- function(x) sprintf("%.3f", x)
ci3 <- function(r, lo, hi) sprintf("%.3f [%.3f, %.3f]", r, lo, hi)
rule <- "------------------------------------------------------------------"

optimum_burnin <- 0       # optimum before the shift (SLiM: optimum = 0.0 until shift_time)

# Bootstrap a cell mean over runs: returns 90% and 95% percentile CIs.
boot_mean_ci <- function(x, n_boot) {
  bm <- colMeans(matrix(x[sample.int(length(x), length(x) * n_boot, replace = TRUE)], nrow = length(x)))
  c(lo90 = quantile(bm, 0.05, names = FALSE), hi90 = quantile(bm, 0.95, names = FALSE),
    lo   = quantile(bm, 0.025, names = FALSE), hi  = quantile(bm, 0.975, names = FALSE))
}

agree_one <- function(obs, pred) {
  ok <- obs > 0 & pred > 0 & is.finite(obs / pred)
  r  <- obs[ok] / pred[ok]; ci <- boot_mean_ci(r, n_boot)
  f  <- lm(log(obs[ok]) ~ log(pred[ok]))
  data.frame(n_runs = sum(ok), obs_over_pred = mean(r), lo = ci[["lo"]], hi = ci[["hi"]],
             loglog_slope = coef(f)[[2]], r_squared = summary(f)$r.squared)
}

grid_txt <- function(d, value_fmt = "%.2f%s", with_class = TRUE) {
  lab <- if (with_class) sprintf(value_fmt, d$ratio, class_code[d$class]) else sprintf(value_fmt, d$ratio)
  g <- tapply(lab, list(d$psi, d$n), identity); g[is.na(g)] <- "-"
  paste0("     ", capture.output(print(noquote(g))))
}
