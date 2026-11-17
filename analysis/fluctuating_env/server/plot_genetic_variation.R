#!/usr/bin/env Rscript
# Is genetic variation maintained or eroded under a fluctuating optimum, and does it sit in
# large-effect alleles? Each measure is shown relative to stabilizing selection with the same
# copying rule and psi (constant environment before its shift, ticks 9950-9999; n = 8).
# Tile = log2(fluctuating / stabilizing): 0 = same as under stabilizing selection, +1 = twice as
# much, -1 = half.
#   from the data table (mean over ticks 10000-10999, 20 reps; baseline: 60 reps x 3 shift folders):
#     g_sd          genetic sd of the trait (includes linkage disequilibrium)
#     H             heterozygosity per genome: SLiM's per-site H x maxL. The two datasets use
#                   different maxL (fluc 2500, cons 10000; same total mutation rate U and free
#                   recombination), so per-site H is not comparable between them.
#     seg_site_pop  number of segregating sites
#   from allele_effect_summary.R (10 reps; snapshots):
#     Va            genic variance sum 2 a^2 p (1-p)
#     big_share     share of Va carried by alleles with |a| > 1.5
#     va_abs_a      Va-weighted mean |a| (typical effect size of the alleles carrying the variance)
#     big_maf       mean minor-allele frequency of segregating large-effect alleles
# Also writes the stabilizing baselines (absolute values) to genetic_variation_baseline_n8.csv.
#
# Inputs:
#   figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_track_gen9950-11000.rds
#   figure/r_data/cons_env_data_with_adjusted_shift_size/data_summary_n8_baseline_gen9950-9999.rds
#     (load_data.R -d cons_env_data_with_adjusted_shift_size -t data --n 8 --psi 0,0.2,0.5,0.8
#      --cols z_mean,g_mean,z_sd,g_sd,H,seg_site_pop --gens 9950:9999 --tag n8_baseline_gen9950-9999)
#   figure/r_data/fluc_env_adjusted_shift_size/allele_effects_n8_rep1_10.rds
# Output: analysis/fluctuating_env/server/plots/genetic_variation/variation_<measure>_n8_8.7cm.pdf
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/plot_genetic_variation.R

source("analysis/fluctuating_env/server/fluc_heatmap_helpers.R")

fluc_file  <- "figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_track_gen9950-11000.rds"
base_file  <- "figure/r_data/cons_env_data_with_adjusted_shift_size/data_summary_n8_baseline_gen9950-9999.rds"
eff_file   <- "figure/r_data/fluc_env_adjusted_shift_size/allele_effects_n8_rep1_10.rds"
out_dir    <- "analysis/fluctuating_env/server/plots/genetic_variation"
maxL_fluc  <- 2500         # fluc_*.slim
maxL_cons  <- 10000        # cons_*.slim
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- data-table measures ----
dm <- c("g_sd", "H", "seg_site_pop")
fl <- as.data.table(readRDS(fluc_file))[tick >= 10000 & tick <= 10999]
fl[, H := H * maxL_fluc]
fl <- fl[, lapply(.SD, mean), by = .(reg, fluc_opt, fluc_interval, psi, ID), .SDcols = dm][
  , lapply(.SD, mean), by = .(reg, fluc_opt, fluc_interval, psi), .SDcols = dm]
bs <- as.data.table(readRDS(base_file))
bs[, H := H * maxL_cons]
bs <- bs[, lapply(.SD, mean), by = .(reg, psi, ID), .SDcols = dm][, lapply(.SD, mean), by = .(reg, psi), .SDcols = dm]

# ---- allele measures ----
am <- c("Va", "big_share", "va_abs_a", "big_maf")
ef <- as.data.table(readRDS(eff_file))
efl <- ef[env == "fluctuating", lapply(.SD, mean), by = .(reg, fluc_opt, fluc_interval, psi), .SDcols = am]
ebs <- ef[env == "stabilizing", lapply(.SD, mean), by = .(reg, psi), .SDcols = am]

base <- bs[ebs, on = .(reg, psi)]
fwrite(base[order(reg, psi)], file.path(out_dir, "genetic_variation_baseline_n8.csv"))
cat("stabilizing baseline:\n"); print(base[order(reg, psi)], digits = 3)

all_f <- fl[efl, on = .(reg, fluc_opt, fluc_interval, psi)]
titles <- c(g_sd = "Genetic sd", H = "Heterozygosity per genome", seg_site_pop = "Number of segregating sites",
            Va = "Genic variance", big_share = "Share of genic variance in large-effect alleles (|a| > 1.5)",
            va_abs_a = "Typical effect size carrying the variance (Va-weighted mean |a|)",
            big_maf = "Mean minor-allele frequency of large-effect alleles")
rat <- list()
for (k in names(titles)) {
  d <- all_f[, .(reg, fluc_opt, fluc_interval, psi, f = get(k))][base[, .(reg, psi, b = get(k))], on = .(reg, psi)]
  d[, val := log2(f / b)]
  rat[[k]] <- d[, .(measure = k, reg, fluc_opt, fluc_interval, psi, fluctuating = f, stabilizing = b, log2_ratio = val)]
  p <- hm_plot(d, sprintf("%s: fluctuating vs stabilizing selection", titles[[k]]), fill = "div",
               legend_title = "log2(fluctuating / stabilizing, same rule and psi)")
  f <- file.path(out_dir, sprintf("variation_%s_n8_8.7cm.pdf", k))
  ggsave(f, p, width = fig_w, height = fig_h)
  cat(sprintf("  %-12s log2 ratio %+.2f to %+.2f  %s\n", k, min(d$val), max(d$val), f))
}
fwrite(rbindlist(rat), file.path(out_dir, "genetic_variation_n8.csv"))
