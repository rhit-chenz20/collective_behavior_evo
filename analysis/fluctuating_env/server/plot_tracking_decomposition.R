#!/usr/bin/env Rscript
# Who tracks the fluctuating optimum: genes or copying?
# Mean phenotype z = mean genotype g + social part s (s = z - g). By the model's phenotype rules
# (fluc_*.slim, calculate_phenotype):
#   ave: z_i = ((n-1)(1-psi) g_i + psi sum_group g) / (n-1+psi)  ->  mean z = mean g exactly (s = 0)
#   fit: z_i = (1-psi) g_i + psi g_best (group member closest to the optimum)  -> s pulled toward
#        the optimum
#   ext: z_i = (1-psi) g_i + psi g_extreme (member furthest from the group mean, either side)
#        -> s ~ 0 on average
# For each sim, over a whole number of periods from tick 10000 (window <= 1000 ticks), the
# fundamental (wave at the forcing period, see allele_periodic_summary.R) of z_mean, g_mean and
# the optimum. Per sim:
#   gain_z        |C_z| / |C_opt|: how far the mean phenotype swings relative to the optimum
#   gain_g        |C_g| / |C_opt|: same for the mean genotype
#   social_share  Re((C_z - C_g) / C_z): share of the phenotype's swing supplied by copying
#   lag_z         delay of the phenotype's swing behind the optimum, fraction of the period
# Tiles = mean over the 20 reps.
#
# Input:  figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_track_gen9950-11000.rds
#         (load_data.R -d fluc_env_adjusted_shift_size -t data --n 8 --psi 0,0.2,0.5,0.8
#          --cols optimum,z_mean,g_mean,z_sd,g_sd,H,seg_site_pop --gens 9950:11000
#          --tag n8_track_gen9950-11000)
# Output: analysis/fluctuating_env/server/plots/tracking/tracking_<measure>_n8_8.7cm.pdf, .csv
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/plot_tracking_decomposition.R

source("analysis/fluctuating_env/server/fluc_heatmap_helpers.R")

in_file <- "figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_track_gen9950-11000.rds"
out_dir <- "analysis/fluctuating_env/server/plots/tracking"
t_from  <- 10000
max_len <- 1000
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s <- as.data.table(readRDS(in_file))
s[, period := 2L * fluc_interval]
s[, len := whole_cycles(period, max_len)]
w <- s[tick >= t_from & tick < t_from + len]
per_sim <- w[, {
  co <- fourier1(optimum, tick, period[1]); cz <- fourier1(z_mean, tick, period[1])
  cg <- fourier1(g_mean, tick, period[1])
  .(gain_z = Mod(cz) / Mod(co), gain_g = Mod(cg) / Mod(co),
    social_share = Re((cz - cg) / cz),
    lag_z = ((Arg(co) - Arg(cz)) %% (2 * pi)) / (2 * pi), n_tick = .N)
}, by = .(reg, fluc_opt, fluc_interval, psi, rep, ID)]
stopifnot(all(per_sim$n_tick == whole_cycles(2L * per_sim$fluc_interval, max_len)))
m <- per_sim[, c(lapply(.SD, mean), list(reps = .N)), by = .(reg, fluc_opt, fluc_interval, psi),
             .SDcols = c("gain_z", "gain_g", "social_share", "lag_z")]
fwrite(m[order(reg, psi, fluc_interval, fluc_opt)], file.path(out_dir, "tracking_n8.csv"))
cat(sprintf("%d sims, %d-%d reps per tile\n", nrow(per_sim), min(m$reps), max(m$reps)))
print(m[, .(social_share_min = min(social_share), social_share_max = max(social_share)), by = reg])

plots <- list(
  gain_z = list("Phenotype tracking: swing of mean phenotype / swing of optimum", "seq", c(0, 1)),
  gain_g = list("Genetic tracking: swing of mean genotype / swing of optimum", "seq", c(0, 1)),
  social_share = list("Share of the phenotype's swing supplied by copying", "seq", c(0, 1)),
  lag_z = list("Lag of mean phenotype behind the optimum (fraction of period)", "seq", NULL))
for (k in names(plots)) {
  pl <- plots[[k]]
  d <- copy(m)[, val := get(k)]
  p <- hm_plot(d, pl[[1]], fill = pl[[2]], limits = pl[[3]])
  f <- file.path(out_dir, sprintf("tracking_%s_n8_8.7cm.pdf", k))
  ggsave(f, p, width = fig_w, height = fig_h)
  cat(sprintf("  %-13s %.3f - %.3f  %s\n", k, min(d$val), max(d$val), f))
}
