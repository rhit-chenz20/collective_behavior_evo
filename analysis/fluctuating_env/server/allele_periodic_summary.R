#!/usr/bin/env Rscript
# Per-sim summary of how allele frequencies follow the fluctuating optimum
# (fluc_env_adjusted_shift_size, n = 8, all fluc_opt x fluc_interval, psi_set, reps 1-10).
# Reads the allele parquet files directly, one file (= one psi, 10 reps) at a time, and keeps only
# per-sim numbers.
#
# Window: allele records every 10 ticks from 10005; the first N records with N * 10 a whole number
# of periods (period = 2 * fluc_interval: N = 100, 100, 96, 88 for periods 40, 100, 160, 220).
# Over a whole number of cycles the Fourier terms below are exact.
#
# For each allele (frequency p(t), 0 when absent, 1 when fixed), omega = 2 pi / period:
#   c = (2 / N) * sum_t p(t) exp(-i omega t)   (fundamental, complex; |c| = amplitude of the
#                                               periodic swing of p, in frequency units)
#   v = temporal variance of p over the window
# Per sim:
#   swing      mean |c| over alleles whose frequency reaches swing_min_freq in the window and that
#              are not fixed throughout (the alleles drawn in allele_traj.pdf)
#   coherence  sum(|c|^2 / 2) / sum(v) over all alleles: share of allele-frequency variance over
#              time that is at the forcing period (1 = pure periodic swings, ~0 = drift only)
#   gain       |C_g| / |C_opt|, with C_g = sum 2 * a * c (fundamental of the genetic mean
#              g = sum 2 a p) and C_opt the fundamental of the optimum: amplitude of the genetic
#              mean's periodic swing relative to the optimum's (1 = genetic mean swings as much)
#   lag        delay of the genetic mean's swing behind the optimum's, in generations (0..period)
#   gain_inphase  Re(C_g / C_opt): gain counting only the part in phase with the optimum
#
# Output: figure/r_data/fluc_env_adjusted_shift_size/allele_periodic_n8_rep1_10.rds
#
# Run from the repo root (reads 288 files; ~15-20 min on one core):
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/allele_periodic_summary.R

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

# ---- settings ----
root           <- "/fs/cbsubscb10/storage/zc524/fluc_env_adjusted_shift_size/parquet/allele"
out_file       <- "figure/r_data/fluc_env_adjusted_shift_size/allele_periodic_n8_rep1_10.rds"
n_val          <- 8
psi_set        <- c("0.0", "0.2", "0.5", "0.8")    # as in the file names
rep_block      <- "reps_1"                         # file part holding reps 1-10
first_rec      <- 10005
dt             <- 10
max_rec        <- 100                              # records 10005 ... 10995
swing_min_freq <- 0.05

# optimum at tick t: flips at 50 + k * interval, the first (tick 50) to +fluc_opt / 2
opt_at <- function(t, fluc_opt, interval) {
  k <- floor((t - 50) / interval) + 1
  ifelse(k %% 2 == 1, fluc_opt / 2, -fluc_opt / 2)
}

summarise_file <- function(path, fluc_opt, interval) {
  period <- 2 * interval
  N <- max(which((seq_len(max_rec) * dt) %% period == 0))
  times <- first_rec + (seq_len(N) - 1) * dt
  om <- 2 * pi / period
  last_rec <- max(times)
  a <- as.data.table(open_dataset(path) |>
    dplyr::filter(gen >= first_rec, gen <= last_rec) |>
    dplyr::select(gen, mut_id, effect_size, frequencies, psi, rep, ID) |>
    dplyr::collect())
  stopifnot(all(a$gen %in% times))
  a[, `:=`(cs = cos(om * gen), sn = sin(om * gen))]
  al <- a[, .(a = effect_size[1],
              re = 2 / N * sum(frequencies * cs), im = -2 / N * sum(frequencies * sn),
              m1 = sum(frequencies) / N, m2 = sum(frequencies^2) / N,
              maxp = max(frequencies), n_fixed = sum(frequencies >= 1)),
          by = .(psi, rep, ID, mut_id)]
  al[, `:=`(amp = sqrt(re^2 + im^2), v = pmax(m2 - m1^2, 0))]
  o <- opt_at(times, fluc_opt, interval)
  C_opt <- complex(real = 2 / N * sum(o * cos(om * times)), imaginary = -2 / N * sum(o * sin(om * times)))
  s <- al[, {
    Cg <- complex(real = sum(2 * a * re), imaginary = sum(2 * a * im))
    r <- Cg / C_opt
    .(swing = mean(amp[maxp >= swing_min_freq & n_fixed < N]),
      n_swing = sum(maxp >= swing_min_freq & n_fixed < N),
      coherence = sum(amp^2 / 2) / sum(v),
      gain = Mod(r), gain_inphase = Re(r),
      lag = ((-Arg(r)) %% (2 * pi)) / om)
  }, by = .(psi, rep, ID)]
  s[, `:=`(fluc_opt = fluc_opt, fluc_interval = interval, period = period, n_rec = N)]
  s
}

dirs <- CJ(fluc_opt = c("9.4", "28.2", "47.0", "65.8", "84.6", "103.4"),
           fluc_interval = c(20L, 50L, 80L, 110L), reg = c("ave", "fit", "ext"), psi = psi_set)
res <- vector("list", nrow(dirs))
t_start <- Sys.time()
for (i in seq_len(nrow(dirs))) {
  d <- dirs[i]
  path <- file.path(root, sprintf("fluc_opt=%s/fluc_interval=%d/reg=%s/n=%d/part-psi_%s-%s-0.parquet",
                                  d$fluc_opt, d$fluc_interval, d$reg, n_val, d$psi, rep_block))
  if (!file.exists(path)) { cat("missing:", path, "\n"); next }
  s <- summarise_file(path, as.numeric(d$fluc_opt), d$fluc_interval)
  s[, reg := d$reg]
  res[[i]] <- s
  cat(sprintf("[%3d/%d] opt %s int %d %s psi %s: %d sims, %.1f min\n", i, nrow(dirs), d$fluc_opt,
              d$fluc_interval, d$reg, d$psi, nrow(s), as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
res <- rbindlist(res)
setcolorder(res, c("fluc_opt", "fluc_interval", "period", "reg", "psi", "rep", "ID"))
saveRDS(res, out_file)
cat("saved", out_file, ":", nrow(res), "sims\n")
