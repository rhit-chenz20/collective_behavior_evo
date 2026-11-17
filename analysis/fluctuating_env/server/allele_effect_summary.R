#!/usr/bin/env Rscript
# Which alleles carry the genetic variance under a fluctuating optimum vs under stabilizing
# selection? Per-sim summaries of the segregating alleles, averaged over snapshots.
#   fluctuating: fluc_env_adjusted_shift_size, all fluc_opt x fluc_interval, snapshots
#                fluc_snap (every 100 ticks from 10005)
#   stabilizing: cons_env_data_with_adjusted_shift_size, before the shift (ticks cons_snap,
#                optimum still 0); the three shift_size folders are pooled (identical treatment
#                before the shift)
#   both: n = 8, psi_set, reps 1-10
# Per snapshot, over segregating alleles (0 < p < 1), effect a, frequency p:
#   n_seg       number of segregating alleles
#   Va          genic variance sum 2 a^2 p (1 - p)
#   big_share   share of Va from alleles with |a| > big_a  (|a| > 1.5: ~10% of new mutations,
#               since a^2 ~ Exp(1))
#   va_abs_a    Va-weighted mean |a|: typical effect size of the alleles that carry the variance
#   big_maf     mean minor-allele frequency of segregating alleles with |a| > big_a
#   small_maf   same for |a| <= big_a
#
# Output: figure/r_data/fluc_env_adjusted_shift_size/allele_effects_n8_rep1_10.rds
#   (columns env = fluctuating / stabilizing; fluc_opt, fluc_interval NA for stabilizing)
# Run from the repo root (~10 min):
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/allele_effect_summary.R

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

# ---- settings ----
fluc_root <- "/fs/cbsubscb10/storage/zc524/fluc_env_adjusted_shift_size/parquet/allele"
cons_root <- "/fs/cbsubscb10/storage/zc524/cons_env_data_with_adjusted_shift_size/parquet/allele"
out_file  <- "figure/r_data/fluc_env_adjusted_shift_size/allele_effects_n8_rep1_10.rds"
psi_set   <- c(0, 0.2, 0.5, 0.8)
reps      <- 1:10
fluc_snap <- seq(10005, 10905, by = 100)
cons_snap <- c(9950, 9960, 9970, 9980, 9990)
big_a     <- 1.5

summ <- function(a) {
  a <- a[frequencies > 0 & frequencies < 1]
  s <- a[, {
    h <- 2 * effect_size^2 * frequencies * (1 - frequencies)
    big <- abs(effect_size) > big_a
    maf <- pmin(frequencies, 1 - frequencies)
    .(n_seg = .N, Va = sum(h), big_share = sum(h[big]) / sum(h),
      va_abs_a = sum(abs(effect_size) * h) / sum(h),
      big_maf = mean(maf[big]), small_maf = mean(maf[!big]))
  }, by = .(psi, rep, ID, gen)]
  s[, lapply(.SD, mean), by = .(psi, rep, ID), .SDcols = !"gen"]
}
read_sel <- function(path, gens, psi_val = NULL) {
  d <- open_dataset(path)
  g <- gens
  # OR-of-== on gen (see the arrow %in% quirk in data_storage_description.txt)
  filt <- Reduce(function(x, y) call("|", x, y), lapply(g, function(v) call("==", quote(gen), v)))
  max_rep <- max(reps)
  q <- d |> dplyr::filter(!!filt, rep <= !!max_rep)
  if (!is.null(psi_val)) q <- q |> dplyr::filter(psi == !!psi_val)
  as.data.table(q |> dplyr::select(gen, effect_size, frequencies, psi, rep, ID) |> dplyr::collect())
}

res <- list()
t_start <- Sys.time()
# stabilizing baseline
for (sz in c("9.4", "18.8", "28.2")) for (r in c("ave", "fit", "ext")) for (ps in psi_set) {
  a <- read_sel(file.path(cons_root, sprintf("shift_size=%s/reg=%s/n=8", sz, r)), cons_snap, ps)
  s <- summ(a)
  s[, `:=`(env = "stabilizing", reg = r, fluc_opt = NA_real_, fluc_interval = NA_integer_)]
  res[[length(res) + 1]] <- s
  cat(sprintf("stabilizing shift folder %s %s psi %s: %d sims, %.1f min\n", sz, r, ps, nrow(s),
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
# fluctuating
dirs <- CJ(fluc_opt = c("9.4", "28.2", "47.0", "65.8", "84.6", "103.4"),
           fluc_interval = c(20L, 50L, 80L, 110L), reg = c("ave", "fit", "ext"),
           psi = c("0.0", "0.2", "0.5", "0.8"))
for (i in seq_len(nrow(dirs))) {
  d <- dirs[i]
  path <- file.path(fluc_root, sprintf("fluc_opt=%s/fluc_interval=%d/reg=%s/n=8/part-psi_%s-reps_1-0.parquet",
                                       d$fluc_opt, d$fluc_interval, d$reg, d$psi))
  s <- summ(read_sel(path, fluc_snap))
  s[, `:=`(env = "fluctuating", reg = d$reg, fluc_opt = as.numeric(d$fluc_opt), fluc_interval = d$fluc_interval)]
  res[[length(res) + 1]] <- s
  if (i %% 12 == 0) cat(sprintf("[%d/%d] fluctuating, %.1f min\n", i, nrow(dirs),
                                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
res <- rbindlist(res, use.names = TRUE)
saveRDS(res, out_file)
cat("saved", out_file, ":", nrow(res), "sims\n")
print(res[, .(sims = .N), by = .(env, reg, psi)][order(env, reg, psi)])
