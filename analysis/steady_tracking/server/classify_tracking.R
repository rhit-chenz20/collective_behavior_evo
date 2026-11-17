#!/usr/bin/env Rscript
# Classify every steady_tracking sim as tracking (success) or too unfit (failed) from its mean
# fitness at the end, instead of from single ticks (the metadata's `failed` is set by any one tick
# below 10%, so momentary dips mark sims that end up tracking as failed).
#
#   fit_500   mean over the last 500 ticks of (population mean fitness / fitness at the optimum)
#
#   tracking    fit_500 >= cutoff (10%)                                -> success
#   too_unfit   fit_500 <  cutoff                                      -> failed
#   aborted     SLiM stopped the sim (fitness underflowed, lag ~1700)  -> failed
#
# Also reported, not used in the rule: lag_500 (mean lag, optimum - mean phenotype, over the same
# ticks) and lag_slope_500 (its least-squares slope per tick; ~0 = steady lag).
#
# Input: figure/r_data/steady_tracking/{data_summary_all,data_summary_fast_speeds,crashed_sims}.*
#        figure/r_data/steady_tracking_long/{data_summary_all,metadata_all,crashed_sims}.*
# Output: analysis/steady_tracking/server/plots_fit500/tracking_class/
#   tracking_class_per_sim.csv       one row per sim (dataset = steady_tracking or long)
#   tracking_class_summary.csv       counts per combination
#   tracking_class_long_vs_12000.csv long runs: class at 15000 vs class at 12000 of the same sims
#
# Run from the repo root (~1 min):
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/classify_tracking.R

suppressPackageStartupMessages(library(data.table))
source("analysis/steady_tracking/server/steady_data.R")

out_dir <- "analysis/steady_tracking/server/plots_fit500/tracking_class"
cutoff  <- 0.1                        # fitness ratio cutoff
win     <- 500                        # ticks averaged at the end
peak    <- dnorm(0, 0, sqrt(2000))    # fitness at the optimum, Vs = 2N = 2000
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

classify_sims <- function(d, end_tick) {
  d <- d[tick > end_tick - win & tick <= end_tick]
  res <- d[, .(fit_500 = mean(fitness_mean) / peak,
               lag_500 = mean(z_opt),
               lag_slope_500 = coef(lm.fit(cbind(1, tick), z_opt))[2],
               n_ticks = .N),
           by = .(shift_speed, reg, psi, rep, ID)]
  res[, class := fifelse(fit_500 >= cutoff, "tracking", "too_unfit")]
  res[, end_tick := end_tick]
  stopifnot(all(res$n_ticks == win))   # every sim in the data reached end_tick
  res
}

read_cols <- c("shift_speed", "reg", "psi", "rep", "ID", "tick", "z_opt", "fitness_mean")
load_rds <- function(f) as.data.table(readRDS(f))[tick > 10000, ..read_cols]

# steady_tracking dataset (original + fast_speeds), end 12000
st <- rbindlist(lapply(st_files, function(f) {
  x <- classify_sims(load_rds(file.path(st_data_dir, f[1])), 12000); gc(verbose = FALSE); x
}))
st[, dataset := "steady_tracking"]
cr <- read_crashed(speeds = NULL)[, .(shift_speed, reg, psi, rep, ID, crash_tick, class = "aborted",
                         dataset = "steady_tracking")]

# long_borderline (same combinations, to 15000): classify at 15000 and at 12000
long_dir <- "figure/r_data/steady_tracking_long"
dl <- load_rds(file.path(long_dir, "data_summary_all.rds"))
lg   <- classify_sims(dl, 15000)[, dataset := "long"]
lg12 <- classify_sims(dl, 12000)
cr_l <- read_crashed(long_dir, speeds = NULL)
if (nrow(cr_l)) cr <- rbind(cr, cr_l[, .(shift_speed, reg, psi, rep, ID, crash_tick,
                                         class = "aborted", dataset = "long")])

ps <- rbindlist(list(st, lg, cr), use.names = TRUE, fill = TRUE)
ps[, outcome := fifelse(class == "tracking", "success", "failed")]

# the metadata's per-tick flag, for comparison
meta <- rbind(read_meta(keep_partial = TRUE, speeds = NULL)[, dataset := "steady_tracking"],
              as.data.table(readRDS(file.path(long_dir, "metadata_all.rds")))[, dataset := "long"],
              fill = TRUE)[, .(dataset, ID, meta_failed = failed == "T", n_fail_ticks)]
ps <- merge(ps, meta, by = c("dataset", "ID"), all.x = TRUE)
setcolorder(ps, c("dataset", "shift_speed", "reg", "psi", "rep", "ID", "outcome", "class", "fit_500"))
fwrite(ps[order(dataset, shift_speed, reg, psi, rep)], file.path(out_dir, "tracking_class_per_sim.csv"))

sm <- ps[, .(n = .N, success = sum(outcome == "success"), failed = sum(outcome == "failed"),
             too_unfit = sum(class == "too_unfit"), aborted = sum(class == "aborted"),
             meta_failed = sum(meta_failed, na.rm = TRUE),
             fit_500_med = signif(median(fit_500, na.rm = TRUE), 3),
             lag_500_med = round(median(lag_500, na.rm = TRUE), 1)),
         by = .(dataset, shift_speed, reg, psi)]
fwrite(sm[order(dataset, reg, psi, shift_speed)], file.path(out_dir, "tracking_class_summary.csv"))

cmp <- merge(lg12[, .(ID, class_12000 = class, fit_500_12000 = fit_500)],
             lg[, .(ID, shift_speed, reg, psi, rep, class_15000 = class, fit_500_15000 = fit_500)], by = "ID")
fwrite(cmp[order(shift_speed, reg, psi, rep)], file.path(out_dir, "tracking_class_long_vs_12000.csv"))

cat("per sim:\n");            print(ps[, .N, by = .(dataset, outcome, class)][order(dataset, outcome, class)])
cat("\nvs metadata failed:\n"); print(ps[dataset == "steady_tracking", .N, by = .(outcome, meta_failed)][order(outcome, meta_failed)])
cat("\nlong runs, class at 12000 vs 15000:\n"); print(cmp[, .N, by = .(class_12000, class_15000)])
cat("\nwrote", out_dir, "\n")
