# psi = x vs psi = 0 across selection regimes, constant + fluctuating environment
# (section 10 of stat_analysis.Rmd, with fluctuating-environment panels added).
# Run from the project root:
#   source("analysis/stat_analysis/10_psi_vs_psi0_fluc.R")
# Parameters and data file paths are in 00_setup.R; data is loaded by 00_load_data_px.R.
#
# Data (per-tick summaries from load_data.R, saved in analysis/rds/<dataset>/):
#   constant env (px_data_df):     tick, psi, n, reg, rep, shift_size, ID, fitness_mean, z_mean, g_mean
#   fluctuating env (fluc_data_df): tick, psi, n, reg, rep, fluc_opt, fluc_interval, ID, fitness_mean
#                    (every tick: do not thin generations, the half-cycle min / max need them)
# (the parameter columns, rep, ID and tick are always kept by load_data.R)
#
# For each model and n, and each metric, compare the mean over runs at every psi = x > 0
# with the mean at psi = 0. A cell is "higher" if psi = x is above psi = 0 for every x,
# "lower" if below for every x, and "mixed" otherwise. Colors and labels are by "better"
# (blue) vs "worse" (red): higher is better for fitness, lower is better for the
# generations to reach the shift (px_lower_better); the
# label gives how many x are better.
# Metrics, one value per run:
#   stabilizing: fitness_mean averaged over the burn-in (ticks < shift_time; recorded from 9950)
#   directional: fitness_mean at shift_time (generation 10000);
#                for each fraction f in px_reach_frac, generations after the shift until
#                z_mean / g_mean first reach optimum_burnin + f * shift_size (assumes an
#                upward shift). Runs that have not reached it by px_reach_deadline count as
#                slower than every run that did (Inf), so the mean of a psi with any such
#                run is Inf: a psi = x with one is worse than psi = 0 (not better when
#                psi = 0 has one too), and every psi = x without one is better than a
#                psi = 0 with one. These psi are listed in psi_unreached in the csv.
#                phenotype lag time for each T in px_lag_T: sum over generations
#                10000 .. 10000 + T - 1 of (optimum_burnin + shift_size - z_mean) / shift_size
#                (lower = faster; compared by the mean over runs).
#   fluctuating: one panel per (fluc_opt, fluc_interval). The optimum flips at ticks
#                50 + k * fluc_interval; a half cycle is the fluc_interval ticks from one
#                flip up to the next (the tick of the flip already has the new optimum,
#                see px_flip_lag). Only half cycles entirely inside px_fluc_window are used.
#                cycle average: fitness_mean averaged over whole cycles (an even number of
#                  half cycles, so both optima count equally)
#                peak:  max fitness_mean within each half cycle, averaged over half cycles
#                       (expected at the last tick before the next flip)
#                worst: min fitness_mean within each half cycle, averaged over half cycles
#                       (expected at the flip tick)
#                peak / worst (edge): mean fitness_mean over the last / first px_edge_k ticks
#                       of each half cycle, averaged over half cycles (less sensitive to
#                       tick-to-tick noise than max / min)
#                Where the max / min actually fall is written to
#                psi_vs_psi0_fluc_peak_worst_phase.csv.
# The fluctuating dataset may be incomplete (rerun in progress): cells use the runs present
# (n_runs in the csv).

source("analysis/stat_analysis/00_setup.R")
px_load_fluc <- TRUE
source("analysis/stat_analysis/00_load_data_px.R")
library(dplyr)   # loaded here only: it masks stats::filter / lag, which sections 1-9 do not need

# ---- settings ----
shift_time     <- 10000
px_psi         <- NULL            # psi = x values compared with psi = 0; NULL = all in the data
px_psi_exclude <- c()            # psi values dropped before anything is computed; c() = none
px_reach_frac  <- c(0.5)     # fractions of shift_size for the "generations to reach" metrics
px_reach_deadline <- 11000        # every run must reach each fraction by this generation
px_lag_T       <- c(25, 100)      # windows (generations after the shift) for the phenotype lag time
px_fluc_window <- c(9950, 11000)  # ticks (inclusive) used for the fluctuating metrics
# ticks from a flip (50 + k * fluc_interval) to the first logged fitness under the new optimum.
# 0: the flip late() event is defined before the late() that computes abs_fit, so the flip
# tick is already logged with the new optimum. If the phase csv shows the min at phase 1 and
# the max at phase 0, set this to 1.
px_flip_lag    <- 0
# ticks averaged for the edge versions of peak / worst, which avoid picking up tick-to-tick
# noise the way max / min do: peak = mean of the last px_edge_k ticks of each half cycle,
# worst = mean of the first px_edge_k ticks (from the flip tick)
px_edge_k      <- 5
# which settings are plotted, and colors: see 10_psi_vs_psi0_fluc_plot.R


# keep the needed columns, round psi, drop excluded psi, add env and a run key
# (ID is unique only within a dataset, so the run key includes env)
px_prep <- function(df, env, cols) {
  names(df)[names(df) == "ticks"] <- "tick"
  df <- df[df$reg %in% models, intersect(cols, names(df))]
  df$psi <- round(df$psi, 2)  # psi may be stored as float32 (0.2 -> 0.2000000029)
  if (length(px_psi_exclude)) df <- df[!df$psi %in% round(px_psi_exclude, 2), ]
  df$env <- env
  df$run <- if ("ID" %in% names(df)) paste(env, df$ID) else
    paste(env, df$shift_size, df$fluc_opt, df$fluc_interval, df$reg, df$n, df$psi, df$rep)
  df
}
if (0 %in% px_psi_exclude) stop("psi = 0 is the reference and cannot be excluded")
info_cols <- c("tick", "psi", "n", "reg", "rep", "shift_size", "fluc_opt", "fluc_interval", "ID")
px_df  <- px_prep(px_data_df, "cons", c(info_cols, "fitness_mean", "z_mean", "g_mean"))
fl_df  <- px_prep(fluc_data_df, "fluc", c(info_cols, "fitness_mean"))
px_info <- do.call(rbind, lapply(list(px_df, fl_df), function(d) {
  d <- d[!duplicated(d$run), ]
  for (k in c("shift_size", "fluc_opt", "fluc_interval")) if (is.null(d[[k]])) d[[k]] <- NA
  d[c("run", "env", "shift_size", "fluc_opt", "fluc_interval", "reg", "n", "psi", "rep")]
}))

# per-run values (named by run) -> rows of px_run
px_rows <- function(values, selection, metric, by_shift = TRUE) {
  info <- px_info[match(names(values), px_info$run), ]
  data.frame(selection = selection, metric = metric,
             shift_size = if (by_shift) info$shift_size else NA,
             fluc_opt = info$fluc_opt, fluc_interval = info$fluc_interval,
             reg = info$reg, n = info$n, psi = info$psi, rep = info$rep,
             value = as.numeric(values))
}

# stabilizing: the burn-in is shared across shift sizes, so keep one shift_size per run
# (as in sections 2-9). fitness_mean is 0 (not recorded) before shift_time - 50.
st <- px_df[px_df$tick < shift_time & !is.na(px_df$fitness_mean) & px_df$fitness_mean > 0, ]
st_key <- interaction(st$psi, st$n, st$reg, st$rep, drop = TRUE)
st <- st[st$shift_size == ave(st$shift_size, st_key, FUN = function(x) min(x, na.rm = TRUE)), ]
px_list <- list(px_rows(tapply(st$fitness_mean, st$run, mean), "Stabilizing",
                        "Mean fitness\n(gen < 10000)", by_shift = FALSE))

# directional: fitness right at the shift
d0 <- px_df[px_df$tick == shift_time & !is.na(px_df$fitness_mean), ]
px_list[[2]] <- px_rows(tapply(d0$fitness_mean, d0$run, mean), "Directional",
                        "Fitness at\ngen 10000")

# directional: generations until the population mean reaches a fraction of the shift
post <- px_df[px_df$tick >= shift_time & px_df$tick <= px_reach_deadline, ]
runs <- as.character(unique(post$run))
for (f in px_reach_frac) for (col in c("z_mean", "g_mean")) {
  tr  <- c(z_mean = "phenotype", g_mean = "genotype")[[col]]
  hit <- !is.na(post[[col]]) & post[[col]] >= optimum_burnin + f * post$shift_size
  t_hit <- tapply(post$tick[hit], post$run[hit], min)
  t_reach <- setNames(rep(Inf, length(runs)), runs)       # Inf = not reached by the deadline
  t_reach[names(t_hit)] <- t_hit - shift_time
  px_list[[length(px_list) + 1]] <- px_rows(t_reach, "Directional",
    paste0("Gens to ", 100 * f, "% of shift\n(", tr, ")"))
}

# directional: phenotype lag time, the sum over the first T generations after the shift of
# the fraction of the shift still to go, (opt - z_mean) / shift_size with
# opt = optimum_burnin + shift_size. In generations; lower = faster. Unlike the crossing
# time it is not biased by generation-to-generation noise in z_mean. NA if a generation of
# the window is missing.
for (TT in px_lag_T) {
  w <- px_df[px_df$tick >= shift_time & px_df$tick < shift_time + TT & !is.na(px_df$z_mean), ]
  s <- rowsum(cbind(1, (optimum_burnin + w$shift_size - w$z_mean) / w$shift_size), w$run)
  lag <- setNames(ifelse(s[, 1] == TT, s[, 2], NA), rownames(s))
  px_list[[length(px_list) + 1]] <- px_rows(lag, "Directional",
    paste0("Phenotype lag behind optimum\n(first ", TT, " gens)"), by_shift = TRUE)
}

# fluctuating: per half cycle min / max / mean / edge means of fitness_mean, then averaged per run
if (px_edge_k < 1 || 2 * px_edge_k > min(fl_df$fluc_interval))
  stop("px_edge_k must be between 1 and half the shortest fluc_interval (", min(fl_df$fluc_interval), ")")
fw <- fl_df[fl_df$tick >= px_fluc_window[1] & fl_df$tick <= px_fluc_window[2] &
              !is.na(fl_df$fitness_mean) & fl_df$fitness_mean > 0,
            c("run", "fluc_opt", "fluc_interval", "tick", "fitness_mean")]
fw$half  <- (fw$tick - 50 - px_flip_lag) %/% fw$fluc_interval    # half-cycle index
fw$phase <- (fw$tick - 50 - px_flip_lag) %% fw$fluc_interval     # ticks since the flip
fw <- fw[order(fw$run, fw$tick), ]
hc <- fw |>
  group_by(run, fluc_opt, fluc_interval, half) |>
  summarise(n_ticks = n(), hc_mean = mean(fitness_mean),
            peak = max(fitness_mean), worst = min(fitness_mean),
            peak_edge  = mean(fitness_mean[phase >= fluc_interval[1] - px_edge_k]),
            worst_edge = mean(fitness_mean[phase < px_edge_k]),
            peak_phase = phase[which.max(fitness_mean)],
            worst_phase = phase[which.min(fitness_mean)], .groups = "drop") |>
  filter(n_ticks == fluc_interval)                               # whole half cycles only
if (!nrow(hc)) stop("No complete half cycle in the fluctuating data: it needs every tick",
                    " (load without --gen-values / --gens thinning).")
# keep an even number of half cycles per run so the cycle average covers whole cycles
hc <- hc |> group_by(run) |> filter(rank(half) <= 2 * (n() %/% 2)) |> ungroup()
fl_run <- hc |> group_by(run) |>
  summarise(cycle_avg = mean(hc_mean), peak = mean(peak), worst = mean(worst),
            peak_edge = mean(peak_edge), worst_edge = mean(worst_edge), n_half = n())
px_list[[length(px_list) + 1]] <- px_rows(setNames(fl_run$cycle_avg, fl_run$run),
  "Fluctuating", "Cycle average\nfitness")
px_list[[length(px_list) + 1]] <- px_rows(setNames(fl_run$peak, fl_run$run),
  "Fluctuating", "Peak fitness\n(half-cycle max)")
px_list[[length(px_list) + 1]] <- px_rows(setNames(fl_run$peak_edge, fl_run$run),
  "Fluctuating", paste0("Peak fitness\n(last ", px_edge_k, " ticks)"))
px_list[[length(px_list) + 1]] <- px_rows(setNames(fl_run$worst, fl_run$run),
  "Fluctuating", "Worst fitness\n(half-cycle min)")
px_list[[length(px_list) + 1]] <- px_rows(setNames(fl_run$worst_edge, fl_run$run),
  "Fluctuating", paste0("Worst fitness\n(first ", px_edge_k, " ticks)"))

# where in the half cycle the max / min fall (checks px_flip_lag and the expectation that
# the peak is the last tick and the worst the flip tick)
fl_phase <- hc |>
  group_by(fluc_opt, fluc_interval) |>
  summarise(n_half_cycles = n(),
            frac_peak_at_last  = mean(peak_phase == fluc_interval - 1),
            frac_worst_at_flip = mean(worst_phase == 0),
            median_peak_phase  = median(peak_phase),
            median_worst_phase = median(worst_phase), .groups = "drop") |>
  arrange(fluc_opt, fluc_interval) |> as.data.frame()
write.csv(fl_phase, px_out("results", "psi_vs_psi0_fluc_peak_worst_phase.csv"), row.names = FALSE)
print(fl_phase)

px_run <- do.call(rbind, px_list)
px_lower_better <- grep("^Gens to |^Phenotype lag", unique(px_run$metric), value = TRUE)  # fewer gens / less lag
rm(st, d0, post, fw, hc)

# panels: one for stabilizing, one per shift size for directional,
# one per (fluc_opt, fluc_interval) for fluctuating
px_run$panel <- ifelse(px_run$selection == "Fluctuating",
                       paste0("Fluctuating\nopt = ", px_run$fluc_opt, ", int = ", px_run$fluc_interval),
                       ifelse(is.na(px_run$shift_size), px_run$selection,
                              paste0(px_run$selection, "\nshift = ", px_run$shift_size)))
fl_set <- unique(px_info[px_info$env == "fluc", c("fluc_opt", "fluc_interval")])
fl_set <- fl_set[order(fl_set$fluc_opt, fl_set$fluc_interval), ]
px_panels <- c("Stabilizing",
               paste0("Directional\nshift = ", sort(unique(px_info$shift_size))),
               paste0("Fluctuating\nopt = ", fl_set$fluc_opt, ", int = ", fl_set$fluc_interval))

# compare each psi = x with psi = 0 within panel x metric x model x n
px_all_psi   <- sort(unique(px_run$psi))
px_cell_list <- list()
for (g in split(px_run, list(px_run$panel, px_run$metric, px_run$reg, px_run$n), drop = TRUE)) {
  # NA if any run of that psi has no value; Inf if any run did not reach the threshold
  m    <- tapply(g$value, g$psi, mean)
  psis <- as.numeric(names(m))
  unreached <- is.infinite(m)
  if (!any(psis == 0)) next
  m0 <- m[[which(psis == 0)]]
  tested <- psis != 0 & (if (is.null(px_psi)) TRUE else psis %in% px_psi)
  if (!any(tested)) next
  dropped <- tested & is.na(m)                # psi = x with unsatisfied runs: not compared
  xs <- if (is.na(m0)) tested & FALSE else tested & !is.na(m)
  above <- m[xs] > m0
  better <- if (g$metric[1] %in% px_lower_better) m[xs] < m0 else above
  classify <- function(v, yes, no) if (!length(v)) NA else if (all(v)) yes else if (!any(v)) no else "mixed"
  # mean over runs at psi = 0 and at each tested psi = x (NA if not tested or not all runs have a value)
  means <- setNames(rep(NA_real_, length(px_all_psi)), paste0("mean_psi_", px_all_psi))
  means[paste0("mean_psi_", psis[tested | psis == 0])] <- m[tested | psis == 0]
  px_cell_list[[length(px_cell_list) + 1]] <- data.frame(
    panel = g$panel[1], metric = g$metric[1], model = g$reg[1], n = g$n[1],
    fluc_opt = g$fluc_opt[1], fluc_interval = g$fluc_interval[1],
    n_x = sum(xs), n_higher = sum(above),
    class = classify(above, "higher", "lower"),
    psi_lower = paste(psis[xs][!above], collapse = ", "),   # x values with psi = 0 >= psi = x
    n_better = sum(better),
    color_class = classify(better, "better", "worse"),
    psi0_complete = !is.na(m0),
    psi_dropped = paste(psis[dropped], collapse = ", "),    # x values left out (missing values)
    psi_unreached = paste(psis[unreached], collapse = ", "), # psi with runs not reaching the threshold
    n_runs = nrow(g), n_runs_na = sum(is.na(g$value)), as.list(means), check.names = FALSE)
}
px_cell <- do.call(rbind, px_cell_list)
px_cell$panel  <- factor(px_cell$panel, levels = px_panels)
px_cell$metric <- factor(px_cell$metric, levels = unique(px_run$metric))
px_cell$model  <- factor(px_cell$model, levels = models)
px_cell$color_class <- factor(px_cell$color_class, levels = c("better", "worse", "mixed"))
px_csv <- px_cell[order(px_cell$panel, px_cell$metric, px_cell$model, px_cell$n), ]
px_csv$panel  <- gsub("\n", " ", px_csv$panel)    # one line per row in the csv
px_csv$metric <- gsub("\n", " ", px_csv$metric)
write.csv(px_csv, px_out("results", "psi_vs_psi0_table_with_fluc.csv"), row.names = FALSE)

# settings used, for the plot captions
write.csv(data.frame(
  setting = c("fluc_window", "reach_deadline", "reach_frac", "flip_lag", "edge_k", "psi_exclude", "psi"),
  value   = c(paste(px_fluc_window, collapse = "-"), px_reach_deadline,
              paste(px_reach_frac, collapse = ", "), px_flip_lag, px_edge_k,
              paste(px_psi_exclude, collapse = ", "),
              if (is.null(px_psi)) "all" else paste(px_psi, collapse = ", "))),
  px_out("results", "psi_vs_psi0_table_with_fluc_settings.csv"), row.names = FALSE)

# plots, drawn from the csv (run 10_psi_vs_psi0_fluc_plot.R alone to redo only the plots)
source("analysis/stat_analysis/10_psi_vs_psi0_fluc_plot.R")
