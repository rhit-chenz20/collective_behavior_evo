#!/usr/bin/env Rscript
# Outcome of every sim of the steady_tracking dataset (speeds 0.5-2.9, 60 reps; fast_speeds run
# 3.2-5.0, 20 reps) as a phase diagram over speed x psi, one panel per reg.
#   aborted   SLiM stopped the sim: the population fell ~1700 behind the optimum and every
#             fitness underflowed to 0 (crashed_sims.csv, from extract_crashed_sims.sh)
#   for the sims that reached tick 12000: the response ratio over ticks 11500-12000,
#             1 - (lag gain per tick) / speed. 1 = the mean keeps pace with the optimum (steady
#             lag); 0 = it no longer moves at all.
#
# Input: figure/r_data/steady_tracking/{data_summary_all,data_summary_fast_speeds,crashed_sims}.*
# Output: analysis/steady_tracking/server/plots_fit500/outcomes/
#   outcomes.pdf   page 1  fraction of reps aborted
#                  page 2  median response ratio of the reps that reached 12000
#                  page 3  lag at tick 12000 vs speed (reps that reached 12000; combinations with
#                          aborted reps drawn hollow: the rest are the ones that kept up best)
#                  page 4  when the aborted sims stopped (ticks after the shift) vs speed
#   outcomes_per_sim.csv, outcomes_summary.csv
#
# Run from the repo root (~1 min):
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_outcomes.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")

out_dir    <- "analysis/steady_tracking/server/plots_fit500/outcomes"
shift_tick <- 10000
end_win    <- c(11500, 12000)   # response ratio window
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_summary(cols = "z_opt", keep_partial = TRUE)
d <- d[tick %in% c(end_win[1]:(end_win[1] + 100), (end_win[2] - 100):end_win[2])]
cr <- read_crashed()

# per surviving sim: lag gain per tick between the first and last 100 ticks of the window
ps <- d[, .(lag_end = mean(z_opt[tick > end_win[2] - 100]),
            gain = (mean(z_opt[tick > end_win[2] - 100]) - mean(z_opt[tick <= end_win[1] + 100])) / (end_win[2] - end_win[1] - 100)),
        by = .(shift_speed, reg, psi, rep, ID)]
ps[, `:=`(response = 1 - gain / shift_speed, outcome = "reached 12000")]
cr_rows <- cr[, .(shift_speed, reg, psi, rep, ID, crash_tick, outcome = "aborted")]
ps <- rbindlist(list(ps, cr_rows), use.names = TRUE, fill = TRUE)
fwrite(ps[order(shift_speed, reg, psi, rep)], file.path(out_dir, "outcomes_per_sim.csv"))

sm <- ps[, .(n = .N, n_aborted = sum(outcome == "aborted"),
             response = median(response, na.rm = TRUE), lag_end = median(lag_end, na.rm = TRUE),
             crash_after = median(crash_tick - shift_tick, na.rm = TRUE)),
         by = .(shift_speed, reg, psi)]
sm[, frac_aborted := n_aborted / n]
fwrite(sm[order(reg, psi, shift_speed)], file.path(out_dir, "outcomes_summary.csv"))

speeds <- sort(unique(sm$shift_speed))
fmt <- function(x) {
  x[, `:=`(speed_f = factor(shift_speed, levels = speeds), psi_f = factor(psi),
           reg = reg_factor(reg))]
  x
}
sm <- fmt(sm); ps <- fmt(ps)
base_theme <- list(theme_bw(base_size = 11),
                   theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92")))
sub_txt <- "speeds up to 2.9: 60 reps (steady_tracking); 3.2 and above: 20 reps (fast_speeds)"

pdf(file.path(out_dir, "outcomes.pdf"), width = 13, height = 5.5)
print(ggplot(sm, aes(psi_f, speed_f, fill = frac_aborted)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = ifelse(n_aborted > 0, sprintf("%d/%d", n_aborted, n), "")), size = 2.6) +
        facet_wrap(~ reg, nrow = 1) +
        scale_fill_gradient("fraction\naborted", low = "grey95", high = "#b2182b", limits = c(0, 1)) +
        labs(x = "psi", y = "speed of the optimum",
             title = "Sims aborted because the population fell ~1700 behind the optimum (all fitness underflowed to 0)",
             subtitle = sub_txt) + base_theme)
print(ggplot(sm, aes(psi_f, speed_f, fill = response)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = ifelse(is.na(response), "all\naborted", sprintf("%.2f", response))), size = 2.4) +
        facet_wrap(~ reg, nrow = 1) +
        scale_fill_gradientn("response\nratio", colours = c("#b2182b", "#f4a582", "#f7f7f7", "#92c5de", "#2166ac"),
                             values = c(0, 0.5, 0.8, 0.95, 1), limits = c(0, 1.05), na.value = "grey60", oob = scales::squish) +
        labs(x = "psi", y = "speed of the optimum",
             title = sprintf("Response ratio over ticks %d-%d (1 = keeps pace, steady lag; < 1 = lag still growing), median of the reps that reached 12000",
                             end_win[1], end_win[2]),
             subtitle = paste0(sub_txt, "; a combination with aborted reps shows only the reps that kept up best")) + base_theme)
le <- sm[!is.na(lag_end)]
print(ggplot(le, aes(shift_speed, lag_end, colour = psi_num(psi_f), group = psi_f)) +
        geom_line(data = le[n_aborted == 0], linewidth = 0.5) +
        geom_point(aes(shape = n_aborted > 0), size = 2) +
        scale_shape_manual(NULL, values = c(`FALSE` = 16, `TRUE` = 1),
                           labels = c("all reps reached 12000", "some reps aborted (survivors only)")) +
        facet_wrap(~ reg, nrow = 1) +
        scale_y_log10() +
        psi_scale(aesthetics = "colour") +
        labs(x = "speed of the optimum", y = "lag at tick 12000 (median, log scale)",
             title = "Lag at the end of the run vs speed", subtitle = sub_txt) + base_theme)
if (nrow(cr)) {
  pc <- ps[outcome == "aborted"]
  print(ggplot(pc, aes(shift_speed, crash_tick - shift_tick, colour = psi_num(psi_f), group = psi_f)) +
          geom_point(position = position_jitter(width = 0.05, height = 0), size = 0.9, alpha = 0.6) +
          stat_summary(fun = median, geom = "line", linewidth = 0.6) +
          geom_function(fun = function(v) 1690 / v, colour = "black", linetype = "dashed", inherit.aes = FALSE) +
          facet_wrap(~ reg, nrow = 1) +
          psi_scale(aesthetics = "colour") +
          labs(x = "speed of the optimum", y = "ticks after the shift when SLiM stopped",
               title = "When the aborted sims stopped (dashed: a population that never moved would stop at 1690 / speed)",
               subtitle = "each point a sim (last recorded tick; the abort is within 20 ticks after); line = median") + base_theme)
}
invisible(dev.off())

print(dcast(sm, reg + psi ~ shift_speed, value.var = "n_aborted"))
print(dcast(sm, reg + psi ~ shift_speed, value.var = "response", fun.aggregate = function(x) round(x, 2)))
cat("plots written to", out_dir, "\n")
