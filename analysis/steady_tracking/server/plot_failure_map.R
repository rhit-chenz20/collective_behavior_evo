#!/usr/bin/env Rscript
# Which parameter combinations fail to track the moving optimum? Outcome per sim from
# classify_tracking.R: failed = mean fitness over the last 500 ticks < 10% of the fitness at the
# optimum ("too unfit"), or aborted by SLiM (fitness underflowed, lag ~1700).
#
# Input: plots_fit500/tracking_class/tracking_class_per_sim.csv (run classify_tracking.R first)
# Output: analysis/steady_tracking/server/plots_fit500/failure_map/
#   failure_map.pdf
#     page 1  fraction of reps failed, speed x psi, panel per copying rule (k/n in each tile)
#     page 2  the same as curves: fraction failed vs speed, line per psi, panel per rule
#     page 3  where the transition is: per rule x psi, the highest speed at which every rep
#             tracked and the lowest at which every rep failed (bar), and the first speed at which
#             at least half failed (point)
#     page 4  mean fitness over the last 500 ticks of every sim vs speed (log scale), panel per
#             rule x psi; dashed: the 10% cutoff; fitness < 1e-4 and aborted sims drawn at 1e-4
#     page 5  median end fitness per combination, speed x psi, panel per rule
#   failure_map_summary.csv    per combination: n, failed, too_unfit, aborted, fraction failed,
#                              median fit_500 and lag_500
#   failure_transition.csv     the speeds of page 3
#
# Only the steady_tracking dataset (the grid, speeds within plot_speeds in steady_data.R); the long runs repeat 11 of its
# combinations to tick 15000 and are compared in tracking_class_long_vs_12000.csv.
#
# Run from the repo root (a few seconds):
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_failure_map.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")   # plot_root, class_file, reg_factor()

out_dir   <- file.path(plot_root, "failure_map")
cutoff    <- 0.1
fit_floor <- 1e-4      # page 4/5: fitness below this (and aborted sims) drawn at the floor
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p <- keep_speeds(fread(class_file)[dataset == "steady_tracking"])   # plot_speeds
p[, failed := outcome == "failed"]

sm <- p[, .(n = .N, failed = sum(failed), too_unfit = sum(class == "too_unfit"),
            aborted = sum(class == "aborted"),
            fit_500_med = median(fifelse(class == "aborted", 0, fit_500)),
            lag_500_med = median(lag_500, na.rm = TRUE)),
        by = .(shift_speed, reg, psi)]
sm[, frac_failed := failed / n]
fwrite(sm[order(reg, psi, shift_speed)], file.path(out_dir, "failure_map_summary.csv"))

speeds <- sort(unique(sm$shift_speed))
fmt <- function(x) {
  x[, `:=`(rule = reg_factor(reg), speed_f = factor(shift_speed, levels = speeds),
           psi_f = factor(sprintf("%.1f", psi)))]
  x
}
sm <- fmt(sm); p <- fmt(p)

base_theme <- list(
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        plot.title = element_text(face = "bold"))
)
out_cols  <- c(tracked = "#2166ac", failed = "#b2182b")
n_note    <- "60 reps per combination at speeds up to 2.9, 20 at 3.2-3.5"

# ---- transition speeds per rule x psi ----
tr <- sm[, {
  ok  <- shift_speed[frac_failed == 0]; bad <- shift_speed[frac_failed == 1]
  half <- shift_speed[frac_failed >= 0.5]
  .(last_all_tracked = if (length(ok)) max(ok[ok < min(c(bad, Inf))], -Inf) else NA_real_,
    first_all_failed = if (length(bad)) min(bad) else NA_real_,
    first_half_failed = if (length(half)) min(half) else NA_real_)
}, by = .(reg, rule, psi, psi_f)]
tr[!is.finite(last_all_tracked), last_all_tracked := NA]
fwrite(tr[order(reg, psi), .(reg, psi, last_all_tracked, first_half_failed, first_all_failed)],
       file.path(out_dir, "failure_transition.csv"))

pdf(file.path(out_dir, "failure_map.pdf"), width = 13, height = 5.5)

# page 1: heatmap of fraction failed
print(ggplot(sm, aes(psi_f, speed_f, fill = frac_failed)) +
        geom_tile(colour = "white", linewidth = 0.6) +
        geom_text(aes(label = sprintf("%d/%d", failed, n),
                      colour = frac_failed > 0.55), size = 2.6, show.legend = FALSE) +
        facet_wrap(~ rule, nrow = 1) +
        scale_fill_gradient("fraction\nfailed", low = "#f4f4f4", high = "#b2182b", limits = c(0, 1)) +
        scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white")) +
        labs(x = "psi", y = "shift speed",
             title = "Which combinations fail to track the moving optimum",
             subtitle = paste0("failed = mean fitness over the last 500 ticks < 10% of the fitness at the optimum, ",
                               "or aborted by SLiM; label: failed / reps (", n_note, ")")) +
        base_theme + theme(panel.grid = element_blank()))

# page 2: fraction failed vs speed
print(ggplot(sm, aes(shift_speed, frac_failed, colour = psi, group = psi_f)) +
        geom_hline(yintercept = 0.5, colour = "grey70", linetype = "dotted") +
        geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
        facet_wrap(~ rule, nrow = 1) +
        psi_scale(aesthetics = "colour") +
        scale_x_continuous(breaks = speeds[seq(1, length(speeds), 2)]) +
        scale_y_continuous(labels = scales::percent) +
        labs(x = "shift speed", y = "reps failed",
             title = "Fraction of reps that fail to track, by speed",
             subtitle = paste0("failed = mean fitness over the last 500 ticks < 10% of optimum, or aborted; ", n_note)) +
        base_theme)

# page 3: transition speeds
print(ggplot(tr, aes(psi_f, colour = rule)) +
        geom_linerange(aes(ymin = last_all_tracked, ymax = first_all_failed),
                       linewidth = 1.2, alpha = 0.5, position = position_dodge(width = 0.5)) +
        geom_point(aes(y = first_half_failed), size = 2.8, position = position_dodge(width = 0.5)) +
        scale_colour_manual(NULL, values = rule_cols) +
        scale_y_continuous(breaks = speeds) +
        labs(x = "psi", y = "shift speed",
             title = "Where tracking breaks down",
             subtitle = paste("bar: from the highest speed at which every rep tracked to the lowest at which every rep failed;",
                              "point: lowest speed at which at least half the reps failed", sep = "\n")) +
        base_theme + theme(legend.position = "top"))

# page 4: end fitness of every sim
pp <- copy(p)
pp[, fit_plot := pmax(fifelse(class == "aborted", fit_floor, fit_500), fit_floor)]
pp[, shape := fifelse(class == "aborted", "aborted", "reached 12000")]
pp[, outcome_f := factor(fifelse(failed, "failed", "tracked"), levels = c("tracked", "failed"))]
print(ggplot(pp, aes(shift_speed, fit_plot, colour = outcome_f, shape = shape)) +
        geom_hline(yintercept = cutoff, linetype = "dashed", colour = "grey30") +
        geom_jitter(width = 0.06, height = 0, size = 0.7, alpha = 0.5) +
        facet_grid(psi_f ~ rule, labeller = labeller(psi_f = function(x) paste("psi", x))) +
        scale_y_log10(breaks = c(1e-4, 1e-3, 0.01, 0.1, 1), labels = c("<=1e-4", "0.001", "0.01", "0.1", "1")) +
        scale_colour_manual(NULL, values = out_cols) +
        scale_shape_manual(NULL, values = c(`reached 12000` = 16, aborted = 4)) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
               shape = guide_legend(override.aes = list(size = 3))) +
        labs(x = "shift speed", y = "mean fitness / fitness at optimum (last 500 ticks)",
             title = "End fitness of every sim",
             subtitle = "dashed: 10% cutoff; aborted sims (x) and fitness below 1e-4 drawn at 1e-4") +
        base_theme + theme(legend.position = "top"))

# page 5: median end fitness per combination
sm[, lf := log10(pmax(fit_500_med, fit_floor))]
print(ggplot(sm, aes(psi_f, speed_f, fill = lf)) +
        geom_tile(colour = "white", linewidth = 0.6) +
        geom_text(aes(label = ifelse(fit_500_med >= 0.01, sub("^0", "", sprintf("%.2f", fit_500_med)), "<.01"),
                      colour = lf < -2.5), size = 2.3, show.legend = FALSE) +
        facet_wrap(~ rule, nrow = 1) +
        scale_fill_gradient2("median\nend fitness", low = "#b2182b", mid = "#f4f4f4", high = "#2166ac",
                             midpoint = log10(cutoff), limits = c(log10(fit_floor), 0),
                             breaks = -4:0, labels = c("<=1e-4", "0.001", "0.01", "0.1", "1")) +
        scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white")) +
        labs(x = "psi", y = "shift speed",
             title = "Median end fitness per combination",
             subtitle = "mean fitness / fitness at optimum over the last 500 ticks, median over reps (aborted = 0); white = the 10% cutoff; label < .01: below 1%") +
        base_theme + theme(panel.grid = element_blank()))

invisible(dev.off())

cat("transition speeds:\n"); print(tr[order(reg, psi), .(reg, psi, last_all_tracked, first_half_failed, first_all_failed)])
cat("plots written to", out_dir, "\n")
