#!/usr/bin/env Rscript
# Does the copying rule change how unevenly reproduction is spread? In a Wright-Fisher model with
# fitnessScaling, an individual's expected offspring number is proportional to its fitness, so the
# variance in offspring number is set by
#   I = var(w) / mean(w)^2   (opportunity for selection) = (fitness_sd / fitness_mean)^2
# Prediction: at the same lag, Copy-the-extreme / Copy-the-best have a higher I than
# Copy-the-average. (The comparison with heterozygosity was dropped: every QTL is under selection
# and there are no neutral sites, so H does not measure drift.)
#
# Input: figure/r_data/steady_tracking/data_summary_all.rds (fitness_mean, fitness_sd, z_opt, z_sd)
# Output: analysis/steady_tracking/server/plots_fit500/fitness_variance/
#   I_vs_lag.pdf        I against lag, one point per sim (end window); page 1 colour = rule,
#                       page 2 panel per psi
#   I_over_time.pdf     mean I over time per combination (panel per rule x speed, colour = psi)
#   zsd_over_time.pdf   phenotypic SD over time: the compression of phenotypes by each rule
#   fitness_variance_end.csv  per-combination means over the end window
#   I_missing_underflow.csv   sims left out of I_vs_lag: I is used only when fitness_mean >= 1e-140 on
#                             every tick of the end window (below that var(fitness) underflows to 0);
#                             in I_over_time a bin is left blank if any tick of any rep is below it
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_fitness_variance.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")   # read_summary(), read_meta(): both runs

# ---- settings ----
data_dir   <- "figure/r_data/steady_tracking"
out_dir    <- "analysis/steady_tracking/server/plots_fit500/fitness_variance"
shift_tick <- 10000
end_window <- c(11500, 12000)   # end state: mean over these ticks
time_bin   <- 20                # ticks per point in the over-time plots

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_summary()
W <- speed_width(uniqueN(d$shift_speed))   # page width: one facet column per speed
d <- d[, .(shift_speed, reg, psi, ID, tick, lag = z_opt, z_sd,
           I = fifelse(fitness_mean > 1e-140, (fitness_sd / fitness_mean)^2, NA_real_))]   # fitness_sd^2 underflows to 0 below ~1e-150: I unknown

speeds <- sort(unique(d$shift_speed))
fmt <- function(x) {
  x[, `:=`(speed_f = factor(shift_speed, levels = speeds, labels = paste("speed", speeds)),
           reg = reg_factor(reg),
           psi_f = factor(psi, labels = paste("psi", sort(unique(psi)))))]
  x
}
reg_cols <- scale_colour_manual(NULL, values = rule_cols)
base_theme <- list(theme_bw(base_size = 11),
                   theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92")))
big_points <- guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

# ---- end state per sim ----
e <- d[tick >= end_window[1] & tick <= end_window[2],
       .(I = mean(I), lag = mean(lag), z_sd = mean(z_sd)), by = .(shift_speed, reg, psi, ID)]
# I only where every tick of the window is usable: a mean over the few ticks above the underflow
# cutoff is noisy and biased to the moments the population came closest to the optimum
no_I <- e[!is.finite(I), .(n_no_I = .N), by = .(shift_speed, reg, psi)]
fwrite(no_I, file.path(out_dir, "I_missing_underflow.csv"))
if (nrow(no_I)) { cat("sims without I (fitness underflowed in the end window):\n"); print(no_I) }
e <- e[is.finite(I) & I > 0]
e <- fmt(e)
miss_txt <- if (nrow(no_I)) paste0("left out (fitness < 1e-140 on some tick, var(fitness) underflows): ",
                                     paste(no_I[, sprintf("%s psi %.1f speed %.1f (%d sims)", reg_names[reg], psi, shift_speed, n_no_I)],
                                           collapse = ", ")) else ""
if (nzchar(miss_txt)) miss_txt <- paste0("\n", paste(strwrap(miss_txt, 110), collapse = "\n"))
sub_e <- sprintf("one point per sim: mean over ticks %d-%d; I = var(fitness) / mean(fitness)^2%s", end_window[1], end_window[2], miss_txt)

pdf(file.path(out_dir, "I_vs_lag.pdf"), width = 10, height = 7)
print(ggplot(e, aes(lag, I, colour = reg)) +
        geom_point(size = 0.7, alpha = 0.5) +
        scale_x_log10() + scale_y_log10() + reg_cols + big_points +
        labs(x = "lag, log scale", y = "variance in relative fitness I, log scale",
             title = "Unevenness of reproduction vs lag, all combinations",
             subtitle = paste(sub_e, "\nprediction: at the same lag, Copy-the-extreme / Copy-the-best above Copy-the-average")) +
        base_theme)
print(ggplot(e, aes(lag, I, colour = reg)) +
        geom_point(size = 0.7, alpha = 0.5) +
        facet_wrap(~ psi_f, nrow = 1) +
        scale_x_log10() + scale_y_log10() + reg_cols + big_points +
        labs(x = "lag, log scale", y = "variance in relative fitness I, log scale",
             title = "Unevenness of reproduction vs lag, by psi (points within a panel: the speeds)", subtitle = sub_e) +
        base_theme)
invisible(dev.off())

# ---- over time ----
tt <- d[tick >= shift_tick - 50][, tbin := tick %/% time_bin * time_bin + time_bin / 2]
tt <- tt[, .(I = mean(I), z_sd = mean(z_sd)), by = .(shift_speed, reg, psi, tbin)]   # NA (gap) if I underflowed on any tick of any rep
tt <- fmt(tt)
over_time <- function(v, ylab, title, log_y = FALSE) {
  p <- ggplot(tt, aes(tbin, .data[[v]], colour = psi_num(psi_f), group = psi_f)) +
    geom_vline(xintercept = shift_tick, linetype = "dotted", colour = "grey40") +
    geom_line(linewidth = 0.5) +
    facet_grid(reg ~ speed_f) + psi_scale(aesthetics = "colour") +
    labs(x = "tick", y = ylab, title = title,
         subtitle = sprintf("mean of 60 reps per %d-tick bin; dotted: optimum starts moving", time_bin)) +
    base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (log_y) p <- p + scale_y_log10()
  p
}
ggsave(file.path(out_dir, "I_over_time.pdf"),
       over_time("I", "variance in relative fitness I (log scale)", "Unevenness of reproduction over time", log_y = TRUE),
       width = W, height = 7)
ggsave(file.path(out_dir, "zsd_over_time.pdf"),
       over_time("z_sd", "SD of phenotype", "Phenotypic SD over time: how much each rule compresses phenotypes"),
       width = W, height = 7)

# ---- table ----
tab <- e[, .(lag = mean(lag), I = mean(I), z_sd = mean(z_sd)), by = .(shift_speed, reg, psi)]
setorder(tab, psi, shift_speed, reg)
fwrite(tab, file.path(out_dir, "fitness_variance_end.csv"))
cat("plots written to", out_dir, "\n")
