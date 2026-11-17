#!/usr/bin/env Rscript
# Each flip of the fluctuating optimum as an abrupt shift: does a population that lives with
# fluctuation respond to a flip differently from a population shifted once after a long constant
# environment?
#   fluctuating (fluc_env_adjusted_shift_size): every flip at tick f in [10000, 11000 - half
#     period]; progress(tau) = (z_mean(f + tau) - old optimum) / (new optimum - old optimum),
#     tau = 0 ... half period (until the next flip); both flip directions pooled.
#   abrupt (cons_env_data_with_adjusted_shift_size): optimum 0 -> shift_size at tick 10000;
#     progress(tau) = z_mean(10000 + tau) / shift_size.
#   Matched distances: shift_size 9.4 and 28.2 = fluc_opt 9.4 and 28.2.
#   progress 0 = still at the old optimum, 1 = at the new optimum. A fluctuating population can
#   start above 0 (it had not reached the old optimum when the flip came).
# Outputs (analysis/fluctuating_env/server/plots/flip_response/):
#   flip_response_dist<d>_n8.pdf   rows psi, columns rule; mean progress over flips and reps vs
#                                  generations since the flip; black = abrupt shift, colours =
#                                  fluctuation period
#   flip_t50_n8.pdf                generations to reach progress 0.5 (from the mean curve), all
#                                  fluctuating amplitudes x periods as a heatmap; NA (grey) = not
#                                  reached before the next flip
#   flip_response_n8.csv           the mean curves
#
# Inputs:
#   figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_track_gen9950-11000.rds
#   figure/r_data/cons_env_data_with_adjusted_shift_size/data_summary_n2_8_20_50_allpsi_gen9950-11000.rds
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/plot_flip_response.R

source("analysis/fluctuating_env/server/fluc_heatmap_helpers.R")

fluc_file <- "figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_track_gen9950-11000.rds"
cons_file <- "figure/r_data/cons_env_data_with_adjusted_shift_size/data_summary_n2_8_20_50_allpsi_gen9950-11000.rds"
out_dir   <- "analysis/fluctuating_env/server/plots/flip_response"
psi_set   <- c(0, 0.2, 0.5, 0.8)
dists     <- c(9.4, 28.2)
t_from    <- 10000
t_end     <- 11000
max_tau   <- 110            # longest half period
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- fluctuating: cut every flip ----
s <- as.data.table(readRDS(fluc_file))[, .(reg, fluc_opt, fluc_interval, psi, rep, ID, tick, optimum, z_mean)]
setkey(s, ID, tick)
s[, prev_opt := shift(optimum), by = ID]
fl <- s[tick >= t_from & !is.na(prev_opt) & optimum != prev_opt & tick + fluc_interval - 1 <= t_end,
        .(ID, f = tick, new = optimum, old = prev_opt)]
cut <- s[fl, on = .(ID), allow.cartesian = TRUE][tick >= f & tick < f + fluc_interval]
cut[, `:=`(tau = tick - f, progress = (z_mean - old) / (new - old))]
fc <- cut[, .(progress = mean(progress), flips = .N), by = .(reg, fluc_opt, fluc_interval, psi, tau)]
cat(sprintf("fluctuating: %d flips cut (%d-%d per sim)\n", nrow(fl), min(fl[, .N, by = ID]$N), max(fl[, .N, by = ID]$N)))

# ---- abrupt shift ----
cn <- as.data.table(readRDS(cons_file))[n == 8 & psi %in% psi_set & tick >= t_from & tick <= t_from + max_tau]
cn[, `:=`(tau = tick - t_from, progress = z_mean / shift_size)]
ac <- cn[, .(progress = mean(progress), reps = uniqueN(ID)), by = .(reg, shift_size, psi, tau)]
cat(sprintf("abrupt: %d-%d reps per curve\n", min(ac$reps), max(ac$reps)))

fwrite(rbind(fc[, .(env = "fluctuating", reg, distance = fluc_opt, period = 2L * fluc_interval, psi, tau, progress)],
             ac[, .(env = "abrupt", reg, distance = shift_size, period = NA_integer_, psi, tau, progress)]),
       file.path(out_dir, "flip_response_n8.csv"))

# ---- curves at the matched distances ----
periods <- sort(unique(2L * fc$fluc_interval))
per_cols <- setNames(c("#F2A541", "#E0402B", "#9B2D86", "#3B1F6B"), paste("period", periods))
lab_psi <- function(d) d[, psi_lab := factor(sprintf("psi == %s", format(psi)),
                                             levels = sprintf("psi == %s", format(psi_set)))]
for (dd in dists) {
  f1 <- fc[abs(fluc_opt - dd) < 1e-9][, `:=`(curve = paste("period", 2L * fluc_interval), reg = reg_factor(reg))]
  a1 <- ac[abs(shift_size - dd) < 1e-9][, `:=`(curve = "abrupt shift", reg = reg_factor(reg))]
  lab_psi(f1); lab_psi(a1)
  p <- ggplot(mapping = aes(tau, progress, colour = curve)) +
    geom_hline(yintercept = c(0, 1), colour = "grey80", linewidth = 0.3) +
    geom_hline(yintercept = 0.5, colour = "grey80", linewidth = 0.3, linetype = "dotted") +
    geom_line(data = a1, linewidth = 0.5) +
    geom_line(data = f1, linewidth = 0.45) +
    facet_grid(psi_lab ~ reg, switch = "y", labeller = labeller(psi_lab = label_parsed, reg = label_value)) +
    scale_colour_manual(NULL, values = c("abrupt shift" = "black", per_cols),
                        breaks = c("abrupt shift", names(per_cols))) +
    scale_x_continuous(limits = c(0, max_tau), expand = expansion(0.01)) +
    labs(x = "Generations since the flip / shift",
         y = "Progress toward the new optimum (0 = old, 1 = new)",
         title = sprintf("Response to a change of the optimum by %s (%d sd): fluctuating flips vs one abrupt shift",
                         format(dd), as.integer(round(dd / sd_unit)))) +
    theme_bw(base_size = 7) +
    theme(panel.grid = element_blank(), strip.background = element_blank(), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90), legend.position = "bottom",
          plot.title = element_text(size = 7, hjust = 0.5))
  f <- file.path(out_dir, sprintf("flip_response_dist%s_n8.pdf", format(dd)))
  ggsave(f, p, width = 6, height = 6)
  cat("  written", f, "\n")
}

# ---- t50 heatmap over all amplitudes and periods ----
t50 <- function(tau, pr) { i <- which(pr >= 0.5); if (length(i)) tau[min(i)] else NA_real_ }
ft <- fc[order(tau), .(val = t50(tau, progress), start = progress[tau == 0]), by = .(reg, fluc_opt, fluc_interval, psi)]
at <- ac[order(tau), .(t50 = t50(tau, progress)), by = .(reg, shift_size, psi)]
cat("abrupt-shift t50 (generations):\n"); print(dcast(at, reg + psi ~ shift_size, value.var = "t50"))
p <- hm_plot(ft, "Generations for the mean phenotype to get halfway to the new optimum after a flip",
             fill = "seq", legend_title = NULL) +
  scale_fill_viridis_c(NULL, option = "magma", direction = -1, na.value = "grey80",
                       guide = guide_colourbar(direction = "horizontal", barwidth = unit(4, "cm"),
                                               barheight = unit(0.13, "cm")))
suppressMessages(ggsave(file.path(out_dir, "flip_t50_n8.pdf"), p, width = fig_w, height = fig_h))
cat(sprintf("t50 range %s - %s; %d tiles never reached 0.5\n", min(ft$val, na.rm = TRUE),
            max(ft$val, na.rm = TRUE), sum(is.na(ft$val))))
