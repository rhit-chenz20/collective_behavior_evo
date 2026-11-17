#!/usr/bin/env Rscript
# Fitness load over time under the fluctuating optimum (fluc_env_adjusted_shift_size, n = 8).
#   load = 1 - mean fitness / maximum fitness, maximum = dnorm(0, 0, sqrt(Vs)) (an individual
#   exactly at the optimum); 0 = no loss, 1 = everything lost.
# Rows = the four corners of the amplitude x period grid (as plot_fluc_dynamics.R), columns =
# copying rule; x = cycles since the window start (generations / period, so every row spans 0-2);
# window = n_cycles full cycles from the
# first optimum flip at or after win_from (same window as plot_fitness_heatmap.R "cycles").
# Dots = mean over the 20 reps per generation, ribbon = 95% CI; colour = psi; dotted verticals =
# optimum flips. Each row has its own y range.
#
# Input:  figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_fitness_gen10000-10500.rds
#         (load_data.R -d fluc_env_adjusted_shift_size -t data --n 8 --psi 0,0.2,0.5,0.8
#          --cols optimum,fitness_mean --gens 10000:10500 --tag n8_fitness_gen10000-10500)
# Output: analysis/fluctuating_env/server/plots/fitness_load/fitness_load_n8_corners_8.7cm.pdf
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/plot_fitness_load.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(MetBrewer)
})

# ---- settings ----
data_file <- "figure/r_data/fluc_env_adjusted_shift_size/data_summary_n8_fitness_gen10000-10500.rds"
out_dir   <- "analysis/fluctuating_env/server/plots/fitness_load"
corners   <- data.table(fluc_opt = c(9.4, 9.4, 103.4, 103.4), fluc_interval = c(20L, 110L, 20L, 110L))
win_from  <- 10000
n_cycles  <- 2
psi_show  <- c(0, 0.2, 0.5, 0.8)
psi_grid  <- seq(0, 1, by = 0.1)
Vs        <- 2000
sd_unit   <- 9.4
reg_names <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")

# ---- figure size ----
fig_w <- 8.7 / 2.54
fig_h <- 4.2
txt_pt <- 5
tick_pt <- 4.5
point_size <- 0.35

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fmax <- dnorm(0, 0, sqrt(Vs))
s <- as.data.table(readRDS(data_file))[corners, on = .(fluc_opt, fluc_interval)][psi %in% psi_show]
s[, load := 1 - fitness_mean / fmax]

# window per corner: first flip at or after win_from, n_cycles full cycles
flip_ticks <- s[, .(ID1 = ID[1]), by = .(fluc_opt, fluc_interval)][s, on = .(fluc_opt, fluc_interval)][
  ID == ID1][order(tick)][, .(f = tick[which(diff(optimum) != 0) + 1]), by = .(fluc_opt, fluc_interval)]
win <- flip_ticks[f >= win_from, .(t0 = min(f)), by = .(fluc_opt, fluc_interval)]
win[, t1 := t0 + 2 * n_cycles * fluc_interval - 1]
stopifnot(all(win$t1 <= max(s$tick)))
d <- s[win, on = .(fluc_opt, fluc_interval)][tick >= t0 & tick <= t1]
d[, x := (tick - t0) / (2 * fluc_interval)]
m <- d[, .(mean = mean(load), se = sd(load) / sqrt(.N), reps = .N), by = .(fluc_opt, fluc_interval, reg, psi, x)]
fl <- flip_ticks[win, on = .(fluc_opt, fluc_interval)][f > t0 & f <= t1, .(fluc_opt, fluc_interval, x = (f - t0) / (2 * fluc_interval))]

corner_lab <- function(o, i) factor(sprintf("%d sd, period %d", as.integer(round(o / sd_unit)), 2L * i),
                                    levels = sprintf("%d sd, period %d", as.integer(round(corners$fluc_opt / sd_unit)),
                                                     2L * corners$fluc_interval))
m[, `:=`(corner = corner_lab(fluc_opt, fluc_interval), reg = factor(reg_names[reg], levels = reg_names))]
fl[, corner := corner_lab(fluc_opt, fluc_interval)]
cat(sprintf("%d reps per point; load range %.3f - %.3f\n", min(m$reps), min(m$mean), max(m$mean)))
print(m[, .(mean_load = round(mean(mean), 3)), by = .(corner, reg, psi)][order(corner, reg, psi)])

pal <- as.character(met.brewer("Hokusai2", length(psi_grid)))
p <- ggplot(m[order(psi)], aes(x, mean, colour = psi, fill = psi, group = psi)) +
  geom_vline(data = fl, aes(xintercept = x), linetype = "dotted", colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), colour = NA, alpha = 0.3) +
  geom_point(size = point_size, shape = 16) +
  facet_grid(corner ~ reg, scales = "free_y", switch = "y") +
  scale_colour_gradientn(expression(psi), colours = pal, values = psi_grid, limits = range(psi_grid),
                         breaks = psi_show, labels = format(psi_show, drop0trailing = TRUE),
                         aesthetics = c("colour", "fill"),
                         guide = guide_colourbar(nbin = length(psi_grid), display = "rectangles",
                                                 direction = "horizontal", title.position = "left",
                                                 barwidth = unit(2.5, "cm"), barheight = unit(0.13, "cm"))) +
  scale_x_continuous(expand = expansion(0.01), breaks = c(0, 0.5, 1, 1.5, 2)) +
  labs(x = "Cycles since window start (1 cycle = period)",
       y = "Fitness load (1 - mean fitness / maximum fitness)") +
  theme_bw(base_size = txt_pt) +
  theme(panel.grid = element_blank(), panel.spacing = unit(3, "pt"),
        axis.text = element_text(size = tick_pt), axis.ticks.length = unit(1.2, "pt"),
        strip.background = element_blank(), strip.placement = "outside",
        strip.text.x = element_text(size = 5.5, margin = margin(1, 0, 1.5, 0)),
        strip.text.y.left = element_text(size = txt_pt, angle = 90, margin = margin(0, 1, 0, 0)),
        legend.position = "bottom", legend.title = element_text(size = 6),
        legend.text = element_text(size = tick_pt), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"), plot.margin = margin(2, 4, 2, 2))
f <- file.path(out_dir, "fitness_load_n8_corners_8.7cm.pdf")
ggsave(f, p, width = fig_w, height = fig_h)
cat("written", f, "\n")
