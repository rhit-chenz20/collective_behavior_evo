#!/usr/bin/env Rscript
# How alleles change during selection on the moving optimum (steady_tracking run).
#
# Input (from analysis/load_data/server/load_data.R -t allele ...):
#   figure/r_data/steady_tracking/allele_rep1_10_to10310_every20_selectedSpeed.rds
#   one row per mutation x recorded generation x sim; fixed mutations stay in the table with
#   frequency 1 and isFixed "T"; a lost mutation simply has no more rows.
# Output: analysis/steady_tracking/server/plots_fit500/alleles/
#   allele_fate.pdf          alleles segregating at t0: fraction fixed / lost / still segregating
#                            by t1, for positive- and negative-effect alleles (panel per reg x speed,
#                            x = psi)
#   allele_dp_vs_effect.pdf  dot plot, one page per reg (rows psi, columns speed): each allele
#                            segregating at t0, x = effect size, y = frequency change t0 -> t1,
#                            colour = fate (all reps pooled)
#   allele_fixations.pdf     mutations fixed since t0 over time, positive vs negative effect
#                            (mean over reps)
#   allele_traj.pdf          one page per reg: frequency trajectories of rep traj_rep, alleles
#                            that reach traj_min_freq at some point, colour = effect size
#
# Run from the repo root, optionally with the tag of the input file (default selectedSpeed:
# speeds 0.5, 1.4, 2.3, 2.9; fastSpeed: 3.5, 3.8, 4.1, written to plots/alleles/fastSpeed/):
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_alleles.R [fastSpeed]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("analysis/steady_tracking/server/steady_data.R")   # reg_factor(): display names of the rules

# ---- settings ----
speed_set     <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(speed_set)) speed_set <- "selectedSpeed"
in_file       <- sprintf("figure/r_data/steady_tracking/allele_rep1_10_to10310_every20_%s.rds", speed_set)
out_dir       <- file.path("analysis/steady_tracking/server/plots_fit500/alleles",
                           if (speed_set == "selectedSpeed") "" else speed_set)
t0            <- 10000   # optimum starts moving (reference generation)
t1            <- NULL    # end of the window; NULL = last generation in the file
traj_rep      <- 1       # rep shown in the trajectory plot
traj_min_freq <- 0.05    # trajectory plot: alleles whose frequency reaches this at some point
fix_ylim      <- c(0, 360) # y range of allele_fixations.pdf, the same for every speed set so they compare

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

a <- keep_speeds(as.data.table(readRDS(in_file)))   # plot_speeds (steady_data.R)
if (!nrow(a)) stop(in_file, ": no speed within plot_speeds")
if (is.null(t1)) t1 <- max(a$gen)
a[, `:=`(fixed = isFixed == "T",
         sign = factor(ifelse(effect_size > 0, "positive effect", "negative effect"),
                       levels = c("positive effect", "negative effect")))]
speeds <- sort(unique(a$shift_speed))
a[, `:=`(psi_f = factor(psi),
         speed_f = factor(shift_speed, levels = speeds, labels = paste("speed", speeds)),
         reg = reg_factor(reg))]
cat(sprintf("%s rows, %d sims, generations %d-%d; window t0 = %d -> t1 = %d\n",
            format(nrow(a), big.mark = ","), uniqueN(a$ID), min(a$gen), max(a$gen), t0, t1))

fate_cols <- c(fixed = "#b2182b", segregating = "grey55", lost = "#2166ac")
base_theme <- list(
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"))
)
psi_lab <- labeller(psi_f = function(x) paste("psi", x))

# ---- fate of the alleles segregating at t0 ----
start <- a[gen == t0 & !fixed, .(shift_speed, speed_f, reg, psi_f, rep, ID, mut_id, effect_size, sign, p0 = frequencies)]
end   <- a[gen == t1, .(ID, mut_id, p1 = frequencies, fixed1 = fixed)]
f <- merge(start, end, by = c("ID", "mut_id"), all.x = TRUE)
f[, fate := factor(fifelse(is.na(p1), "lost", fifelse(fixed1, "fixed", "segregating")),
                   levels = c("fixed", "segregating", "lost"))]
f[is.na(p1), p1 := 0]
f[, dp := p1 - p0]

fs <- f[, .N, by = .(reg, speed_f, psi_f, sign, fate)]
fs[, frac := N / sum(N), by = .(reg, speed_f, psi_f, sign)]
p <- ggplot(fs, aes(psi_f, frac, fill = fate)) +
  geom_col(width = 0.8) +
  facet_grid(reg + sign ~ speed_f) +
  scale_fill_manual(NULL, values = fate_cols) +
  labs(x = "psi", y = "fraction of alleles segregating at the shift",
       title = sprintf("Fate by tick %d of the alleles segregating at tick %d (10 reps pooled)", t1, t0)) +
  base_theme
ggsave(file.path(out_dir, "allele_fate.pdf"), p, width = 12, height = 13)

# ---- dot plot: frequency change vs effect size ----
pdf(file.path(out_dir, "allele_dp_vs_effect.pdf"), width = 14, height = 11)
for (g in levels(a$reg)) {
  p <- ggplot(f[reg == g], aes(effect_size, dp, colour = fate)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(size = 0.4, alpha = 0.4, shape = 16) +
    facet_grid(psi_f ~ speed_f, labeller = psi_lab) +
    scale_colour_manual(NULL, values = fate_cols,
                        guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(x = "effect size", y = sprintf("frequency change, tick %d -> %d", t0, t1),
         title = sprintf("Frequency change vs effect size of alleles segregating at the shift, %s", g),
         subtitle = "every allele of 10 reps; lost alleles end at frequency 0, fixed at 1") +
    base_theme
  print(p)
}
invisible(dev.off())

# ---- fixations since t0 over time, by effect sign ----
fixed_t0 <- a[gen == t0 & fixed, .(ID, mut_id)]
fx <- a[fixed == TRUE & gen >= t0][!fixed_t0, on = .(ID, mut_id)]
fx <- fx[, .N, by = .(reg, speed_f, psi_f, ID, gen, sign)]
# fill in zeros (sims / generations / signs with no new fixations yet)
grid <- CJ(ID = unique(a$ID), gen = sort(unique(a[gen >= t0]$gen)), sign = levels(a$sign))
grid <- merge(grid, unique(a[, .(ID, reg, speed_f, psi_f)]), by = "ID")
grid[, sign := factor(sign, levels = levels(a$sign))]
fx <- merge(grid, fx, by = c("ID", "gen", "sign", "reg", "speed_f", "psi_f"), all.x = TRUE)
fx[is.na(N), N := 0]
fxs <- fx[, .(mean = mean(N)), by = .(reg, speed_f, psi_f, gen, sign)]
p <- ggplot(fxs, aes(gen, mean, colour = psi_num(psi_f), linetype = sign, group = interaction(psi_f, sign))) +
  geom_line(linewidth = 0.6) +
  facet_grid(reg ~ speed_f) +
  psi_scale(aesthetics = "colour") +
  scale_linetype_manual(NULL, values = c("positive effect" = "solid", "negative effect" = "dashed")) +
  coord_cartesian(ylim = fix_ylim) +
  labs(x = "generation", y = sprintf("mutations fixed since tick %d (mean of reps)", t0),
       title = "New fixations after the shift: positive- vs negative-effect alleles") +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "allele_fixations.pdf"), p, width = 14, height = 8)

# ---- frequency trajectories of one rep ----
keep <- a[rep == traj_rep, .(maxf = max(frequencies), minf = min(frequencies)), by = .(ID, mut_id)][
  maxf >= traj_min_freq & minf < 1]                  # skip alleles fixed throughout the window
tr <- a[rep == traj_rep][keep, on = .(ID, mut_id)]
lim <- max(abs(tr$effect_size))
pdf(file.path(out_dir, "allele_traj.pdf"), width = 14, height = 11)
for (g in levels(a$reg)) {
  p <- ggplot(tr[reg == g], aes(gen, frequencies, group = mut_id, colour = effect_size)) +
    geom_vline(xintercept = t0, linetype = "dotted", colour = "grey40") +
    geom_line(linewidth = 0.3, alpha = 0.7) +
    facet_grid(psi_f ~ speed_f, labeller = psi_lab) +
    scale_colour_gradient2("effect size", low = "#2166ac", mid = "grey85", high = "#b2182b",
                           midpoint = 0, limits = c(-lim, lim)) +
    labs(x = "generation", y = "allele frequency",
         title = sprintf("Allele frequency trajectories, rep %d, %s", traj_rep, g),
         subtitle = sprintf("alleles reaching frequency %.2f at some point; a line that stops = allele lost; dotted: optimum starts moving",
                            traj_min_freq)) +
    base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}
invisible(dev.off())
cat("plots written to", out_dir, "\n")
