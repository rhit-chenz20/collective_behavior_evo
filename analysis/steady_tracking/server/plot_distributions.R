#!/usr/bin/env Rscript
# Heatmaps of the phenotype and genotype distributions over time, one PDF per parameter
# combination (shift_speed x reg x psi), phenotype above genotype on the same y axis.
# Each column of tiles is one recorded generation; tile shade = fraction of that generation's
# individuals (all reps pooled) whose value falls in the bin. Dashed line: the optimum.
#
# Input (from analysis/load_data/server/load_data.R, e.g.
#   -t phenotype,genotype --reps 1:10 --gen-values 9970,9990,10000,10010:10510:20 --shift-speed 0.5,1.4,2.3,2.9):
#   figure/r_data/steady_tracking/{phenotype,genotype}_rep1_10_to10510_every20_selectedSpeed.rds
#   plots_fit500/tracking_class/tracking_class_per_sim.csv   (how many reps failed: mean fitness
#   over the last 500 ticks < 10% of optimum, classify_tracking.R)
# Output: analysis/steady_tracking/server/plots_fit500/distributions/dist_sp_<speed>_reg_<reg>_psi_<psi>.pdf
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_distributions.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- settings ----
data_dir   <- "figure/r_data/steady_tracking"
tag        <- "rep1_10_to10510_every20_selectedSpeed"
out_dir    <- "analysis/steady_tracking/server/plots_fit500/distributions"
shift_tick <- 10000
n_bins     <- 200        # bins over the y range of each combination
pool_reps  <- TRUE       # FALSE: one page per rep in each file instead of all reps pooled
opt_colour <- "#d7301f"  # dashed optimum line (red so it shows on dark tiles)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_type <- function(type) {
  x <- as.data.table(readRDS(file.path(data_dir, sprintf("%s_%s.rds", type, tag))))
  setnames(x, type, "value")
  x[, type := type]
}
d <- rbind(read_type("phenotype"), read_type("genotype"))
d[, type := factor(type, levels = c("phenotype", "genotype"))]
source("analysis/steady_tracking/server/steady_data.R")   # read_class(), keep_speeds()
d <- keep_speeds(d)

m <- read_class()[, .(ID, failed)]

# the optimum moves by shift_speed every tick from shift_tick: (tick - shift_tick + 1) x speed
optimum_at <- function(tick, speed) ifelse(tick >= shift_tick, (tick - shift_tick + 1) * speed, 0)

gens <- sort(unique(d$generation))
gen_levels <- as.character(gens)
# label every other generation so the axis stays readable
gen_breaks <- gen_levels[seq(1, length(gen_levels), by = 2)]

heatmap <- function(x, speed, title, subtitle) {
  lo <- min(x$value); hi <- max(x$value)
  breaks <- seq(lo, hi, length.out = n_bins + 1)
  bw <- breaks[2] - breaks[1]
  x[, bin := pmin(findInterval(value, breaks), n_bins)]
  h <- x[, .N, by = .(type, generation, bin)]
  h[, frac := N / sum(N), by = .(type, generation)]
  h[, `:=`(y = breaks[bin] + bw / 2, gen_f = factor(generation, levels = gen_levels))]
  opt <- data.table(gen_f = factor(gen_levels, levels = gen_levels),
                    optimum = optimum_at(gens, speed))
  opt <- opt[CJ(type = levels(d$type), i = seq_len(nrow(opt)))[, .(type, i)]$i][
    , type := factor(rep(levels(d$type), each = length(gens)), levels = levels(d$type))]

  ggplot(h, aes(gen_f, y)) +
    geom_tile(aes(fill = frac), height = bw) +
    geom_line(data = opt, aes(gen_f, optimum, group = type), colour = opt_colour,
              linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = match(shift_tick, gens), linetype = "dotted", colour = "grey50") +
    facet_wrap(~ type, ncol = 1) +
    scale_fill_gradient("fraction of\nindividuals", low = "grey95", high = "black") +
    scale_x_discrete(breaks = gen_breaks) +
    labs(x = "generation", y = "value", title = title, subtitle = subtitle) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = "grey95"))
}

combos <- unique(d[, .(shift_speed, reg, psi)])[order(shift_speed, reg, psi)]
for (i in seq_len(nrow(combos))) {
  cb <- combos[i]
  x <- d[shift_speed == cb$shift_speed & reg == cb$reg & psi == cb$psi]
  ids <- unique(x$ID)
  n_failed <- m[ID %in% ids, sum(failed)]
  title <- sprintf("speed %s, %s, psi %s", cb$shift_speed, reg_names[cb$reg], cb$psi)
  sub <- sprintf("%d/%d reps failed; dashed: optimum; dotted: optimum starts moving", n_failed, length(ids))
  f <- file.path(out_dir, sprintf("dist_sp_%s_reg_%s_psi_%s.pdf", cb$shift_speed, cb$reg, cb$psi))
  pdf(f, width = 10, height = 9)
  if (pool_reps) {
    print(heatmap(copy(x), cb$shift_speed, title, paste(sprintf("%d reps pooled;", length(ids)), sub)))
  } else {
    for (r in sort(unique(x$rep))) {
      xr <- x[rep == r]
      print(heatmap(copy(xr), cb$shift_speed, sprintf("%s, rep %d", title, r),
                    sprintf("%s; dashed: optimum", if (m[ID == xr$ID[1], failed]) "failed" else "tracked")))
    }
  }
  invisible(dev.off())
  cat(sprintf("[%d/%d] %s\n", i, nrow(combos), basename(f)))
}
