#!/usr/bin/env Rscript
# How much of the variance in relative fitness is due to an individual's own genotype, to the
# genotypes of its group mates on average, and to which member the group copied?
#
# For every sim and recorded generation (1000 individuals, 125 random groups of 8), with
# w = relative fitness, g = own genotypic value, g_bar = group mean genotypic value and
# c = the value the rule copies (ave: mean of the others; fit: the member closest to the
# optimum; ext: the member furthest from the group mean), nested linear regressions give
#   R2_own   = R^2 of  w ~ g                  own genotype
#   R2_mean  = R^2 of  w ~ g + g_bar          + group mates' genotypes on average
#   R2_copy  = R^2 of  w ~ g + g_bar + c      + which member was copied
# and the variance in relative fitness I = var(w) is split into
#   own genotype         I * R2_own
#   group mean           I * (R2_mean - R2_own)
#   who was copied       I * (R2_copy - R2_mean)       (0 for ave: its c is a function of g and g_bar)
#   curvature            I * (1 - R2_copy)            fitness is Gaussian in phenotype, and phenotype
#                                                      is exactly linear in g, g_bar and c, so this is
#                                                      only the non-linearity of the fitness function
# Only the own-genotype part changes mean genotype: the selection differential cov(w, g) uses
# only the part of w correlated with g. Everything else is variance in reproductive success that
# the individual's own genes do not explain, i.e. luck from the point of view of its alleles.
#
# Input (load_data.R ... --tag rep1_10_to10510_every20_selectedSpeed):
#   figure/r_data/steady_tracking/{phenotype,genotype,group_tag}_<tag>.rds
# Output: analysis/steady_tracking/server/plots_fit500/fitness_components/
#   fitness_components.pdf   page 1: shares of the variance in relative fitness (bars sum to 1)
#                            page 2: the same parts in absolute units of I
#                            page 3: own-genotype part vs everything else, over time
#   fitness_components_per_gen.csv, fitness_components_summary.csv
#   example_groups.txt       one group per rule, to follow the computation by hand
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_fitness_components.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- settings ----
data_dir   <- "figure/r_data/steady_tracking"
tag        <- "rep1_10_to10510_every20_selectedSpeed"
out_dir    <- "analysis/steady_tracking/server/plots_fit500/fitness_components"
shift_tick <- 10000
Vs         <- 2000
n_group    <- 8
avg_gens   <- c(10310, 10510)   # generations averaged in the bar plots
example    <- list(speed = 1.4, psi = 0.8, generation = 10310, rep = 1)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("analysis/steady_tracking/server/steady_data.R")   # keep_speeds(): plot_speeds
rd <- function(type) keep_speeds(as.data.table(readRDS(file.path(data_dir, sprintf("%s_%s.rds", type, tag)))))
x <- rd("phenotype")[, .(shift_speed, reg, psi, rep, ID, generation, ind, phenotype)]
x <- merge(x, rd("genotype")[, .(ID, generation, ind, g = genotype)], by = c("ID", "generation", "ind"))
x <- merge(x, rd("group_tag")[, .(ID, generation, ind, group_tag)], by = c("ID", "generation", "ind"))
x[, optimum := fifelse(generation >= shift_tick, (generation - shift_tick + 1) * shift_speed, 0)]
x[, w := { f <- dnorm(phenotype, optimum, sqrt(Vs)); f / mean(f) }, by = .(ID, generation)]
setkey(x, ID, generation, group_tag)
x[, `:=`(g_bar = mean(g),
         c_ave = (sum(g) - g) / (.N - 1),
         c_fit = g[which.max(dnorm(g, optimum[1], sqrt(Vs)))],
         c_ext = g[which.max(abs(g - mean(g)))]),
  by = .(ID, generation, group_tag)]
x[, copied := fcase(reg == "ave", c_ave, reg == "fit", c_fit, reg == "ext", c_ext)]

# ---- worked example: one group per rule ----
ex <- x[shift_speed == example$speed & psi == example$psi & generation == example$generation & rep == example$rep]
sink(file.path(out_dir, "example_groups.txt"))
cat(sprintf("speed %s, psi %s, generation %d, rep %d: optimum %.1f\n",
            example$speed, example$psi, example$generation, example$rep, ex$optimum[1]))
cat("phenotype rule: ave  z = ((n-1)(1-psi) g + psi * sum(group g)) / (n-1+psi)\n",
    "               fit  z = (1-psi) g + psi * g of the member closest to the optimum\n",
    "               ext  z = (1-psi) g + psi * g of the member furthest from the group mean\n")
for (r in c("ave", "fit", "ext")) {
  e <- ex[reg == r]
  grp <- e[group_tag == e$group_tag[1]]
  cat(sprintf("\n== %s: group %d (the population's mean relative fitness is 1)\n", r, grp$group_tag[1]))
  print(grp[, .(g = round(g, 1), group_mean = round(g_bar, 1), copied = round(copied, 1),
                phenotype = round(phenotype, 1), rel_fitness = round(w, 3))], row.names = FALSE)
  # spread of the group-level quantities across the population's 125 groups
  gs <- unique(e[, .(group_tag, g_bar, copied_grp = if (r == "ave") g_bar else copied)])
  cat(sprintf("across the 125 groups: SD of group mean genotype %.1f, SD of the copied value %.1f (SD of g %.1f)\n",
              sd(gs$g_bar), sd(gs$copied_grp), sd(e$g)))
}
sink()

# ---- nested regressions per sim x generation ----
r2 <- function(y, X) {           # R^2 of an OLS fit of y on the columns of X (with intercept)
  X <- cbind(1, X)
  fit <- lm.fit(X, y)
  1 - sum(fit$residuals^2) / sum((y - mean(y))^2)
}
cs <- x[, {
  r_own  <- r2(w, cbind(g))
  r_mean <- r2(w, cbind(g, g_bar))
  r_copy <- if (reg[1] == "ave") r_mean else r2(w, cbind(g, g_bar, copied))
  I <- mean((w - 1)^2)
  .(I = I, R2_own = r_own, R2_mean = r_mean, R2_copy = r_copy,
    S = mean((w - 1) * (g - mean(g))) / sqrt(mean((g - mean(g))^2)))
}, by = .(shift_speed, reg, psi, rep, ID, generation)]
fwrite(cs, file.path(out_dir, "fitness_components_per_gen.csv"))

parts <- cs[, .(own = R2_own, group_mean = R2_mean - R2_own, copied = R2_copy - R2_mean,
                curvature = 1 - R2_copy, I = I, S = S), by = .(shift_speed, reg, psi, ID, generation)]

speeds <- sort(unique(parts$shift_speed))
fmt <- function(d) {
  d[, `:=`(speed_f = factor(shift_speed, levels = speeds, labels = paste("speed", speeds)),
           reg = reg_factor(reg), psi_f = factor(psi))]
  d
}
comp_levels <- c("curvature", "copied", "group_mean", "own")
comp_labels <- c("curvature of the fitness function", "which member was copied",
                 "group mates' genotypes (mean)", "own genotype")
comp_cols <- setNames(c("grey85", "#d95f02", "#7570b3", "#1b9e77"), comp_labels)
base_theme <- list(theme_bw(base_size = 11),
                   theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92")))

w_avg <- parts[generation >= avg_gens[1] & generation <= avg_gens[2]]
# shares: averaged over sim-generations; absolute: share x I, averaged
sh <- w_avg[, .(own = mean(own), group_mean = mean(group_mean), copied = mean(copied), curvature = mean(curvature),
                own_abs = mean(own * I), group_mean_abs = mean(group_mean * I), copied_abs = mean(copied * I),
                curvature_abs = mean(curvature * I), I = mean(I), S = mean(S)),
            by = .(shift_speed, reg, psi)]
fwrite(sh[order(shift_speed, psi, reg)], file.path(out_dir, "fitness_components_summary.csv"))
sh <- fmt(sh)
long <- function(cols) {
  l <- melt(sh, id.vars = c("reg", "speed_f", "psi_f"), measure.vars = cols, variable.name = "part")
  l[, part := factor(sub("_abs$", "", part), levels = comp_levels, labels = comp_labels)]
  l
}
sub_txt <- sprintf("mean over generations %d-%d and 10 reps; nested regressions of relative fitness, added in the order own -> group mean -> copied",
                   avg_gens[1], avg_gens[2])

pdf(file.path(out_dir, "fitness_components.pdf"), width = 13, height = 8)
print(ggplot(long(c("own", "group_mean", "copied", "curvature")), aes(psi_f, value, fill = part)) +
        geom_col(width = 0.8) +
        facet_grid(reg ~ speed_f) +
        scale_fill_manual(NULL, values = comp_cols) +
        labs(x = "psi", y = "share of the variance in relative fitness",
             title = "What fitness differences are made of: own genotype vs group mates vs who was copied",
             subtitle = sub_txt) +
        base_theme + theme(legend.position = "bottom"))
print(ggplot(long(c("own_abs", "group_mean_abs", "copied_abs", "curvature_abs")), aes(psi_f, value, fill = part)) +
        geom_col(width = 0.8) +
        facet_grid(reg ~ speed_f, scales = "free_y") +
        scale_fill_manual(NULL, values = comp_cols) +
        labs(x = "psi", y = "variance in relative fitness I (stacked parts)",
             title = "The same parts in absolute units: the green part is what selection on genotype needs, the rest is extra",
             subtitle = paste0(sub_txt, "; y axis free per row")) +
        base_theme + theme(legend.position = "bottom"))
# over time: own-genotype part vs the rest, absolute
tt <- fmt(parts[, .(own = mean(own * I), rest = mean((1 - own) * I)), by = .(shift_speed, reg, psi, generation)])
tt <- melt(tt, id.vars = c("reg", "speed_f", "psi_f", "generation"), measure.vars = c("own", "rest"), variable.name = "part")
tt[, part := factor(part, labels = c("own genotype", "everything else"))]
print(ggplot(tt, aes(generation, value, colour = psi_num(psi_f), linetype = part, group = interaction(psi_f, part))) +
        geom_vline(xintercept = shift_tick, linetype = "dotted", colour = "grey40") +
        geom_line(linewidth = 0.55) +
        facet_grid(reg ~ speed_f, scales = "free_y") +
        psi_scale(aesthetics = "colour") +
        scale_linetype_manual(NULL, values = c("own genotype" = "solid", "everything else" = "dashed")) +
        labs(x = "generation", y = "variance in relative fitness I (part)",
             title = "Variance in relative fitness explained by own genotype vs everything else, over time",
             subtitle = "mean of 10 reps; y axis free per row") +
        base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
invisible(dev.off())

print(sh[shift_speed %in% c(1.4, 2.9) & psi %in% c(0, 0.8),
         .(shift_speed, psi, reg, own = round(own, 2), group_mean = round(group_mean, 2), copied = round(copied, 2),
           curvature = round(curvature, 2), I = signif(I, 2), I_own = signif(own_abs, 2),
           I_rest = signif(I - own_abs, 2), S = round(S, 2))][order(shift_speed, psi, reg)])
cat("plots written to", out_dir, "\n")
