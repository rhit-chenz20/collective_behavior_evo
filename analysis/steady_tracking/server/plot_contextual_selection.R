#!/usr/bin/env Rscript
# Contextual analysis: how much of the selection on genotypic value acts through the individual's
# own genotype and how much through its group? For every sim and recorded generation:
#
#   model A (standard contextual analysis, Heisler & Damuth 1987):
#     w_rel = a + b_own * g + b_group * g_bar + e
#     g = own genotypic value, g_bar = mean genotypic value of the group (incl. self)
#   model B (rule-specific context):
#     w_rel = a + b_own * g + b_copied * c + e
#     c = the value the rule copies: mean of the other members (ave), genotypic value of the member
#         closest to the optimum (fit), of the member furthest from the group mean (ext)
#
# Both predictors are divided by the population SD of g, so coefficients are the change in
# relative fitness per genetic SD. The selection differential on genotype, S = cov(w_rel, g) / sd(g),
# splits exactly into an individual part b_own * var(g)/sd(g)^2 = b_own and a group part
# b_group * cov(g, g_bar) / sd(g)^2 (Heisler & Damuth 1987); the group share = group part / S.
# Groups are random, so cov(g, g_bar) / var(g) is about 1/8.
#
# Input (load_data.R ... --tag rep1_10_to10510_every20_selectedSpeed):
#   figure/r_data/steady_tracking/{phenotype,genotype,group_tag}_<tag>.rds
# Output: analysis/steady_tracking/server/plots_fit500/contextual_selection/
#   coef_mean_context.pdf      model A: b_own and b_group over time (panel per reg x speed, line per psi)
#   coef_copied_context.pdf    model B: b_own and b_copied over time
#   selection_partition.pdf    S split into individual and group parts (mean over gens in avg_gens)
#   phenotype_rebuild_check.pdf  phenotype recomputed from group genotypes vs the recorded phenotype
#   contextual_per_gen.csv     per sim and generation: all coefficients, S and its parts, R^2
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_contextual_selection.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- settings ----
data_dir   <- "figure/r_data/steady_tracking"
tag        <- "rep1_10_to10510_every20_selectedSpeed"
out_dir    <- "analysis/steady_tracking/server/plots_fit500/contextual_selection"
shift_tick <- 10000
Vs         <- 2000              # = 2 N in base_script.slim
n_group    <- 8
avg_gens   <- c(10310, 10510)   # generations averaged in selection_partition.pdf

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("analysis/steady_tracking/server/steady_data.R")   # keep_speeds(): plot_speeds
rd <- function(type) keep_speeds(as.data.table(readRDS(file.path(data_dir, sprintf("%s_%s.rds", type, tag)))))
x <- rd("phenotype")[, .(shift_speed, reg, psi, rep, ID, generation, ind, phenotype)]
x <- merge(x, rd("genotype")[, .(ID, generation, ind, g = genotype)], by = c("ID", "generation", "ind"))
x <- merge(x, rd("group_tag")[, .(ID, generation, ind, group_tag)], by = c("ID", "generation", "ind"))
cat(sprintf("%s individual records, %d sims\n", format(nrow(x), big.mark = ","), uniqueN(x$ID)))

x[, optimum := fifelse(generation >= shift_tick, (generation - shift_tick + 1) * shift_speed, 0)]
x[, w_rel := { w <- dnorm(phenotype, optimum, sqrt(Vs)); w / mean(w) }, by = .(ID, generation)]

# ---- group context: mean, and the value each rule copies (as in make_model.py) ----
setkey(x, ID, generation, group_tag)
x[, `:=`(g_bar = mean(g),
         c_ave = (sum(g) - g) / (.N - 1),
         c_fit = g[which.max(dnorm(g, optimum[1], sqrt(Vs)))],
         c_ext = g[which.max(abs(g - mean(g)))]),
  by = .(ID, generation, group_tag)]
x[, copied := fcase(reg == "ave", c_ave, reg == "fit", c_fit, reg == "ext", c_ext)]
x[, pheno_rebuilt := fifelse(reg == "ave",
                             ((n_group - 1) * (1 - psi) * g + psi * g_bar * n_group) / (n_group - 1 + psi),
                             (1 - psi) * g + psi * copied)]
x[, c("c_ave", "c_fit", "c_ext") := NULL]
err <- x[, abs(pheno_rebuilt - phenotype)]
cat(sprintf("phenotype rebuilt from group genotypes: median |error| %.2g, 99.9%% quantile %.2g (phenotypes are logged to ~6 digits)\n",
            median(err), quantile(err, 0.999)))

# ---- regressions per sim x generation (OLS with two predictors, closed form) ----
ols2 <- function(y, x1, x2) {
  v <- cov(cbind(y, x1, x2))
  b <- tryCatch(solve(v[2:3, 2:3], v[2:3, 1]), error = function(e) c(NA_real_, NA_real_))
  r2 <- if (anyNA(b)) NA_real_ else sum(b * v[2:3, 1]) / v[1, 1]
  list(b1 = b[[1]], b2 = b[[2]], r2 = r2)
}
cs <- x[, {
  sg <- sqrt(mean((g - mean(g))^2))
  gs <- g / sg; gb <- g_bar / sg; cc <- copied / sg
  A <- ols2(w_rel, gs, gb); B <- ols2(w_rel, gs, cc)
  S <- mean((w_rel - mean(w_rel)) * (gs - mean(gs)))      # selection differential per genetic SD
  cov_g_gb <- mean((gs - mean(gs)) * (gb - mean(gb)))     # ~1/8 with random groups
  var_gs <- mean((gs - mean(gs))^2)                       # = 1
  .(S = S, b_own_A = A$b1, b_group = A$b2, r2_A = A$r2,
    S_ind = A$b1 * var_gs, S_group = A$b2 * cov_g_gb, cov_g_gbar = cov_g_gb,
    b_own_B = B$b1, b_copied = B$b2, r2_B = B$r2)
}, by = .(shift_speed, reg, psi, rep, ID, generation)]
cs[, group_share := S_group / S]
fwrite(cs, file.path(out_dir, "contextual_per_gen.csv"))

speeds <- sort(unique(cs$shift_speed))
fmt <- function(d) {
  d[, `:=`(speed_f = factor(shift_speed, levels = speeds, labels = paste("speed", speeds)),
           reg = reg_factor(reg), psi_f = factor(psi))]
  d
}
cs <- fmt(cs)
base_theme <- list(
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_text(angle = 45, hjust = 1))
)
shift_line <- geom_vline(xintercept = shift_tick, linetype = "dotted", colour = "grey40")
zero_line  <- geom_hline(yintercept = 0, colour = "grey60")

coef_plot <- function(own, ctx, ctx_label, title, sub) {
  m <- cs[, .(own = mean(get(own)), ctx = mean(get(ctx))), by = .(reg, speed_f, psi_f, generation)]
  m <- melt(m, id.vars = c("reg", "speed_f", "psi_f", "generation"), variable.name = "coef")
  m[, coef := factor(coef, levels = c("own", "ctx"), labels = c("own genotype", ctx_label))]
  ggplot(m, aes(generation, value, colour = psi_num(psi_f), linetype = coef, group = interaction(psi_f, coef))) +
    zero_line + shift_line +
    geom_line(linewidth = 0.55) +
    facet_grid(reg ~ speed_f) + psi_scale(aesthetics = "colour") +
    scale_linetype_manual(NULL, values = setNames(c("solid", "dashed"), c("own genotype", ctx_label))) +
    labs(x = "generation", y = "change in relative fitness per genetic SD", title = title, subtitle = sub) +
    base_theme
}
ggsave(file.path(out_dir, "coef_mean_context.pdf"),
       coef_plot("b_own_A", "b_group", "group mean genotype",
                 "Contextual analysis: relative fitness ~ own genotype + group mean genotype",
                 "partial regression coefficients, mean of 10 reps; both predictors in units of the population's genetic SD"),
       width = 14, height = 8)
ggsave(file.path(out_dir, "coef_copied_context.pdf"),
       coef_plot("b_own_B", "b_copied", "copied value",
                 "Relative fitness ~ own genotype + the value the rule copies",
                 "copied value: mean of the others (Copy-the-average), member closest to the optimum (Copy-the-best), member furthest from the group mean (Copy-the-extreme)"),
       width = 14, height = 8)

# ---- partition of the selection differential ----
pt <- cs[generation >= avg_gens[1] & generation <= avg_gens[2],
         .(S = mean(S), S_ind = mean(S_ind), S_group = mean(S_group)), by = .(reg, speed_f, psi_f)]
ptl <- melt(pt, id.vars = c("reg", "speed_f", "psi_f"), measure.vars = c("S_ind", "S_group"), variable.name = "part")
ptl[, part := factor(part, levels = c("S_group", "S_ind"), labels = c("through the group", "through own genotype"))]
pdf(file.path(out_dir, "selection_partition.pdf"), width = 12, height = 8)
print(ggplot(ptl, aes(psi_f, value, fill = part)) +
        zero_line +
        geom_col(width = 0.75) +
        geom_point(data = pt, aes(psi_f, S), inherit.aes = FALSE, shape = 4, size = 2) +
        facet_grid(reg ~ speed_f, scales = "free_y") +
        scale_fill_manual(NULL, values = c("through own genotype" = "grey60", "through the group" = "#b2182b")) +
        labs(x = "psi", y = "selection differential on genotype (relative fitness per genetic SD)",
             title = "Selection on genotype split into an individual and a group part (model A)",
             subtitle = sprintf("mean over generations %d-%d and 10 reps; x = total S; y axis free per row",
                                avg_gens[1], avg_gens[2])) +
        base_theme)
sh <- cs[generation >= avg_gens[1] & generation <= avg_gens[2],
         .(share = sum(S_group) / sum(S)), by = .(reg, speed_f, psi_f)]
print(ggplot(sh, aes(psi_f, share, colour = reg, group = reg)) +
        geom_hline(yintercept = c(0, 1), colour = "grey60") +
        geom_line() + geom_point(size = 2) +
        facet_wrap(~ speed_f, nrow = 1) +
        scale_colour_manual(NULL, values = rule_cols) +
        labs(x = "psi", y = "share of S acting through the group",
             title = "Share of the selection on genotype that acts through the group (model A)",
             subtitle = sprintf("sum of group parts / sum of S over generations %d-%d and 10 reps", avg_gens[1], avg_gens[2])) +
        base_theme)
invisible(dev.off())

# ---- check ----
set.seed(1)
smp <- x[sample(.N, min(.N, 2e5))]
p <- ggplot(smp, aes(phenotype, pheno_rebuilt, colour = reg)) +
  geom_abline(linetype = "dashed") + geom_point(size = 0.3, alpha = 0.3) +
  labs(x = "recorded phenotype", y = "phenotype rebuilt from group genotypes",
       title = "Check: copying rules reproduced from genotypes and groups (200,000 random individuals)") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "phenotype_rebuild_check.pdf"), p, width = 8, height = 7)

tab <- cs[generation >= avg_gens[1] & generation <= avg_gens[2],
          .(S = signif(mean(S), 2), b_own = signif(mean(b_own_A), 2), b_group = signif(mean(b_group), 2),
            group_share = round(sum(S_group) / sum(S), 2),
            b_own_B = signif(mean(b_own_B), 2), b_copied = signif(mean(b_copied), 2),
            r2_A = round(mean(r2_A), 2), r2_B = round(mean(r2_B), 2)),
          by = .(shift_speed, reg, psi)][order(shift_speed, psi, reg)]
fwrite(tab, file.path(out_dir, "contextual_summary.csv"))
print(tab[shift_speed %in% c(1.4, 2.9)])
cat("plots written to", out_dir, "\n")
