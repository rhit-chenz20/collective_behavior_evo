#!/usr/bin/env Rscript
# Is selection acting on whole groups under fit / ext ("one individual sets the group's fate")
# and on individuals under ave? For every sim and recorded generation, each individual's fitness
# is recomputed from its phenotype, w = dnorm(phenotype, optimum, sqrt(Vs)), and the variance in
# relative fitness (w / mean w) is split into a between-group and a within-group part, using the
# groups of that generation (group_tag; 125 groups of 8, reshuffled every generation):
#   I = var(w_rel) = var_between (of group means) + var_within (around group means)
# The same split is done for the phenotype and for the genotypic value.
#
# Input (load_data.R ... --reps 1:10 --shift-speed 0.5,1.4,2.3,2.9
#        --gen-values 9970,9990,10000,10010:10510:20 --tag rep1_10_to10510_every20_selectedSpeed):
#   figure/r_data/steady_tracking/{phenotype,genotype,group_tag}_<tag>.rds
#   figure/r_data/steady_tracking/data_summary_all.rds   (only to check I against fitness_sd/mean)
# Output: analysis/steady_tracking/server/plots_fit500/fitness_between_within/
#   fitness_between_share.pdf  share of the variance in relative fitness that is between groups,
#                              over time (panel per reg x speed, line per psi)
#   fitness_var_parts.pdf      I split into between- and within-group parts (stacked), panel per
#                              reg x psi, one page per speed
#   pheno_geno_between_share.pdf  the same between-group share for phenotype and genotypic value
#   between_within_check.pdf   I from individuals vs I from the data file (should lie on y = x)
#   between_within_per_gen.csv per sim and generation: I, parts, shares
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/steady_tracking/server/plot_fitness_between_within.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- settings ----
data_dir   <- "figure/r_data/steady_tracking"
tag        <- "rep1_10_to10510_every20_selectedSpeed"
out_dir    <- "analysis/steady_tracking/server/plots_fit500/fitness_between_within"
shift_tick <- 10000
Vs         <- 2000     # = 2 N in base_script.slim

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("analysis/steady_tracking/server/steady_data.R")   # keep_speeds(): plot_speeds
rd <- function(type) keep_speeds(as.data.table(readRDS(file.path(data_dir, sprintf("%s_%s.rds", type, tag)))))
keys <- c("shift_speed", "reg", "psi", "rep", "ID", "generation", "ind")
x <- rd("phenotype")[, c(keys, "phenotype"), with = FALSE]
x <- merge(x, rd("genotype")[, .(ID, generation, ind, genotype)], by = c("ID", "generation", "ind"))
x <- merge(x, rd("group_tag")[, .(ID, generation, ind, group_tag)], by = c("ID", "generation", "ind"))
cat(sprintf("%s individual records, %d sims\n", format(nrow(x), big.mark = ","), uniqueN(x$ID)))

# optimum at each recorded tick: moves by shift_speed per tick from shift_tick
x[, optimum := fifelse(generation >= shift_tick, (generation - shift_tick + 1) * shift_speed, 0)]
x[, w := dnorm(phenotype, optimum, sqrt(Vs))]
x[, w_rel := w / mean(w), by = .(ID, generation)]

# variance split (population variances, as for fitness_sd in SLiM's log)
pvar <- function(v) mean((v - mean(v))^2)
split_var <- function(v, grp) {
  gm <- ave(v, grp)
  c(total = pvar(v), between = pvar(gm), within = mean((v - gm)^2))
}
s <- x[, {
  f <- split_var(w_rel, group_tag); z <- split_var(phenotype, group_tag); g <- split_var(genotype, group_tag)
  .(I = f[["total"]], I_between = f[["between"]], I_within = f[["within"]],
    z_var = z[["total"]], z_between = z[["between"]], g_var = g[["total"]], g_between = g[["between"]])
}, by = .(shift_speed, reg, psi, rep, ID, generation)]
s[, `:=`(share_fit = I_between / I, share_z = z_between / z_var, share_g = g_between / g_var)]
fwrite(s, file.path(out_dir, "between_within_per_gen.csv"))

speeds <- sort(unique(s$shift_speed))
s[, `:=`(speed_f = factor(shift_speed, levels = speeds, labels = paste("speed", speeds)),
         reg = reg_factor(reg),
         psi_f = factor(psi))]
psi_cols <- psi_scale(aesthetics = "colour")
base_theme <- list(
  geom_vline(xintercept = shift_tick, linetype = "dotted", colour = "grey40"),
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_text(angle = 45, hjust = 1))
)
m <- s[, lapply(.SD, mean), by = .(reg, speed_f, psi_f, generation),
       .SDcols = c("I", "I_between", "I_within", "share_fit", "share_z", "share_g")]

# ---- 1. between-group share of fitness variance ----
# with random groups of 8 and no copying, the between-group share of any individual trait is 1/8
p <- ggplot(m, aes(generation, share_fit, colour = psi_num(psi_f), group = psi_f)) +
  geom_hline(yintercept = 1 / 8, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.6) +
  facet_grid(reg ~ speed_f) + psi_cols +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "generation", y = "between-group share of the variance in relative fitness",
       title = "How much of the fitness differences are between whole groups",
       subtitle = "mean of 10 reps; dashed: 1/8, the share expected with random groups and no copying") +
  base_theme
ggsave(file.path(out_dir, "fitness_between_share.pdf"), p, width = 14, height = 8)

# ---- 2. I split into its parts ----
parts <- melt(m, id.vars = c("reg", "speed_f", "psi_f", "generation"),
              measure.vars = c("I_between", "I_within"), variable.name = "part", value.name = "var")
parts[, part := factor(part, levels = c("I_within", "I_between"), labels = c("within groups", "between groups"))]
pdf(file.path(out_dir, "fitness_var_parts.pdf"), width = 13, height = 8)
for (sp in levels(parts$speed_f)) {
  print(ggplot(parts[speed_f == sp], aes(generation, var, fill = part)) +
          geom_area(position = "stack", colour = NA) +
          facet_grid(reg ~ psi_f, labeller = labeller(psi_f = function(v) paste("psi", v)), scales = "free_y") +
          scale_fill_manual(NULL, values = c("within groups" = "grey70", "between groups" = "#b2182b")) +
          labs(x = "generation", y = "variance in relative fitness I (stacked parts)",
               title = sprintf("Variance in relative fitness split into between- and within-group parts, %s", sp),
               subtitle = "mean of 10 reps; y axis free per row") +
          base_theme)
}
invisible(dev.off())

# ---- 3. the same split for phenotype and genotype ----
pg <- melt(m, id.vars = c("reg", "speed_f", "psi_f", "generation"),
           measure.vars = c("share_z", "share_g"), variable.name = "trait", value.name = "share")
pg[, trait := factor(trait, levels = c("share_z", "share_g"), labels = c("phenotype", "genotypic value"))]
p <- ggplot(pg, aes(generation, share, colour = psi_num(psi_f), linetype = trait, group = interaction(psi_f, trait))) +
  geom_hline(yintercept = 1 / 8, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.5) +
  facet_grid(reg ~ speed_f) + psi_cols +
  scale_linetype_manual(NULL, values = c(phenotype = "solid", "genotypic value" = "dotted")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "generation", y = "between-group share of the variance",
       title = "Between-group share of the variance in phenotype (solid) and genotypic value (dotted)",
       subtitle = "genotypic value stays near 1/8 (groups are random); copying moves the phenotype share") +
  base_theme
ggsave(file.path(out_dir, "pheno_geno_between_share.pdf"), p, width = 14, height = 8)

# ---- 4. check against the data file ----
ids <- unique(s$ID)   # outside the [ ]: the data file has its own column called s
dd <- keep_speeds(as.data.table(readRDS(file.path(data_dir, "data_summary_all.rds"))))[
  ID %in% ids, .(ID, generation = tick, I_data = (fitness_sd / fitness_mean)^2)]
ck <- merge(s[, .(ID, generation, I, reg)], dd, by = c("ID", "generation"))
cat(sprintf("check: %d sim-generations, median |I - I_data| / I_data = %.2g\n",
            nrow(ck), median(abs(ck$I - ck$I_data) / ck$I_data)))
p <- ggplot(ck, aes(I_data, I, colour = reg)) +
  geom_abline(linetype = "dashed") + geom_point(size = 0.5, alpha = 0.4) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "I from the data file (fitness_sd / fitness_mean)^2", y = "I recomputed from individuals",
       title = "Check: fitness recomputed from phenotype and optimum") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "between_within_check.pdf"), p, width = 7, height = 6)

print(s[generation >= 10310, .(share_fit = round(mean(share_fit), 2), share_z = round(mean(share_z), 2),
                               I = signif(mean(I), 2)), by = .(shift_speed, reg, psi)][order(shift_speed, psi, reg)][
  shift_speed %in% c(1.4, 2.9)])
cat("plots written to", out_dir, "\n")
