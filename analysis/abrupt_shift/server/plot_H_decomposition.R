#!/usr/bin/env Rscript
# Abrupt shift (cons_env_data_with_adjusted_shift_size): the optimum jumps once, by shift_size,
# at tick 10000; the run ends at 11000. Where does heterozygosity go after the jump?
#   H = sum_i 2 p_i (1 - p_i) / maxL      (maxL = 10000)
# so H splits into the contributions of groups of alleles, and
#   H = (number of segregating alleles) x (mean 2pq per allele) / maxL.
# Same analysis as analysis/steady_tracking/server/plot_H_decomposition.R (moving optimum).
#
# Input (load_data.R -d cons_env_data_with_adjusted_shift_size -t allele --n 8 --psi 0,0.2,0.5,0.8
#        --reps 1:10 --gen-values 9990,10000:10004,10005:10985:20,10995 --tag n8_rep1_10_every20):
#   figure/r_data/cons_env_data_with_adjusted_shift_size/allele_n8_rep1_10_every20.rds
# Output: analysis/abrupt_shift/server/plots/H_decomposition/shift_<size>/
#         (every plot: rows = copying rule, columns = psi; file names end in _shift_<size>)
#   H_by_origin_sign_shift_<size>.pdf  H stacked by allele group: standing (arose before the
#                                      shift) vs new (arose after), positive vs negative effect
#   H_by_freq_class_shift_<size>.pdf   H stacked by the allele's current frequency class
#   H_count_x_mean_shift_<size>.pdf    number of segregating alleles, and mean 2pq per allele
#   fixation_rate_shift_<size>.pdf     fixations per tick between recorded generations, by sign
#   effect_of_segregating_shift_<size>.pdf  mean |effect| of segregating alleles, standing vs new
#   H_parts_start_end_shift_<size>.csv H and its parts at the shift and at the last tick
# All curves are means over the reps of each combination.
#
# Run from the repo root (~1 min):
#   /programs/R-4.5.2/bin/Rscript analysis/abrupt_shift/server/plot_H_decomposition.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- settings ----
in_file    <- "figure/r_data/cons_env_data_with_adjusted_shift_size/allele_n8_rep1_10_every20.rds"
out_root   <- "analysis/abrupt_shift/server/plots/H_decomposition"
shift_tick <- 10000
maxL       <- 10000     # H = sum 2pq / maxL (maxL in cons_*.slim)
freq_breaks <- c(0, 0.05, 0.2, 0.5, 0.8, 0.95, 1)
# display names of the copying rules (the data keep ave / fit / ext)
reg_names <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")

all_a <- as.data.table(readRDS(in_file))
for (size in sort(unique(all_a$shift_size))) {
a <- all_a[shift_size == size]
size_str <- format(size)
out_dir <- file.path(out_root, sprintf("shift_%s", size_str))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(name, ext = "pdf") file.path(out_dir, sprintf("%s_shift_%s.%s", name, size_str, ext))
cat(sprintf("shift size %s: %s rows, %d sims\n", size_str, format(nrow(a), big.mark = ","), uniqueN(a$ID)))
a[, reg := factor(reg_names[as.character(reg)], levels = reg_names)]
a[, psi_f := factor(psi, labels = paste("psi", sort(unique(psi))))]
nrep <- a[, uniqueN(ID), by = .(reg, psi_f)]

seg <- a[isFixed == "F"]
seg[, `:=`(h = 2 * frequencies * (1 - frequencies) / maxL,
           origin = ifelse(originTick < shift_tick, "standing", "new"),
           sign = ifelse(effect_size > 0, "positive", "negative"))]
seg[, group := factor(paste(origin, sign),
                      levels = c("standing positive", "standing negative", "new positive", "new negative"))]
seg[, fclass := cut(frequencies, freq_breaks, include.lowest = TRUE)]

base_theme <- list(
  geom_vline(xintercept = shift_tick, linetype = "dotted", colour = "grey40"),
  facet_grid(reg ~ psi_f),
  theme_bw(base_size = 11),
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_text(angle = 45, hjust = 1))
)
title_sfx <- sprintf("abrupt shift of %s, n = 8, mean of %d reps", size_str, max(nrep$V1))

# per sim sums, with zeros filled in for groups that are empty in a sim/generation, then the
# mean over reps
sum_by <- function(dt, by_col, levels) {
  s <- dt[, .(H = sum(h)), by = c("reg", "psi_f", "ID", "gen", by_col)]
  grid <- CJ(ID = unique(a$ID), gen = sort(unique(a$gen)), g = levels)
  setnames(grid, "g", by_col)
  grid <- merge(grid, unique(a[, .(ID, reg, psi_f)]), by = "ID")
  s <- merge(grid, s, by = c("ID", "gen", by_col, "reg", "psi_f"), all.x = TRUE)
  s[is.na(H), H := 0]
  s <- s[, .(H = mean(H)), by = c("reg", "psi_f", "gen", by_col)]
  s[, (by_col) := factor(get(by_col), levels = levels)]
  s
}

# ---- 1. H by origin x effect sign ----
s1 <- sum_by(seg, "group", levels(seg$group))
p <- ggplot(s1, aes(gen, H, fill = group)) +
  geom_area(position = "stack", colour = NA) +
  scale_fill_manual(NULL, values = c("standing positive" = "#b2182b", "standing negative" = "#2166ac",
                                     "new positive" = "#f4a582", "new negative" = "#92c5de")) +
  labs(x = "tick", y = "heterozygosity H (stacked contributions)",
       title = sprintf("Who carries H: alleles present at the shift (standing) vs arisen after it (new), %s", title_sfx)) +
  base_theme
ggsave(out("H_by_origin_sign"), p, width = 13, height = 8)

# ---- 2. H by current frequency class ----
s2 <- sum_by(seg, "fclass", levels(seg$fclass))
p <- ggplot(s2, aes(gen, H, fill = fclass)) +
  geom_area(position = "stack", colour = NA) +
  scale_fill_brewer("allele frequency", palette = "YlGnBu") +
  labs(x = "tick", y = "heterozygosity H (stacked contributions)",
       title = sprintf("H by the current frequency of the alleles, %s", title_sfx)) +
  base_theme
ggsave(out("H_by_freq_class"), p, width = 13, height = 8)

# ---- 3. H = count x mean 2pq ----
cm <- seg[, .(n = .N, mean_2pq = mean(2 * frequencies * (1 - frequencies))), by = .(reg, psi_f, ID, gen, origin)]
cm_all <- seg[, .(n = .N, mean_2pq = mean(2 * frequencies * (1 - frequencies)), origin = "all"),
              by = .(reg, psi_f, ID, gen)]
cm <- rbind(cm, cm_all)[, .(n = mean(n), mean_2pq = mean(mean_2pq)), by = .(reg, psi_f, gen, origin)]
cm <- melt(cm, measure.vars = c("n", "mean_2pq"), variable.name = "what")
cm[, what := factor(what, levels = c("n", "mean_2pq"),
                    labels = c("number of segregating alleles", "mean 2pq per segregating allele"))]
cm[, origin := factor(origin, levels = c("all", "standing", "new"))]
pdf(out("H_count_x_mean"), width = 13, height = 8)
for (w in levels(cm$what)) {
  print(ggplot(cm[what == w], aes(gen, value, colour = origin)) +
          geom_line(linewidth = 0.6) +
          scale_colour_manual(NULL, values = c(all = "black", standing = "#b2182b", new = "#2166ac")) +
          labs(x = "tick", y = w, title = sprintf("H = count x mean 2pq / %d: %s, %s", maxL, w, title_sfx)) +
          base_theme)
}
invisible(dev.off())

# ---- 4. fixation rate by effect sign ----
fx <- a[isFixed == "T", .(first_fixed = min(gen)), by = .(reg, psi_f, ID, mut_id, effect_size)]
gens <- sort(unique(a$gen))
fx <- fx[first_fixed > min(gens)]               # fixed at the first recorded generation: unknown when
fx[, sign := ifelse(effect_size > 0, "positive", "negative")]
fr <- fx[, .N, by = .(reg, psi_f, ID, first_fixed, sign)]
grid <- CJ(ID = unique(a$ID), first_fixed = gens[-1], sign = c("positive", "negative"))
grid <- merge(grid, unique(a[, .(ID, reg, psi_f)]), by = "ID")
fr <- merge(grid, fr, by = c("ID", "first_fixed", "sign", "reg", "psi_f"), all.x = TRUE)
fr[is.na(N), N := 0]
fr[, dt := first_fixed - gens[match(first_fixed, gens) - 1]]
fr <- fr[, .(rate = mean(N / dt)), by = .(reg, psi_f, gen = first_fixed, sign)]
p <- ggplot(fr, aes(gen, rate, colour = sign)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual("effect", values = c(positive = "#b2182b", negative = "#2166ac")) +
  labs(x = "tick", y = "fixations per tick",
       title = sprintf("Fixation rate between recorded generations, %s", title_sfx)) +
  base_theme
ggsave(out("fixation_rate"), p, width = 13, height = 8)

# ---- 5. effect sizes of the segregating alleles ----
ef <- seg[, .(mean_abs = mean(abs(effect_size)), mean_eff = mean(effect_size)), by = .(reg, psi_f, ID, gen, origin)]
ef <- ef[, .(mean_abs = mean(mean_abs), mean_eff = mean(mean_eff)), by = .(reg, psi_f, gen, origin)]
ef <- melt(ef, measure.vars = c("mean_abs", "mean_eff"), variable.name = "what")
ef[, what := factor(what, labels = c("mean |effect|", "mean effect"))]
p <- ggplot(ef, aes(gen, value, colour = origin, linetype = what)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(NULL, values = c(standing = "#b2182b", new = "#2166ac")) +
  scale_linetype_manual(NULL, values = c("mean |effect|" = "solid", "mean effect" = "dashed")) +
  labs(x = "tick", y = "effect size of segregating alleles (unweighted mean)",
       title = sprintf("Effect sizes of the alleles still segregating, %s", title_sfx)) +
  base_theme
ggsave(out("effect_of_segregating"), p, width = 13, height = 8)

# numbers for the write-up: H and its parts at the shift and at the end
tab <- merge(s1[, .(H = sum(H)), by = .(reg, psi_f, gen)],
             dcast(s1, reg + psi_f + gen ~ group, value.var = "H"), by = c("reg", "psi_f", "gen"))
fwrite(tab[gen %in% c(shift_tick, max(gens))], out("H_parts_start_end", "csv"))
cat("plots written to", out_dir, "\n")
}
