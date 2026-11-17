#!/usr/bin/env Rscript
# Worked example of the measures in allele_periodic_summary.R, on real data from one sim
# (fluc_opt ex_opt, fluc_interval ex_int, reg ex_reg, psi ex_psi, rep ex_rep).
# For a trajectory x(t) recorded at N times spanning a whole number of periods, the fitted wave is
#   x_fit(t) = mean(x) + |c| cos(omega t + arg c),  c = (2 / N) sum_t x(t) exp(-i omega t)
# i.e. the least-squares sine wave with the period of the optimum. |c| is the wave's height
# (swing), and the share of x's variance that the wave explains is (|c|^2 / 2) / var(x).
# Panels:
#   A  the optimum (square wave) and its own fitted wave
#   B  a rhythmic allele: the allele with the largest share of variance explained by the wave
#   C  a drifting allele: the allele with variance closest to B's but the smallest share explained
#   D  genetic mean g = sum 2 a p, with its fitted wave; the optimum's wave scaled for comparison;
#      gain = height of g's wave / height of the optimum's wave, lag = shift between the two waves
#
# Output: analysis/fluctuating_env/server/plots/allele_periodic/allele_fit_example.pdf
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/fluctuating_env/server/plot_allele_fit_example.R

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(ggplot2)
  library(patchwork)
})

# ---- settings ----
root    <- "/fs/cbsubscb10/storage/zc524/fluc_env_adjusted_shift_size/parquet/allele"
out_f   <- "analysis/fluctuating_env/server/plots/allele_periodic/allele_fit_example.pdf"
ex_opt  <- "103.4"
ex_int  <- 110L
ex_reg  <- "ave"
ex_psi  <- "0.0"
ex_rep  <- 1L
first_rec <- 10005
dt      <- 10
max_rec <- 100
min_var <- 0.005        # candidate alleles for B and C: temporal variance of p at least this
pt_col  <- "grey20"
fit_col <- "#E0402B"
opt_col <- "#2B4BE0"

fo <- as.numeric(ex_opt)
period <- 2 * ex_int
N <- max(which((seq_len(max_rec) * dt) %% period == 0))
times <- first_rec + (seq_len(N) - 1) * dt
last_rec <- max(times)
om <- 2 * pi / period
opt_at <- function(t) { k <- floor((t - 50) / ex_int) + 1; ifelse(k %% 2 == 1, fo / 2, -fo / 2) }

path <- file.path(root, sprintf("fluc_opt=%s/fluc_interval=%d/reg=%s/n=8/part-psi_%s-reps_1-0.parquet",
                                ex_opt, ex_int, ex_reg, ex_psi))
a <- as.data.table(open_dataset(path) |>
  dplyr::filter(rep == ex_rep, gen >= first_rec, gen <= last_rec) |>
  dplyr::select(gen, mut_id, effect_size, frequencies) |> dplyr::collect())

# full N x allele trajectories (absent = 0)
full <- CJ(mut_id = unique(a$mut_id), gen = times)
full <- a[full, on = .(mut_id, gen)]
full[is.na(frequencies), frequencies := 0]
full[, effect_size := effect_size[!is.na(effect_size)][1], by = mut_id]

fourier <- function(x, t) complex(real = 2 / N * sum(x * cos(om * t)), imaginary = -2 / N * sum(x * sin(om * t)))
wave <- function(x, t, tt) { c <- fourier(x, t); mean(x) + Mod(c) * cos(om * tt + Arg(c)) }

al <- full[, { c <- fourier(frequencies, gen); v <- mean(frequencies^2) - mean(frequencies)^2
               .(a = effect_size[1], amp = Mod(c), v = v, share = (Mod(c)^2 / 2) / v) }, by = mut_id]
cand <- al[v >= min_var]
B <- cand[which.max(share)]
C <- cand[a * B$a > 0][order(abs(v - B$v) + share * 10)][1]   # same effect sign, similar variance, low share
cat(sprintf("B: mut %d, a = %.2f, swing %.3f, share %.2f\nC: mut %d, a = %.2f, swing %.3f, share %.2f\n",
            B$mut_id, B$a, B$amp, B$share, C$mut_id, C$a, C$amp, C$share))

tt <- seq(min(times), max(times), by = 1)
flips <- seq(50, last_rec, by = ex_int); flips <- flips[flips > min(times)]
base <- list(theme_bw(base_size = 8),
             theme(panel.grid = element_blank(), plot.title = element_text(size = 8, face = "bold"),
                   plot.subtitle = element_text(size = 7)),
             geom_vline(xintercept = flips, colour = "grey85", linewidth = 0.3),
             scale_x_continuous(limits = range(times), expand = expansion(0.01)))

# A: optimum
o <- opt_at(times)
pA <- ggplot() + base +
  geom_step(aes(tt, opt_at(tt)), colour = opt_col, linewidth = 0.5) +
  geom_line(aes(tt, wave(o, times, tt)), colour = opt_col, linetype = "dashed", linewidth = 0.4) +
  labs(x = NULL, y = "optimum", title = "A  Optimum",
       subtitle = sprintf("flips every %d generations (period %d); dashed = its fitted wave", ex_int, period))

allele_panel <- function(r, ttl) {
  x <- full[mut_id == r$mut_id][order(gen)]
  fit <- wave(x$frequencies, x$gen, tt)
  m <- mean(x$frequencies)
  ggplot() + base +
    geom_hline(yintercept = m, colour = "grey60", linetype = "dotted") +
    geom_line(aes(tt, fit), colour = fit_col, linewidth = 0.6) +
    geom_point(aes(x$gen, x$frequencies), colour = pt_col, size = 0.7) +
    annotate("segment", x = tt[which.max(fit)], xend = tt[which.max(fit)], y = m, yend = max(fit),
             colour = fit_col, arrow = arrow(length = unit(3, "pt"), ends = "both")) +
    annotate("text", x = tt[which.max(fit)], y = max(fit), hjust = -0.1, vjust = 1, size = 2.4,
             colour = fit_col, label = sprintf("swing %.3f", r$amp)) +
    labs(x = NULL, y = "allele frequency", title = ttl,
         subtitle = sprintf("effect size %+.2f (%s)\npoints = recorded frequency, red = fitted wave; the wave explains %.0f%% of its variance",
                            r$a, if (r$a < 0) "lowers the trait, so favoured when the optimum is low"
                                 else "raises the trait, so favoured when the optimum is high",
                            100 * r$share))
}
pB <- allele_panel(B, "B  Rhythmic allele: frequency follows the optimum")
pC <- allele_panel(C, "C  Drifting allele: moves as much, but not in rhythm")

# D: genetic mean (fixed substitutions only add a constant, ignored)
g <- full[, .(g = sum(2 * effect_size * frequencies)), by = gen][order(gen)]
cg <- fourier(g$g, g$gen); co <- fourier(o, times)
gain <- Mod(cg) / Mod(co)
lag <- ((Arg(co) - Arg(cg)) %% (2 * pi)) / om
gfit <- wave(g$g, g$gen, tt)
ofit_scaled <- mean(g$g) + Mod(co) * cos(om * tt + Arg(co))     # optimum's wave, centred on g's mean
mid <- tt > quantile(tt, 0.35)                                      # a peak away from the edges
peak_o <- tt[mid][which.max(ofit_scaled[mid])]; peak_g <- peak_o + lag
pD <- ggplot() + base +
  geom_line(aes(tt, ofit_scaled), colour = opt_col, linetype = "dashed", linewidth = 0.4) +
  geom_line(aes(tt, gfit), colour = fit_col, linewidth = 0.6) +
  geom_point(aes(g$gen, g$g), colour = pt_col, size = 0.7) +
  annotate("segment", x = peak_o, xend = peak_g, y = max(ofit_scaled), yend = max(ofit_scaled),
           colour = "black", arrow = arrow(length = unit(3, "pt"))) +
  annotate("text", x = (peak_o + peak_g) / 2, y = max(ofit_scaled), vjust = -0.6, size = 2.4,
           label = sprintf("lag %.0f gen", lag)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "generation", y = "genetic mean (sum 2ap)",
       title = "D  Genetic mean: all alleles combined",
       subtitle = sprintf("red = fitted wave of g; blue dashed = optimum's wave; gain = red height / blue height = %.2f",
                          gain))

p <- pA / pB / pC / pD + plot_layout(heights = c(0.6, 1, 1, 1))
ggsave(out_f, p, width = 7, height = 8)
cat(sprintf("gain %.3f, lag %.1f generations; written %s\n", gain, lag, out_f))
