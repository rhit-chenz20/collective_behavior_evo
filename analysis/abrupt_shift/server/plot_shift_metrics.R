#!/usr/bin/env Rscript
# Summary metrics of the response to an abrupt shift of the optimum
# (cons_env_data_with_adjusted_shift_size; optimum jumps from 0 to shift_size at tick 10000).
# Columns = copying rule. For every shift size in shift_sizes it writes:
#
#   cum_load_n8_sz<s>     row 1: cumulative fitness load since the shift,
#                                 sum_{t >= 0} (1 - w_bar(t) / w_max), w_max = dnorm(0, 0, sqrt(Vs))
#                          row 2: cumulative excess load = same sum after subtracting each sim's
#                                 pre-shift load (mean over ticks 9950-9999): the extra cost of the shift
#                          n = 8, psi_show; dots = mean over reps (every load_step ticks), ribbon = 95% CI
#   t50_n2_8_20_50_sz<s>   generations until the population mean first reaches half the new optimum
#                          (>= shift_size / 2), row 1 phenotype z_bar, row 2 genotype g_bar; x = psi,
#                          colour = n; mean +- 95% CI over reps (log y). Sims that never reach it by
#                          tick 11000 count as 1000 (a lower bound); points with any such sim are hollow.
#   variance_n8_sz<s>      row 1: phenotypic variance V_z = z_sd^2, row 2: genetic variance V_g = g_sd^2
#                          over time; n = 8, psi_show; dots = mean over reps, ribbon = 95% CI
#
# Input (analysis/load_data/server/load_data.R -d cons_env_data_with_adjusted_shift_size -t data
#        --n 2,8,20,50 --cols z_mean,g_mean,z_sd,g_sd,fitness_mean --tag n2_8_20_50_allpsi_gen9950-11000):
#   figure/r_data/cons_env_data_with_adjusted_shift_size/data_summary_n2_8_20_50_allpsi_gen9950-11000.rds
#   (all psi, all shift sizes, reps 1-60, every tick 9950-11000)
# Output: analysis/abrupt_shift/server/plots/shift_metrics/
#
# Run from the repo root:
#   /programs/R-4.5.2/bin/Rscript analysis/abrupt_shift/server/plot_shift_metrics.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(MetBrewer)
  library(patchwork)
})

# ---- settings ----
data_dir    <- "figure/r_data/cons_env_data_with_adjusted_shift_size"
data_tag    <- "n2_8_20_50_allpsi_gen9950-11000"
out_dir     <- "analysis/abrupt_shift/server/plots/shift_metrics"
shift_sizes <- c(9.4, 18.8, 28.2)
shift_tick  <- 10000
end_tick    <- 11000
Vs          <- 2000           # width of the Gaussian fitness function (metadata Vs)
pre_window  <- c(9950, 9999)  # baseline load for the excess-load row
n_main      <- 8              # load and variance figures
t50_n       <- c(2, 8, 20, 50)
psi_show    <- c(0, 0.2, 0.5, 0.8)
psi_grid    <- seq(0, 1, by = 0.1)
n_tested    <- c(2, 4, 5, 8, 10, 20, 25, 40, 50, 100)
x_lim_load  <- c(0, 500)
load_step   <- 5              # plot every 5th tick of the (smooth) cumulative load
x_lim_var   <- c(-50, 300)
x_lim_pgl   <- c(-50, 300)    # load per generation (also its start in the 3-row load figure)
fig_h_1row  <- 1.1            # inches, one row of 16:11 panels with strip and x title
fig_h_3row  <- 2.62
# ---- figure size and style (same as plot_shift_response.R / plot_equilibrium_dist.R) ----
fig_w     <- 8.7 / 2.54
tick_pt   <- 5
txt_pt    <- 5
strip_pt  <- 5.5
sym_pt    <- 6.5
gap       <- unit(4, "pt")
panel_ratio <- 11 / 16
mean_lwd  <- 0.55
point_size <- 0.7
ts_point  <- 0.35             # dots of dense time series
ref_lwd   <- 0.3
err_lwd   <- 0.25
fig_h_2row <- 1.86            # inches, two rows of 16:11 panels at 8.7 cm
fig_theme <- theme(plot.margin = margin(2, 4, 2, 0), legend.background = element_blank(),
                   legend.key = element_blank(), legend.margin = margin(0, 1, 1, 0),
                   legend.text = element_text(size = txt_pt - 0.5, margin = margin(t = 0.5)),
                   axis.text = element_text(size = tick_pt),
                   axis.title.y = element_text(margin = margin(r = 1)),
                   axis.ticks.length = unit(1.2, "pt"),
                   aspect.ratio = panel_ratio,
                   panel.spacing.x = gap,
                   strip.text = element_text(size = strip_pt, margin = margin(1, 0, 1.5, 0)),
                   strip.background = element_blank())
reg_names <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")
reg_factor <- function(r) factor(reg_names[as.character(r)], levels = reg_names)
base_theme <- list(theme_bw(base_size = txt_pt), fig_theme, theme(panel.grid = element_blank()))
sym_title  <- theme(axis.title.y = element_text(size = sym_pt, angle = 0, vjust = 0.5,
                                                margin = margin(r = 1)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cap_pt    <- 4.5             # methods caption under each figure
# wrap a paragraph to the figure width, measuring the rendered text (systemfonts, in points)
wrap_measured <- function(txt, width_pt = fig_w * 72 - 4) {
  words <- strsplit(txt, " ", fixed = TRUE)[[1]]
  lines <- character(); cur <- ""
  for (w in words) {
    cand <- if (nzchar(cur)) paste(cur, w) else w
    if (nzchar(cur) && systemfonts::string_width(cand, family = "Helvetica", size = cap_pt, res = 72) > width_pt) {
      lines <- c(lines, cur); cur <- w
    } else cur <- cand
  }
  c(lines, cur)
}
# caption: one string per paragraph; the figure grows by the caption's height so the panels keep
# their size. cairo_pdf so the caption can use unicode (psi, bars, minus signs)
save_fig <- function(p, stem, h, caption) {
  lines <- unlist(lapply(caption, wrap_measured))
  cap_h <- (length(lines) * cap_pt * 1.18 + 5) / 72
  f <- file.path(out_dir, sprintf("%s_8.7cm.pdf", stem))
  p <- p + plot_annotation(caption = paste(lines, collapse = "\n"),
                           theme = theme(plot.margin = margin(0, 0, 0, 0),
                                         plot.caption = element_text(size = cap_pt, hjust = 0,
                                                                     lineheight = 1.05,
                                                                     margin = margin(t = 4)),
                                         plot.caption.position = "plot"))
  ggsave(f, p, width = fig_w, height = h + cap_h, device = cairo_pdf)
  cat("written", f, "\n")
}

psi_pal <- as.character(met.brewer("Hokusai2", length(psi_grid)))
psi_scale <- scale_colour_gradientn(
  expression(psi), colours = psi_pal, values = psi_grid, limits = range(psi_grid),
  breaks = psi_show, labels = format(psi_show, drop0trailing = TRUE),
  aesthetics = c("colour", "fill"),
  guide = guide_colourbar(nbin = length(psi_grid), display = "rectangles",
                          direction = "horizontal", title.position = "left",
                          theme = theme(legend.title = element_text(size = sym_pt, vjust = 0.75,
                                                                    margin = margin(r = 2))),
                          barwidth = unit(1.25, "cm"), barheight = unit(0.11, "cm")))
n_pal <- rev(as.character(met.brewer("Greek", length(n_tested))))
n_scale <- scale_colour_gradientn(
  "n", colours = n_pal, values = scales::rescale(seq_along(n_tested)),
  limits = c(1, length(n_tested)), breaks = match(t50_n, n_tested), labels = t50_n,
  guide = guide_colourbar(nbin = length(n_tested), display = "rectangles",
                          direction = "horizontal", title.position = "left",
                          theme = theme(legend.title = element_text(size = sym_pt, vjust = 0.75,
                                                                    margin = margin(r = 2))),
                          barwidth = unit(1.0, "cm"), barheight = unit(0.11, "cm")))
legend_inside <- function(x, y, jx = x, jy = y)
  theme(legend.position = "inside", legend.position.inside = c(x, y),
        legend.justification.inside = c(jx, jy))
# mean +- 95% CI over reps, as dots + ribbon
ts_layers <- function() list(
  geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), colour = NA, alpha = 0.3),
  geom_point(size = ts_point, shape = 16),
  facet_wrap(~ reg, nrow = 1))
summ <- function(d, v, by) d[, .(mean = mean(get(v)), se = sd(get(v)) / sqrt(.N), reps = .N), by = by]

# ---- data ----
d_all <- as.data.table(readRDS(file.path(data_dir, sprintf("data_summary_%s.rds", data_tag))))
w_max <- dnorm(0, 0, sqrt(Vs))
d_all[, load := 1 - fitness_mean / w_max]
setkey(d_all, shift_size, reg, n, psi, rep, tick)

for (sz in shift_sizes) {
  d <- d_all[abs(shift_size - sz) < 1e-9]
  szs <- format(sz)

  # ---------------- cumulative load (n = 8) ----------------
  l <- d[n == n_main & psi %in% psi_show]
  l[, base := mean(load[tick >= pre_window[1] & tick <= pre_window[2]]), by = .(reg, psi, rep)]
  l <- l[tick >= shift_tick]
  l[, `:=`(x = tick - shift_tick, cum = cumsum(load), cum_ex = cumsum(load - base)),
    by = .(reg, psi, rep)]
  l <- l[x <= x_lim_load[2] & x %% load_step == 0]
  lt <- summ(l, "cum", c("reg", "psi", "x")); le <- summ(l, "cum_ex", c("reg", "psi", "x"))
  lt[, reg := reg_factor(reg)]; le[, reg := reg_factor(reg)]
  xs_load <- scale_x_continuous(limits = x_lim_load, expand = expansion(0.02),
                                breaks = c(0, 250, 500), labels = c("0", "250", ""))
  p_lt <- ggplot(lt[order(psi)], aes(x, mean, colour = psi, fill = psi, group = psi)) +
    ts_layers() + psi_scale + xs_load +
    scale_y_continuous(expand = expansion(c(0.02, 0.05))) +
    labs(x = NULL, y = "Cumulative\nload") + base_theme +
    theme(axis.text.x = element_blank()) + legend_inside(0.03, 0.97, 0, 1)
  p_le <- ggplot(le[order(psi)], aes(x, mean, colour = psi, fill = psi, group = psi)) +
    ts_layers() + psi_scale + xs_load +
    scale_y_continuous(expand = expansion(c(0.02, 0.05))) +
    labs(x = "Generation after optimum shift", y = "Cumulative\nexcess load") + base_theme +
    theme(strip.text = element_blank(), legend.position = "none")
  save_fig(p_lt / p_le, sprintf("cum_load_n%d_sz%s", n_main, szs), fig_h_2row, c(
    sprintf(paste0("Load in a generation = 1 − W / W_max, where W is the population mean fitness and ",
                   "W_max = 1/√(2πV_s) (V_s = %g) is the fitness of an individual exactly at the optimum. ",
                   "Top: load summed over generations from the shift (generation 0) up to x. ",
                   "Bottom: the same sum after subtracting, from every generation, the simulation's mean load over ",
                   "the %d generations before the shift, i.e. the extra load caused by the shift."),
            Vs, diff(pre_window) + 1),
    sprintf(paste0("n = %d, optimum shift = %s; dots = mean over %d replicate simulations (every %dth generation), ",
                   "ribbons = 95%% CI (±1.96 SE)."), n_main, szs, max(lt$reps), load_step)))
  cat(sprintf("  load: %d reps per point; cumulative load at %d (psi 0 / 0.8):\n", max(lt$reps), x_lim_load[2]))
  print(dcast(lt[x == x_lim_load[2] & psi %in% c(0, 0.8)], reg ~ psi, value.var = "mean"))
  print(dcast(le[x == x_lim_load[2] & psi %in% c(0, 0.8)], reg ~ psi, value.var = "mean"))

  # ---------------- load per generation, before and after the shift (n = 8) ----------------
  lp <- d[n == n_main & psi %in% psi_show & tick - shift_tick >= x_lim_pgl[1] &
            tick - shift_tick <= x_lim_load[2]]
  lp[, x := tick - shift_tick]
  lp <- summ(lp, "load", c("reg", "psi", "x"))
  lp[, reg := reg_factor(reg)]
  vline <- geom_vline(xintercept = 0, colour = "grey60", linetype = "dotted", linewidth = ref_lwd)
  pgl_plot <- function(xl, brk, lab) {
    ggplot(lp[x <= xl[2]][order(psi)], aes(x, mean, colour = psi, fill = psi, group = psi)) +
      vline + ts_layers() + psi_scale +
      scale_x_continuous(limits = xl, expand = expansion(0.01), breaks = brk, labels = lab) +
      scale_y_continuous(expand = expansion(c(0.02, 0.05))) +
      labs(y = "Load per\ngeneration") + base_theme
  }
  pgl_caption <- c(
    sprintf(paste0("Load in a generation = 1 − W / W_max, where W is the population mean fitness and ",
                   "W_max = 1/√(2πV_s) (V_s = %g) is the fitness of an individual exactly at the optimum, ",
                   "i.e. the fraction of the best possible fitness the population misses in that generation. ",
                   "Before the shift it is the cost of phenotypic spread around the optimum (≈ V_z / 2V_s); ",
                   "the optimum jumps at generation 0 (dotted line)."), Vs),
    sprintf(paste0("n = %d, optimum shift = %s; dots = mean over %d replicate simulations (every generation), ",
                   "ribbons = 95%% CI (±1.96 SE)."), n_main, szs, max(lp$reps)))
  p_pgl <- pgl_plot(x_lim_var, c(0, 100, 200, 300), c("0", "", "200", "")) +
    labs(x = "Generation after optimum shift") + legend_inside(0.985, 0.97, 1, 1)
  save_fig(p_pgl, sprintf("load_per_gen_n%d_sz%s", n_main, szs), fig_h_1row, pgl_caption)

  # three rows: load per generation + the two cumulative rows, on one x axis from -50
  xl3 <- c(x_lim_pgl[1], x_lim_load[2])
  xs3 <- scale_x_continuous(limits = xl3, expand = expansion(0.01), breaks = c(0, 250, 500),
                            labels = c("0", "250", ""))
  p3a <- pgl_plot(xl3, c(0, 250, 500), c("0", "250", "")) + labs(x = NULL) +
    theme(axis.text.x = element_blank()) + legend_inside(0.985, 0.97, 1, 1)
  p3b <- p_lt + xs3 + vline + theme(strip.text = element_blank(), legend.position = "none")
  p3c <- p_le + xs3 + vline
  save_fig(p3a / p3b / p3c, sprintf("load_3row_n%d_sz%s", n_main, szs), fig_h_3row, c(
    pgl_caption[1],
    sprintf(paste0("Middle: load summed over generations from the shift up to x; its slope at any x is the load ",
                   "of that generation, and once adapted it settles at the pre-shift (spread) load. Bottom: the ",
                   "same sum after subtracting, from every generation, the simulation's mean load over the %d ",
                   "generations before the shift, i.e. the extra load caused by the shift."), diff(pre_window) + 1),
    sprintf(paste0("n = %d, optimum shift = %s; dots = mean over %d replicate simulations (top: every generation, ",
                   "middle and bottom: every %dth), ribbons = 95%% CI (±1.96 SE)."),
            n_main, szs, max(lp$reps), load_step)))

  # ---------------- time to half the new optimum (n = 2, 8, 20, 50) ----------------
  first_pass <- function(v, t) { i <- which(v >= sz / 2); if (length(i)) t[i[1]] - shift_tick else NA_real_ }
  t50 <- d[n %in% t50_n & tick >= shift_tick,
           .(z = first_pass(z_mean, tick), g = first_pass(g_mean, tick)), by = .(reg, n, psi, rep)]
  t50 <- melt(t50, id.vars = c("reg", "n", "psi", "rep"), variable.name = "what", value.name = "T50")
  t50[, censored := is.na(T50)]
  t50[censored == TRUE, T50 := end_tick - shift_tick]
  ts <- t50[, .(mean = mean(T50), se = sd(T50) / sqrt(.N), any_cens = any(censored),
                n_cens = sum(censored), reps = .N), by = .(what, reg, n, psi)]
  ts[, `:=`(reg = reg_factor(reg), n_idx = match(n, n_tested), n_f = factor(n, levels = t50_n),
            lo = pmax(mean - 1.96 * se, 0), hi = mean + 1.96 * se)]
  y_rng <- range(c(ts$lo, ts$hi))
  dodge <- position_dodge(width = 0.04)
  t50_plot <- function(w, ylab) {
    ggplot(ts[what == w], aes(psi, mean, colour = n_idx, group = n_f)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.03, linewidth = err_lwd, position = dodge) +
      geom_line(linewidth = mean_lwd * 0.6, position = dodge) +
      geom_point(aes(shape = any_cens), size = point_size, fill = "white", stroke = 0.3,
                 position = dodge) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
      facet_wrap(~ reg, nrow = 1) + n_scale +
      scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = expansion(0.06)) +
      # pseudo-log: log-like above ~3, linear near 0 (T50 = 0 when z_bar already sat above shift/2
      # at the shift, which happens at small shifts through pre-shift drift)
      scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10), limits = y_rng,
                         breaks = c(0, 3, 10, 30, 100, 300, 1000),
                         labels = c("0", "3", "10", "30", "100", "300", "1000")) +
      labs(x = NULL, y = ylab) + base_theme
  }
  p_tz <- t50_plot("z", expression(T[50] * "(" * bar(z) * ")")) +
    theme(axis.text.x = element_blank()) + legend_inside(0.03, 0.97, 0, 1)
  p_tg <- t50_plot("g", expression(T[50] * "(" * bar(g) * ")")) + labs(x = expression(psi)) +
    theme(strip.text = element_blank(), legend.position = "none",
          axis.title.x = element_text(size = sym_pt))
  save_fig(p_tz / p_tg, sprintf("t50_n%s_sz%s", paste(t50_n, collapse = "_"), szs), fig_h_2row, c(
    sprintf(paste0("T₅₀ = number of generations after the shift until the population mean phenotype (top) ",
                   "or mean genotypic value (bottom) first reaches half the new optimum (≥ %s), computed for each ",
                   "simulation (0 = already there in the generation of the shift). Points = mean over %d replicate ",
                   "simulations, bars = 95%% CI (±1.96 SE). Hollow points: at least one replicate never reached it ",
                   "by generation %d and was counted as %d, so the mean is a lower bound. Optimum shift = %s."),
            format(sz / 2), max(ts$reps), end_tick - shift_tick, end_tick - shift_tick, szs)))
  cat(sprintf("  T50: %d reps per point; sims not reaching shift/2 by tick %d: z %d, g %d\n",
              max(ts$reps), end_tick, t50[what == "z", sum(censored)], t50[what == "g", sum(censored)]))
  print(dcast(ts[n == n_main & psi %in% psi_show], what + psi ~ reg, value.var = "mean"), digits = 3)

  # ---------------- V_z and V_g over time (n = 8) ----------------
  v <- d[n == n_main & psi %in% psi_show & tick - shift_tick >= x_lim_var[1] & tick - shift_tick <= x_lim_var[2]]
  v[, `:=`(x = tick - shift_tick, Vz = z_sd^2, Vg = g_sd^2)]
  vz <- summ(v, "Vz", c("reg", "psi", "x")); vg <- summ(v, "Vg", c("reg", "psi", "x"))
  vz[, reg := reg_factor(reg)]; vg[, reg := reg_factor(reg)]
  v_lim <- c(0, max(vz$mean + 1.96 * vz$se, vg$mean + 1.96 * vg$se))
  xs_var <- scale_x_continuous(limits = x_lim_var, expand = expansion(0.01),
                               breaks = c(0, 100, 200, 300), labels = c("0", "", "200", ""))
  vline <- geom_vline(xintercept = 0, colour = "grey60", linetype = "dotted", linewidth = ref_lwd)
  p_vz <- ggplot(vz[order(psi)], aes(x, mean, colour = psi, fill = psi, group = psi)) +
    vline + ts_layers() + psi_scale + xs_var +
    scale_y_continuous(limits = v_lim, expand = expansion(c(0, 0.04))) +
    labs(x = NULL, y = expression(V[z])) + base_theme + sym_title +
    theme(axis.text.x = element_blank(), legend.position = "none")
  p_vg <- ggplot(vg[order(psi)], aes(x, mean, colour = psi, fill = psi, group = psi)) +
    vline + ts_layers() + psi_scale + xs_var +
    scale_y_continuous(limits = v_lim, expand = expansion(c(0, 0.04))) +
    labs(x = "Generation after optimum shift", y = expression(V[g])) + base_theme + sym_title +
    # psi colour bar in the empty top of the Copy-the-extreme V_g panel
    theme(strip.text = element_blank()) + legend_inside(0.985, 0.97, 1, 1)
  save_fig(p_vz / p_vg, sprintf("variance_n%d_sz%s", n_main, szs), fig_h_2row, c(
    sprintf(paste0("V_z = variance of individual phenotypes and V_g = variance of individual genotypic values ",
                   "(sum of allelic effects) across the %d individuals of the population, in each generation ",
                   "(squared z_sd and g_sd from the simulation output). Dots = mean over %d replicate simulations, ",
                   "ribbons = 95%% CI (±1.96 SE). n = %d, optimum shift = %s; dotted line = the shift."),
            1000, max(vz$reps), n_main, szs)))
  cat(sprintf("  variance: %d reps per point\n", max(vz$reps)))
}
