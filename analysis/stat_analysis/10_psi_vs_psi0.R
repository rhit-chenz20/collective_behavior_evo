# 10. psi = x vs psi = 0 across selection regimes
# Section 10 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/10_psi_vs_psi0.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
px_load_fluc <- FALSE   # constant env only
source("analysis/stat_analysis/00_load_data_px.R")

# For each model and n, and each metric, compare the mean over runs at every psi = x > 0
# with the mean at psi = 0. A cell is "higher" if psi = x is above psi = 0 for every x,
# "lower" if below for every x, and "mixed" otherwise. Colors and labels are by "better"
# (blue) vs "worse" (red): higher is better for fitness, lower is better for the
# generations to reach the shift (px_lower_better); the label gives how many x are better.
# Metrics, one value per run:
#   stabilizing: fitness_mean averaged over the burn-in (ticks < shift_time; recorded from 9950)
#   directional: fitness_mean at shift_time (generation 10000);
#                for each fraction f in px_reach_frac, generations after the shift until
#                z_mean / g_mean first reach optimum_burnin + f * shift_size (assumes an
#                upward shift). Runs that have not reached it by px_reach_deadline count as
#                slower than every run that did (Inf), so the mean of a psi with any such
#                run is Inf: a psi = x with one is worse than psi = 0 (not better when
#                psi = 0 has one too), and every psi = x without one is better than a
#                psi = 0 with one. These psi are listed in psi_unreached in the csv.
#                phenotype lag time for each T in px_lag_T: sum over generations
#                10000 .. 10000 + T - 1 of (optimum_burnin + shift_size - z_mean) / shift_size
#                (lower = faster; compared by the mean over runs).
# Fluctuating-environment metrics can be added later as more rows of px_run.

shift_time     <- 10000
optimum_burnin <- 0               # optimum before the shift (SLiM: optimum = 0.0 until shift_time)
px_psi         <- NULL            # psi = x values compared with psi = 0; NULL = all in the data
px_psi_exclude <- c()             # psi values dropped before anything is computed; c() = none
px_reach_frac  <- c(0.5, 0.9)     # fractions of shift_size for the "generations to reach" metrics
px_reach_deadline <- 11000        # every run must reach each fraction by this generation
px_lag_T       <- c(25, 100)      # windows (generations after the shift) for the phenotype lag time
px_class_cols  <- c(better = "#2B83BA", worse = "#D7191C")
px_class_cols  <- c(px_class_cols, mixed = "white")   # middle of the gradient (half of x better)

px_df <- px_data_df[px_data_df$reg %in% models,
                      intersect(c("tick", "psi", "n", "reg", "rep", "shift_size", "ID",
                                  "fitness_mean", "z_mean", "g_mean"), names(px_data_df))]
names(px_df)[1] <- "tick"
px_df$psi <- round(px_df$psi, 2)  # psi may be stored as float32 (0.2 -> 0.2000000029)
if (0 %in% px_psi_exclude) stop("psi = 0 is the reference and cannot be excluded")
if (length(px_psi_exclude)) px_df <- px_df[!px_df$psi %in% round(px_psi_exclude, 2), ]
# ID is the unique sim key in data loaded with load_data.R; older data has no ID
px_df$run <- if ("ID" %in% names(px_df)) px_df$ID else
  paste(px_df$shift_size, px_df$reg, px_df$n, px_df$psi, px_df$rep)
px_info <- px_df[!duplicated(px_df$run), c("run", "shift_size", "reg", "n", "psi", "rep")]

# per-run values (named by run) -> rows of px_run
px_rows <- function(values, selection, metric, by_shift) {
  info <- px_info[match(names(values), as.character(px_info$run)), ]
  data.frame(selection = selection, metric = metric,
             shift_size = if (by_shift) info$shift_size else NA,
             reg = info$reg, n = info$n, psi = info$psi, rep = info$rep,
             value = as.numeric(values))
}

# stabilizing: the burn-in is shared across shift sizes, so keep one shift_size per run
# (as in sections 2-9). fitness_mean is 0 (not recorded) before shift_time - 50.
st <- px_df[px_df$tick < shift_time & !is.na(px_df$fitness_mean) & px_df$fitness_mean > 0, ]
st_key <- interaction(st$psi, st$n, st$reg, st$rep, drop = TRUE)
st <- st[st$shift_size == ave(st$shift_size, st_key, FUN = function(x) min(x, na.rm = TRUE)), ]
px_list <- list(px_rows(tapply(st$fitness_mean, st$run, mean), "Stabilizing",
                        "Mean fitness\n(gen < 10000)", by_shift = FALSE))

# directional: fitness right at the shift
d0 <- px_df[px_df$tick == shift_time & !is.na(px_df$fitness_mean), ]
px_list[[2]] <- px_rows(tapply(d0$fitness_mean, d0$run, mean), "Directional",
                        "Fitness at\ngen 10000", by_shift = TRUE)

# directional: generations until the population mean reaches a fraction of the shift
post <- px_df[px_df$tick >= shift_time & px_df$tick <= px_reach_deadline, ]
runs <- as.character(unique(post$run))
for (f in px_reach_frac) for (col in c("z_mean", "g_mean")) {
  tr  <- c(z_mean = "phenotype", g_mean = "genotype")[[col]]
  hit <- !is.na(post[[col]]) & post[[col]] >= optimum_burnin + f * post$shift_size
  t_hit <- tapply(post$tick[hit], post$run[hit], min)
  t_reach <- setNames(rep(Inf, length(runs)), runs)       # Inf = not reached by the deadline
  t_reach[names(t_hit)] <- t_hit - shift_time
  px_list[[length(px_list) + 1]] <- px_rows(t_reach, "Directional",
    paste0("Gens to ", 100 * f, "% of shift\n(", tr, ")"), by_shift = TRUE)
}
# directional: phenotype lag time, the sum over the first T generations after the shift of
# the fraction of the shift still to go, (opt - z_mean) / shift_size with
# opt = optimum_burnin + shift_size. In generations; lower = faster. Unlike the crossing
# time it is not biased by generation-to-generation noise in z_mean. NA if a generation of
# the window is missing.
for (TT in px_lag_T) {
  w <- px_df[px_df$tick >= shift_time & px_df$tick < shift_time + TT & !is.na(px_df$z_mean), ]
  s <- rowsum(cbind(1, (optimum_burnin + w$shift_size - w$z_mean) / w$shift_size), w$run)
  lag <- setNames(ifelse(s[, 1] == TT, s[, 2], NA), rownames(s))
  px_list[[length(px_list) + 1]] <- px_rows(lag, "Directional",
    paste0("Phenotype lag behind optimum\n(first ", TT, " gens)"), by_shift = TRUE)
}
px_run <- do.call(rbind, px_list)
px_lower_better <- grep("^Gens to |^Phenotype lag", unique(px_run$metric), value = TRUE)  # fewer gens / less lag
rm(st, d0, post)

# panels: one for stabilizing, one per shift size for directional
px_run$panel <- ifelse(is.na(px_run$shift_size), px_run$selection,
                       paste0(px_run$selection, "\nshift = ", px_run$shift_size))
px_panels <- c("Stabilizing", paste0("Directional\nshift = ", sort(unique(px_info$shift_size))))
px_panels <- c(px_panels, setdiff(unique(px_run$panel), px_panels))   # e.g. fluctuating later

# compare each psi = x with psi = 0 within panel x metric x model x n
px_all_psi   <- sort(unique(px_run$psi))
px_cell_list <- list()
for (g in split(px_run, list(px_run$panel, px_run$metric, px_run$reg, px_run$n), drop = TRUE)) {
  # NA if any run of that psi has no value; Inf if any run did not reach the threshold
  m    <- tapply(g$value, g$psi, mean)
  psis <- as.numeric(names(m))
  unreached <- is.infinite(m)
  if (!any(psis == 0)) next
  m0 <- m[[which(psis == 0)]]
  tested <- psis != 0 & (if (is.null(px_psi)) TRUE else psis %in% px_psi)
  if (!any(tested)) next
  dropped <- tested & is.na(m)                # psi = x with unsatisfied runs: not compared
  xs <- if (is.na(m0)) tested & FALSE else tested & !is.na(m)
  above <- m[xs] > m0
  better <- if (g$metric[1] %in% px_lower_better) m[xs] < m0 else above
  classify <- function(v, yes, no) if (!length(v)) NA else if (all(v)) yes else if (!any(v)) no else "mixed"
  # mean over runs at psi = 0 and at each tested psi = x (NA if not tested or not all runs have a value)
  means <- setNames(rep(NA_real_, length(px_all_psi)), paste0("mean_psi_", px_all_psi))
  means[paste0("mean_psi_", psis[tested | psis == 0])] <- m[tested | psis == 0]
  px_cell_list[[length(px_cell_list) + 1]] <- data.frame(
    panel = g$panel[1], metric = g$metric[1], model = g$reg[1], n = g$n[1],
    n_x = sum(xs), n_higher = sum(above),
    class = classify(above, "higher", "lower"),
    psi_lower = paste(psis[xs][!above], collapse = ", "),   # x values with psi = 0 >= psi = x
    n_better = sum(better),
    color_class = classify(better, "better", "worse"),
    psi0_complete = !is.na(m0),
    psi_dropped = paste(psis[dropped], collapse = ", "),    # x values left out (missing values)
    psi_unreached = paste(psis[unreached], collapse = ", "), # psi with runs not reaching the threshold
    n_runs = nrow(g), n_runs_na = sum(is.na(g$value)), as.list(means), check.names = FALSE)
}
px_cell <- do.call(rbind, px_cell_list)
px_cell$panel  <- factor(px_cell$panel, levels = px_panels)
px_cell$metric <- factor(px_cell$metric, levels = unique(px_run$metric))
px_cell$model  <- factor(px_cell$model, levels = models)
px_cell$color_class <- factor(px_cell$color_class, levels = names(px_class_cols))
px_csv <- px_cell[order(px_cell$panel, px_cell$metric, px_cell$model, px_cell$n), ]
px_csv$panel  <- gsub("\n", " ", px_csv$panel)    # one line per row in the csv
px_csv$metric <- gsub("\n", " ", px_csv$metric)
write.csv(px_csv, px_out("constant_only", "psi_vs_psi0_table.csv"), row.names = FALSE)

# cell label: n_better / n_x, plus the psi values left out; grey if psi = 0 is incomplete
px_cell$label <- ifelse(!px_cell$psi0_complete, "psi = 0\nnot reached",
                        paste0(px_cell$n_better, "/", px_cell$n_x,
                               ifelse(px_cell$psi_dropped == "", "",
                                      paste0("\nexcl. psi = ", px_cell$psi_dropped))))
px_cell$frac_better <- ifelse(px_cell$n_x > 0, px_cell$n_better / px_cell$n_x, NA)

# table: rows = model x n, columns = selection regime (x shift size) x metric
# fill = fraction of psi = x better than psi = 0: 0 = worse for every x (red),
# 1 = better for every x (blue), mixed cells in between
p <- ggplot(px_cell, aes(x = metric, y = factor(n), fill = frac_better)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3, lineheight = 0.9) +
  facet_grid(model ~ panel, scales = "free", space = "free") +
  scale_fill_gradientn(colors = px_class_cols[c("worse", "mixed", "better")], limits = c(0, 1),
                       na.value = "grey80",
                       breaks = c(0, 0.5, 1),
                       labels = c("psi = 0 better\nfor every x", "half of x\nbetter",
                                  "psi = x better\nfor every x")) +
  labs(x = NULL, y = "n", fill = "Fraction of psi = x better than psi = 0",
       caption = paste("Better = higher mean fitness, fewer generations to reach the % of the shift, or lower phenotype",
                       "lag time (sum over the first T generations of the fraction of the shift still to go).",
                       "Label: number of psi = x values (of those tested) better than psi = 0.",
                       "\nRuns not reaching the shift fraction by generation",
                       px_reach_deadline, "count as slowest (a psi = x with one is worse than psi = 0).")) +
  guides(fill = guide_colorbar(barwidth = 15, title.position = "top", title.hjust = 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        strip.text.y = element_text(angle = 0))
ggsave(px_out("constant_only", "psi_vs_psi0_table.pdf"), p,
       width  = 2.5 + 1.4 * nlevels(droplevels(interaction(px_cell$panel, px_cell$metric))),
       height = 2.5 + 0.5 * nrow(unique(px_cell[c("model", "n")])))
p
