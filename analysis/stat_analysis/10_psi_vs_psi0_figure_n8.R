# Section 10 compact figure for one group size: psi = x vs psi = 0, with the models (reg)
# as columns and one row per selection x metric (Stabilizing, Directional at px_shift, and
# each fluctuating setting). Drawn from the csv written by 10_psi_vs_psi0_fluc.R (nothing is
# recomputed). Run from the project root after that script has run once:
#   source("analysis/stat_analysis/10_psi_vs_psi0_figure_n8.R")
# Reads psi_vs_psi0/results/psi_vs_psi0_table_with_fluc{,_settings}.csv
# Writes psi_vs_psi0/paper/<name>.pdf for each entry of px_n8_figs.

source("analysis/stat_analysis/00_setup.R")

# ---- settings ----
px_n          <- 8      # group size shown
px_shift      <- 18.8   # directional shift size
# one figure per entry: file name -> fluctuating settings as (fluc_opt, interval) rows,
# in row order
px_n8_figs <- list(
  psi_vs_psi0_n8_main.pdf                = data.frame(fluc_opt = c(40, 40), fluc_interval = c(40, 80)),
  psi_vs_psi0_n8_opt_interval_pairs.pdf  = data.frame(fluc_opt = c(20, 80), fluc_interval = c(20, 80))
)
px_reach_keep   <- 50   # "Gens to X% of shift" thresholds shown (%)
# metrics left out (regular expressions): the 5-tick edge versions of peak / worst fitness
px_drop_metrics <- c("\\(last [0-9]+ ticks\\)", "\\(first [0-9]+ ticks\\)")
px_text_scale   <- 3    # all text x this (1 = ggplot default)
px_col_width    <- 2.2  # width of one model column (inches)
px_row_height   <- 0.9  # height of one metric row (inches)
px_wrap_width   <- 30   # metric labels wrap at this many characters
px_class_cols   <- c(worse = "#D7191C", mixed = "white", better = "#2B83BA")

# ---- read ----
px_csv_file <- px_out("results", "psi_vs_psi0_table_with_fluc.csv")
if (!file.exists(px_csv_file)) stop("Not found: ", px_csv_file,
                                    "\n  Run analysis/stat_analysis/10_psi_vs_psi0_fluc.R first.")
cells <- read.csv(px_csv_file, check.names = FALSE,
                  colClasses = c(psi_lower = "character", psi_dropped = "character"))
px_set <- read.csv(px_out("results", "psi_vs_psi0_table_with_fluc_settings.csv"), colClasses = "character")
px_set <- setNames(px_set$value, px_set$setting)

cells <- cells[cells$n == px_n, ]
if (!nrow(cells)) stop("No cells with n = ", px_n)
cells$selection  <- sub(" .*", "", cells$panel)
cells$shift_size <- NA_real_
is_dir <- cells$selection == "Directional"
cells$shift_size[is_dir] <- as.numeric(sub(".*shift = ", "", cells$panel[is_dir]))
cells$metric_order <- match(cells$metric, unique(cells$metric))   # csv order

reach_pct <- suppressWarnings(as.numeric(sub("^Gens to ([0-9.]+)% of shift.*", "\\1", cells$metric)))
cells <- cells[is.na(reach_pct) | reach_pct %in% px_reach_keep, ]
if (length(px_drop_metrics))
  cells <- cells[!Reduce(`|`, lapply(px_drop_metrics, grepl, x = cells$metric)), ]

cells$label <- paste0(cells$n_better, "/", cells$n_x)
cells$frac_better <- ifelse(cells$n_x > 0, cells$n_better / cells$n_x, NA)
wrap <- function(x) vapply(x, function(s) paste(strwrap(s, px_wrap_width), collapse = "\n"), "")

caption <- paste(
  "Better = higher fitness, fewer generations to reach the % of the shift, or lower phenotype lag",
  "behind the optimum (sum over the first generations after the shift of the fraction of the",
  "shift still to go). Runs not reaching the shift fraction by generation",
  px_set[["reach_deadline"]], "count as slowest. Fluctuating: fitness over ticks",
  paste0(px_set[["fluc_window"]], "; peak / worst = highest / lowest fitness of each half cycle,"),
  "averaged over half cycles. Label: number of psi = x values better than psi = 0.")

for (fname in names(px_n8_figs)) {
  pairs <- px_n8_figs[[fname]]
  st <- cells[cells$selection == "Stabilizing", ]
  di <- cells[cells$selection == "Directional" & abs(cells$shift_size - px_shift) < 1e-6, ]
  if (!nrow(di)) stop("No directional data with shift_size = ", px_shift, " and n = ", px_n)
  st$panel_lab <- "Stabilizing"; st$panel_order <- 0
  di$panel_lab <- paste0("Directional\nshift = ", px_shift); di$panel_order <- 1
  d <- rbind(st, di)
  for (i in seq_len(nrow(pairs))) {
    fl <- cells[cells$selection == "Fluctuating" & abs(cells$fluc_opt - pairs$fluc_opt[i]) < 1e-6 &
                  cells$fluc_interval == pairs$fluc_interval[i], ]
    if (!nrow(fl)) stop("No fluctuating data with fluc_opt = ", pairs$fluc_opt[i],
                        ", interval = ", pairs$fluc_interval[i], " and n = ", px_n)
    fl$panel_lab <- paste0("Fluctuating\nopt = ", fl$fluc_opt, ",\ninterval = ", fl$fluc_interval)
    fl$panel_order <- 2 + i / 1e3
    d <- rbind(d, fl)
  }
  d <- d[order(d$panel_order, d$metric_order), ]
  d$panel  <- factor(d$panel_lab, levels = unique(d$panel_lab))
  d$metric <- factor(wrap(d$metric), levels = rev(unique(wrap(d$metric))))   # first metric on top
  d$model  <- factor(d$model, levels = models)

  k <- px_text_scale
  n_rows <- nrow(unique(d[c("panel", "metric")]))
  width  <- 9 * k + px_col_width * nlevels(droplevels(d$model))
  height <- 4.5 * k + px_row_height * n_rows
  cap <- paste(strwrap(caption, floor(width * 72 / (0.5 * 0.8 * 11 * k))), collapse = "\n")
  p <- ggplot(d, aes(x = model, y = metric, fill = frac_better)) +
    geom_tile(color = "white", linewidth = 0.8 * k) +
    geom_text(aes(label = label), size = 3 * k) +
    facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(colors = px_class_cols, limits = c(0, 1), na.value = "grey80",
                         breaks = c(0, 0.5, 1),
                         labels = c("psi = 0 better\nfor every x", "half of x\nbetter",
                                    "psi = x better\nfor every x")) +
    labs(x = NULL, y = NULL, fill = "Fraction of psi = x better than psi = 0",
         title = paste0("psi = x vs psi = 0, n = ", px_n), caption = cap) +
    guides(fill = guide_colorbar(barwidth = 15 * k, barheight = 1.2 * k,
                                 title.position = "top", title.hjust = 0.5)) +
    theme_bw(base_size = 11 * k) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
  file <- px_out("paper", fname)
  ggsave(file, p, width = width, height = height, limitsize = FALSE)
  message("Saved ", file)
}
