# Section 10 figures for the paper: psi = x vs psi = 0 under stabilizing, directional and
# fluctuating selection, drawn from the csv written by 10_psi_vs_psi0_fluc.R (nothing is
# recomputed). Run from the project root after that script has run once:
#   source("analysis/stat_analysis/10_psi_vs_psi0_figure.R")
# The directional shift size and the fluctuating amplitude are matched: a fluctuating
# setting is shown with the directional panel whose shift_size equals its fluc_opt (each
# flip moves the optimum by fluc_opt). Only n values present in both environments are
# shown.
#   main text:  Stabilizing | Directional shift = px_main_shift | Fluctuating opt =
#               px_main_shift (or px_main_fluc_opt), one panel per interval in
#               px_main_interval; metrics in px_main_metrics
#   supplement: one file per directional shift size: Stabilizing | Directional shift = S |
#               Fluctuating opt = S at every interval
# Metrics matching px_drop_metrics are left out of both.
# If no fluctuating data has fluc_opt equal to a shift size (e.g. before the matched runs
# are loaded), the nearest fluc_opt is used with a warning, and the panel label shows the
# value actually used.
# Reads, from analysis/figure/<date>/psi_vs_psi0/results/ (date set in 00_setup.R):
#   psi_vs_psi0_table_with_fluc.csv, psi_vs_psi0_table_with_fluc_settings.csv
# Writes to psi_vs_psi0/paper/: psi_vs_psi0_main.pdf, psi_vs_psi0_supp_shift_<S>.pdf, and one more main-style
# figure per entry of px_pair_figs (fluctuating panels given as fluc_opt x interval pairs)

source("analysis/stat_analysis/00_setup.R")

# ---- settings ----
px_main_shift    <- 18.8       # directional shift size in the main figure
# fluc_opt of the main figure's fluctuating panels: NULL = the one equal to px_main_shift
# (nearest if none is, with a warning), or a number to choose it, e.g. 40
px_main_fluc_opt <- 40
px_main_interval <- c(40, 80)  # fluctuating intervals in the main figure, one panel each
# extra main-style figures: file name -> fluctuating panels as (fluc_opt, interval) pairs, in
# panel order; Stabilizing and Directional at px_main_shift as in the main figure
px_pair_figs <- list(
  psi_vs_psi0_main_opt_interval_pairs.pdf = data.frame(fluc_opt = c(20, 80), fluc_interval = c(20, 80))
)
# metrics in the main figure (regular expressions on the one-line metric names in the csv);
# NULL = all. e.g. c("^Mean fitness", "^Fitness at", "genotype", "^Phenotype lag",
#                    "^Cycle average", "max", "min")
px_main_metrics  <- NULL
# metrics left out of every figure (regular expressions): the 5-tick edge versions of
# peak / worst fitness
px_drop_metrics  <- c("\\(last [0-9]+ ticks\\)", "\\(first [0-9]+ ticks\\)")
px_text_scale    <- 3          # all text x this (1 = ggplot default)
px_col_width     <- 2.6        # width of one metric column (inches)
px_row_height    <- 1.0        # height of one model x n row (inches)
px_common_n      <- TRUE   # only n values present in both environments (FALSE = all n)
px_reach_keep    <- 50     # "Gens to X% of shift" thresholds shown (%), main and supplement
px_class_cols    <- c(worse = "#D7191C", mixed = "white", better = "#2B83BA")
px_wrap_width    <- 14     # metric labels wrap at this many characters


# ---- read ----
px_csv_file <- px_out("results", "psi_vs_psi0_table_with_fluc.csv")
if (!file.exists(px_csv_file)) stop("Not found: ", px_csv_file,
                                    "\n  Run analysis/stat_analysis/10_psi_vs_psi0_fluc.R first.")
# text columns that may be all empty would otherwise be read as logical NA
cells <- read.csv(px_csv_file, check.names = FALSE,
                  colClasses = c(psi_lower = "character", psi_dropped = "character"))
cells$psi_dropped[is.na(cells$psi_dropped)] <- ""
px_set <- read.csv(px_out("results", "psi_vs_psi0_table_with_fluc_settings.csv"), colClasses = "character")
px_set <- setNames(px_set$value, px_set$setting)

cells$selection  <- sub(" .*", "", cells$panel)                 # Stabilizing / Directional / Fluctuating
cells$shift_size <- NA_real_
is_dir <- cells$selection == "Directional"
cells$shift_size[is_dir] <- as.numeric(sub(".*shift = ", "", cells$panel[is_dir]))
cells$metric_order <- match(cells$metric, unique(cells$metric))  # csv order = plot order

# keep only the chosen "Gens to X% of shift" thresholds
reach_pct <- suppressWarnings(as.numeric(sub("^Gens to ([0-9.]+)% of shift.*", "\\1", cells$metric)))
if (!any(reach_pct %in% px_reach_keep))
  warning("px_reach_keep: no 'Gens to ", paste(px_reach_keep, collapse = "/"), "% of shift' metric in the csv",
          call. = FALSE)
cells <- cells[is.na(reach_pct) | reach_pct %in% px_reach_keep, ]
if (length(px_drop_metrics))
  cells <- cells[!Reduce(`|`, lapply(px_drop_metrics, grepl, x = cells$metric)), ]

# the phenotype lag time is added by 10_psi_vs_psi0_fluc.R; an older csv has none
if (!any(grepl("^Phenotype lag", cells$metric)))
  warning("No 'Phenotype lag' metric in ", px_csv_file, ": it was written before the lag was ",
          "added. Rerun analysis/stat_analysis/10_psi_vs_psi0_fluc.R.", call. = FALSE)

# n present in both environments
n_cons <- unique(cells$n[cells$selection != "Fluctuating"])
n_fluc <- unique(cells$n[cells$selection == "Fluctuating"])
n_keep <- if (px_common_n) sort(intersect(n_cons, n_fluc)) else sort(union(n_cons, n_fluc))
if (!length(n_keep)) stop("No n value has data in both environments (constant: ",
                          paste(sort(n_cons), collapse = ", "), "; fluctuating: ",
                          paste(sort(n_fluc), collapse = ", "), ")")
message("n shown: ", paste(n_keep, collapse = ", "),
        if (px_common_n) paste0("  (constant has ", paste(sort(n_cons), collapse = ", "),
                                "; fluctuating has ", paste(sort(n_fluc), collapse = ", "), ")"))
cells <- cells[cells$n %in% n_keep, ]

# fluc_opt matching a shift size: exact if present, otherwise the nearest (with a warning)
fluc_opts <- sort(unique(cells$fluc_opt[cells$selection == "Fluctuating"]))
match_opt <- function(S) {
  hit <- fluc_opts[abs(fluc_opts - S) < 1e-6]
  if (length(hit)) return(hit[1])
  near <- fluc_opts[which.min(abs(fluc_opts - S))]
  warning("No fluctuating data with fluc_opt = ", S, "; using the nearest, fluc_opt = ", near,
          " (available: ", paste(fluc_opts, collapse = ", "), ")", call. = FALSE)
  near
}

# cell label and fill, as in the other section 10 figures
cells$label <- ifelse(!cells$psi0_complete, "psi = 0\nnot reached",
                      paste0(cells$n_better, "/", cells$n_x,
                             ifelse(cells$psi_dropped == "", "",
                                    paste0("\nexcl. psi = ", cells$psi_dropped))))
cells$frac_better <- ifelse(cells$n_x > 0, cells$n_better / cells$n_x, NA)
wrap <- function(x) vapply(x, function(s) paste(strwrap(s, px_wrap_width), collapse = "\n"), "")

caption <- paste(
  "Better = higher fitness, fewer generations to reach the % of the shift, or lower phenotype lag",
  "behind the optimum (sum over the first generations after the shift of the fraction of the",
  "shift still to go). Runs not reaching the shift fraction by generation",
  px_set[["reach_deadline"]], "count as slowest. Fluctuating: fitness over ticks",
  paste0(px_set[["fluc_window"]], "; peak / worst = highest / lowest fitness of each half cycle,"),
  "averaged over half cycles. Directional shift size = fluctuating fluc_opt (distance between",
  "the two optima). Label: number of psi = x values better than psi = 0.")

# one table: rows = model x n, columns = panel x metric
draw_table <- function(d, file, title) {
  d <- d[order(d$panel_order, d$metric_order), ]
  d$panel  <- factor(d$panel_lab, levels = unique(d$panel_lab))
  d$metric <- factor(wrap(d$metric), levels = unique(wrap(d$metric)))
  d$model  <- factor(d$model, levels = models)
  k <- px_text_scale
  width  <- 3 * k + px_col_width * nlevels(droplevels(interaction(d$panel, d$metric)))
  height <- 4 * k + px_row_height * nrow(unique(d[c("model", "n")]))
  # caption (0.8 x base text) wrapped to the figure width; a character is ~0.5 x its size
  cap <- paste(strwrap(caption, floor(width * 72 / (0.5 * 0.8 * 11 * k))), collapse = "\n")
  p <- ggplot(d, aes(x = metric, y = factor(n), fill = frac_better)) +
    geom_tile(color = "white", linewidth = 0.8 * k) +
    geom_text(aes(label = label), size = 3 * k, lineheight = 0.9) +
    facet_grid(model ~ panel, scales = "free", space = "free") +
    scale_fill_gradientn(colors = px_class_cols, limits = c(0, 1), na.value = "grey80",
                         breaks = c(0, 0.5, 1),
                         labels = c("psi = 0 better\nfor every x", "half of x\nbetter",
                                    "psi = x better\nfor every x")) +
    labs(x = NULL, y = "n", fill = "Fraction of psi = x better than psi = 0",
         title = title, caption = cap) +
    guides(fill = guide_colorbar(barwidth = 15 * k, barheight = 1.2 * k,
                                 title.position = "top", title.hjust = 0.5)) +
    theme_bw(base_size = 11 * k) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          strip.text.y = element_text(angle = 0))
  ggsave(file, p, width = width, height = height, limitsize = FALSE)
  message("Saved ", file)
}

# cells of: Stabilizing, Directional shift = S, Fluctuating opt = match_opt(S) at intervals
pick <- function(S, intervals, opt = NULL) {
  if (is.null(opt)) opt <- match_opt(S)
  st <- cells[cells$selection == "Stabilizing", ]
  di <- cells[cells$selection == "Directional" & abs(cells$shift_size - S) < 1e-6, ]
  fl <- cells[cells$selection == "Fluctuating" & abs(cells$fluc_opt - opt) < 1e-6 &
                (is.null(intervals) | cells$fluc_interval %in% intervals), ]
  if (!nrow(di)) stop("No directional data with shift_size = ", S)
  if (!nrow(fl)) stop("No fluctuating data with fluc_opt = ", opt,
                      if (!is.null(intervals)) paste0(" and interval ", paste(intervals, collapse = ", ")))
  st$panel_lab <- "Stabilizing";                       st$panel_order <- 0
  di$panel_lab <- paste0("Directional\nshift = ", S);   di$panel_order <- 1
  fl$panel_lab <- paste0("Fluctuating\nopt = ", fl$fluc_opt, ", interval = ", fl$fluc_interval)
  fl$panel_order <- 2 + fl$fluc_interval / 1e4
  rbind(st, di, fl)
}

# ---- main figure ----
main <- pick(px_main_shift, px_main_interval, px_main_fluc_opt)
if (!is.null(px_main_metrics)) {
  keep <- Reduce(`|`, lapply(px_main_metrics, grepl, x = main$metric))
  main <- main[keep, ]
  if (!nrow(main)) stop("px_main_metrics matches no metric")
}
draw_table(main, px_out("paper", "psi_vs_psi0_main.pdf"),
           paste0("psi = x vs psi = 0: shift size ", px_main_shift,
                  ", fluctuation interval ", paste(px_main_interval, collapse = " / ")))

# ---- extra main-style figures: fluctuating panels as (fluc_opt, interval) pairs ----
for (fname in names(px_pair_figs)) {
  pairs <- px_pair_figs[[fname]]
  d <- pick(px_main_shift, NULL, pairs$fluc_opt[1])        # Stabilizing + Directional
  d <- d[d$selection != "Fluctuating", ]
  for (i in seq_len(nrow(pairs))) {
    fl <- cells[cells$selection == "Fluctuating" & abs(cells$fluc_opt - pairs$fluc_opt[i]) < 1e-6 &
                  cells$fluc_interval == pairs$fluc_interval[i], ]
    if (!nrow(fl)) stop("No fluctuating data with fluc_opt = ", pairs$fluc_opt[i],
                        " and interval = ", pairs$fluc_interval[i])
    fl$panel_lab <- paste0("Fluctuating\nopt = ", fl$fluc_opt, ", interval = ", fl$fluc_interval)
    fl$panel_order <- 2 + i / 1e3
    d <- rbind(d, fl)
  }
  if (!is.null(px_main_metrics))
    d <- d[Reduce(`|`, lapply(px_main_metrics, grepl, x = d$metric)), ]
  draw_table(d, px_out("paper", fname),
             paste0("psi = x vs psi = 0: shift size ", px_main_shift, ", fluctuation (opt, interval) = ",
                    paste0("(", pairs$fluc_opt, ", ", pairs$fluc_interval, ")", collapse = " / ")))
}

# ---- supplement: every shift size, its matched fluctuation at every interval ----
for (S in sort(unique(na.omit(cells$shift_size)))) {
  draw_table(pick(S, NULL), px_out("paper", paste0("psi_vs_psi0_supp_shift_", S, ".pdf")),
             paste0("psi = x vs psi = 0: shift size ", S, ", every fluctuation interval"))
}
