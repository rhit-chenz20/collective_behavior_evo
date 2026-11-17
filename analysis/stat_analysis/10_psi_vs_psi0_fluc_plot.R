# Plots for section 10 (psi = x vs psi = 0, constant + fluctuating environment), drawn
# from the csv written by 10_psi_vs_psi0_fluc.R, so nothing is recomputed. Run from the
# project root after that script has run once:
#   source("analysis/stat_analysis/10_psi_vs_psi0_fluc_plot.R")
# 10_psi_vs_psi0_fluc.R also sources this file at its end.
# Reads, from analysis/figure/<date>/psi_vs_psi0/results/ (date set in 00_setup.R):
#   psi_vs_psi0_table_with_fluc.csv           one row per panel x metric x model x n
#   psi_vs_psi0_table_with_fluc_settings.csv  settings used for it (for the captions)

source("analysis/stat_analysis/00_setup.R")

# ---- plot settings (only affect the figures) ----
# fluctuating settings shown in the combined table (all 30 would make it ~130 in wide);
# NULL = all. fluc_opt is 20-100 by 20, fluc_interval 20-120 by 20. Every setting is shown
# in the per-metric fluctuating overview plots.
px_fluc_opt_main      <- c(20, 60, 120)
px_fluc_interval_main <- c(20, 60, 120)
px_class_cols  <- c(better = "#2B83BA", worse = "#D7191C")
px_class_cols  <- c(px_class_cols, mixed = "white")   # middle of the gradient (half of x better)
px_wrap_width  <- 22                                  # metric labels wrap at this many characters


# ---- read ----
px_csv_file <- px_out("results", "psi_vs_psi0_table_with_fluc.csv")
if (!file.exists(px_csv_file)) stop("Not found: ", px_csv_file,
                                    "\n  Run analysis/stat_analysis/10_psi_vs_psi0_fluc.R first.")
# text columns that may be all empty would otherwise be read as logical NA
px_cell <- read.csv(px_csv_file, check.names = FALSE,
                    colClasses = c(psi_lower = "character", psi_dropped = "character"))
px_cell$psi_dropped[is.na(px_cell$psi_dropped)] <- ""
px_set <- read.csv(px_out("results", "psi_vs_psi0_table_with_fluc_settings.csv"), colClasses = "character")
px_set <- setNames(px_set$value, px_set$setting)
px_fluc_window    <- px_set[["fluc_window"]]
px_reach_deadline <- px_set[["reach_deadline"]]

# the csv has one-line labels; the rows are sorted by panel, then metric, so their order of
# first appearance is the plot order. Panels break after the first word
# ("Directional\nshift = 9.4"); metric labels are wrapped.
px_cell$panel  <- factor(sub(" ", "\n", px_cell$panel), levels = unique(sub(" ", "\n", px_cell$panel)))
wrap <- function(x) vapply(x, function(s) paste(strwrap(s, px_wrap_width), collapse = "\n"), "")
px_cell$metric <- factor(wrap(px_cell$metric), levels = unique(wrap(px_cell$metric)))
px_cell$model  <- factor(px_cell$model, levels = models)

for (k in c("opt", "interval")) {
  want <- get(paste0("px_fluc_", k, "_main"))
  have <- unique(na.omit(px_cell[[paste0("fluc_", k)]]))
  if (length(setdiff(want, have)))
    warning("px_fluc_", k, "_main: no data for ", paste(setdiff(want, have), collapse = ", "),
            " (values in the csv: ", paste(sort(have), collapse = ", "), ")", call. = FALSE)
}

# cell label: n_better / n_x, plus the psi values left out; grey if psi = 0 is incomplete
px_cell$label <- ifelse(!px_cell$psi0_complete, "psi = 0\nnot reached",
                        paste0(px_cell$n_better, "/", px_cell$n_x,
                               ifelse(px_cell$psi_dropped == "", "",
                                      paste0("\nexcl. psi = ", px_cell$psi_dropped))))
px_cell$frac_better <- ifelse(px_cell$n_x > 0, px_cell$n_better / px_cell$n_x, NA)

px_fill <- scale_fill_gradientn(colors = px_class_cols[c("worse", "mixed", "better")],
                                limits = c(0, 1), na.value = "grey80", breaks = c(0, 0.5, 1),
                                labels = c("psi = 0 better\nfor every x", "half of x\nbetter",
                                           "psi = x better\nfor every x"))
px_theme <- list(
  guides(fill = guide_colorbar(barwidth = 15, title.position = "top", title.hjust = 0.5)),
  theme_bw(),
  theme(panel.grid = element_blank(), legend.position = "bottom",
        strip.text.y = element_text(angle = 0)))

# table: rows = model x n, columns = selection regime (x shift size / fluctuation) x metric
# fill = fraction of psi = x better than psi = 0: 0 = worse for every x (red),
# 1 = better for every x (blue), mixed cells in between
main_cell <- px_cell[is.na(px_cell$fluc_opt) |
                       ((is.null(px_fluc_opt_main) | px_cell$fluc_opt %in% px_fluc_opt_main) &
                        (is.null(px_fluc_interval_main) | px_cell$fluc_interval %in% px_fluc_interval_main)), ]
p <- ggplot(main_cell, aes(x = metric, y = factor(n), fill = frac_better)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3, lineheight = 0.9) +
  facet_grid(model ~ panel, scales = "free", space = "free") +
  px_fill +
  labs(x = NULL, y = "n", fill = "Fraction of psi = x better than psi = 0",
       caption = paste("Better = higher fitness, fewer generations to reach the % of the shift, or lower phenotype lag",
                       "time (sum over the first T generations of the fraction of the shift still to go).",
                       "Fluctuating: half-cycle fitness over ticks", px_fluc_window,
                       if ("edge_k" %in% names(px_set)) paste0("(edge versions: mean of the first / last ",
                                                               px_set[["edge_k"]], " ticks of each half cycle)"),
                       ".",
                       "\nLabel: number of psi = x values (of those tested) better than psi = 0.",
                       "Runs not reaching the shift fraction by generation",
                       px_reach_deadline, "count as slowest (a psi = x with one is worse than psi = 0).")) +
  px_theme
# wide with many panels: limitsize = FALSE lets ggsave go past its 50 in limit
ggsave(px_out("overview", "psi_vs_psi0_table_with_fluc.pdf"), p,
       width  = 2.5 + 1.4 * nlevels(droplevels(interaction(main_cell$panel, main_cell$metric))),
       height = 2.5 + 0.5 * nrow(unique(main_cell[c("model", "n")])), limitsize = FALSE)

# fluctuating overview, one file per metric: columns = fluc_interval, x = fluc_opt
fl_cell <- px_cell[!is.na(px_cell$fluc_opt), ]
for (met in levels(droplevels(fl_cell$metric))) {
  d <- fl_cell[fl_cell$metric == met, ]
  q <- ggplot(d, aes(x = factor(fluc_opt), y = factor(n), fill = frac_better)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = 2.5, lineheight = 0.9) +
    facet_grid(model ~ fluc_interval, labeller = labeller(fluc_interval = function(x) paste("interval =", x))) +
    px_fill +
    labs(x = "fluc_opt (distance between the two optima)", y = "n",
         fill = "Fraction of psi = x better than psi = 0",
         title = paste("Fluctuating env:", gsub("\n", " ", met)),
         caption = paste("Ticks", px_fluc_window,
                         "| label: number of psi = x values better than psi = 0")) +
    px_theme
  ggsave(px_out("overview", paste0("psi_vs_psi0_fluc_", gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", met)), ".pdf")),
         q, width = 2.5 + 0.6 * length(unique(d$fluc_opt)) * length(unique(d$fluc_interval)),
         height = 2.5 + 0.4 * nrow(unique(d[c("model", "n")])), limitsize = FALSE)
}
message("Plots saved to ", file.path(px_folder, "overview"))
