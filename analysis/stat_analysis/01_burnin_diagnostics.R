# 1. burn-in diagnostic plots
# Section 1 of analysis/stat_analysis.Rmd. Run from the project root:
#   source("analysis/stat_analysis/01_burnin_diagnostics.R")
# Parameters and shared helper functions are in 00_setup.R.
source("analysis/stat_analysis/00_setup.R")
source("analysis/stat_analysis/00_load_data.R")

# One line per run, colored by model; facets separate psi (columns) and n (rows)
# so that runs from different parameter sets are not overlaid.
plot_burnin <- function(df, yvar, ylab, file) {
  p <- ggplot(df, aes(x = tick, y = .data[[yvar]], group = run, color = model)) +
    geom_line(alpha = 0.4, linewidth = 0.3) +
    facet_grid(n ~ psi, labeller = label_both) +
    scale_color_manual(values = model_cols) +
    labs(x = "Generation", y = ylab, color = "Model") +
    theme_bw()
  ggsave(file, p, width = 10, height = 6, dpi = 300)
  p
}

plot_burnin(burn_df, sd_col,   paste("Population", trait_word, "SD"),
            file.path(base_plot_folder, paste0("burnin_", sd_col, ".pdf")))
plot_burnin(burn_df, mean_col, paste("Population", trait_word, "mean"),
            file.path(base_plot_folder, paste0("burnin_", mean_col, ".pdf")))
