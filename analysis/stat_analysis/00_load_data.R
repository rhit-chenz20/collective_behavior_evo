# Load the constant-environment summary data and prepare the burn-in trajectories.
# Sourced by sections 1-9. The data is read only if cons_data_df is not already loaded.

# ---- load data (skipped if cons_data_df is already in the environment) ----
if (!exists("cons_data_df")) {
  cons_data_df <- read_rds_data(burnin_file)
}

# ---- prepare burn-in data: one trajectory per run ----
# Accept either "ticks" or "tick" as the generation column.
tick_col <- if ("ticks" %in% names(cons_data_df)) "ticks" else "tick"

burn_df <- cons_data_df[cons_data_df[[tick_col]] < burnin_end & cons_data_df$reg %in% models,
                        c(tick_col, "psi", "n", "reg", "rep", "shift_size", mean_col, sd_col)]
names(burn_df)[1] <- "tick"

# A run is one psi/n/reg/rep combination. Before the shift, shift_size does not matter,
# so keep a single shift_size per run; otherwise the same burn-in would be counted
# several times and inflate n.
run_key <- interaction(burn_df$psi, burn_df$n, burn_df$reg, burn_df$rep, drop = TRUE)
keep_sz <- ave(burn_df$shift_size, run_key, FUN = function(x) min(x, na.rm = TRUE))
burn_df <- burn_df[burn_df$shift_size == keep_sz, ]
burn_df <- burn_df[!duplicated(burn_df[c("psi", "n", "reg", "rep", "tick")]), ]

burn_df$model     <- factor(burn_df$reg, levels = models)
burn_df$run       <- paste0(burn_df$reg, "_psi", burn_df$psi, "_n", burn_df$n, "_rep", burn_df$rep)
burn_df$param_set <- paste0("psi = ", burn_df$psi, ", n = ", burn_df$n)
burn_df <- burn_df[order(burn_df$run, burn_df$tick), ]
