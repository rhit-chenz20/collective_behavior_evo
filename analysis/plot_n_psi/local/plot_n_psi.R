# -------------------------------
# Simulation + plotting scaffold
# -------------------------------

set.seed(123)

# ---- Parameters you control ----
n_vals   <- c(2, 3, 4, 5, 8)                 # sweep over these n
psi_vals <- c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)  # sweep over these psi
n_reps   <- 1000                              # repetitions per (n, psi)
mus      <- c(-1000, -500, -100, -10, 0, 10, 100, 500, 1000) # Normal mean
sigma    <- sqrt(2000)                        # Normal sd

# Where to save plots
save_png   <- TRUE
output_dir <- "computed_plots"
if (save_png && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- YOUR computation goes here ----
# Must return a numeric vector of length n given the n sampled points and psi.
compute_fn <- function(n, x, psi) {
  num   <- (n - 1) * (1 - psi) * x + psi * sum(x)
  denom <- (n - 1 + psi) * (1 - psi)
  if (abs(denom) < .Machine$double.eps^0.5) return(rep(NA_real_, n))
  (num / denom)
}

# ---- Colors for psi ----
colors <- rainbow(length(psi_vals))
names(colors) <- as.character(psi_vals)

for (mu in mus) {
  # ---- Loop over n ----
  for (n in n_vals) {
    # Collect all computed values per psi
    all_vals_list <- vector("list", length(psi_vals))
    names(all_vals_list) <- as.character(psi_vals)
    
    for (p in seq_along(psi_vals)) {
      psi <- psi_vals[p]
      vals_all <- numeric(n_reps * n)
      idx <- 1L
      for (rep_i in seq_len(n_reps)) {
        x <- rnorm(n, mu, sigma)
        vals <- compute_fn(n, x, psi)     # length n
        vals_all[idx:(idx + n - 1L)] <- vals
        idx <- idx + n
      }
      all_vals_list[[p]] <- vals_all
    }
    
    # Determine common breaks across all psi for this n
    all_vals_combined <- unlist(all_vals_list, use.names = FALSE)
    num_bins <- 1000
    breaks <- seq(min(all_vals_combined), max(all_vals_combined), length.out = num_bins + 1)
    bin_mids <- 0.5 * (head(breaks, -1) + tail(breaks, -1))
    
    # Open PNG if saving
    if (save_png) {
      png(file.path(output_dir, sprintf("binned_lines_n_%d_mu_%d.png", n, mu)),
          width = 900, height = 650, res = 120)
    }
    
    # Expand right margin and allow drawing outside plot for legend
    op <- par(mar = c(5, 4, 4, 7), xpd = TRUE)
    
    # First pass: compute ylim from all psi (counts)
    ylim_max <- 0
    for (p in seq_along(psi_vals)) {
      h <- hist(all_vals_list[[p]], breaks = breaks, plot = FALSE)
      ylim_max <- max(ylim_max, h$counts)
    }
    
    # Empty plot with correct limits
    plot(NA, xlim = range(bin_mids), ylim = c(0, ylim_max),
         xlab = "Value bin midpoint",
         ylab = "Count",
         main = sprintf("Binned counts as lines (n = %d, mu = %d)", n, mu))
    
    # Add only lines for each psi
    for (p in seq_along(psi_vals)) {
      h <- hist(all_vals_list[[p]], breaks = breaks, plot = FALSE)
      lines(bin_mids, h$counts, lwd = 2, col = colors[p])
    }
    
    # Legend just outside the top-right of plotting area
    legend("topright", inset = c(-0.25, 0), bty = "n",
           legend = paste("psi =", psi_vals),
           col = colors, lty = 1, lwd = 2, cex = 0.9)
    
    # Reset par and close device if needed
    par(op)
    if (save_png) dev.off()
  }
}
