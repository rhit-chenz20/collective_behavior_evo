# Shared loading for the steady_tracking summary analyses: source() it from a plot script.
#
# The steady_tracking dataset was filled by two runs, loaded into separate .rds files:
#   steady_tracking  speeds 0.5-2.9, 60 reps   data_summary_all.rds,         metadata_all.rds
#   fast_speeds      speeds 3.2-5.0, 20 reps   data_summary_fast_speeds.rds, metadata_fast_speeds.rds
# (load_data.R -d steady_tracking -t data,metadata --shift-speed 3.2,...,5.0 --tag fast_speeds)
# Sims that SLiM aborted (fitness underflowed to 0 once the population was ~1700 behind the
# optimum) are not in the data at all; extract_crashed_sims.sh lists them in crashed_sims.csv.
# In a combination where some reps aborted, the reps left are the ones that kept up best, so
# their mean is biased: read_summary() drops every combination with an aborted rep unless
# keep_partial = TRUE. All 20 reps aborted at speed 5.0 and in many runaway combinations at
# 3.2-4.7; plot_outcomes.R shows them.

suppressPackageStartupMessages(library(data.table))

st_data_dir <- "figure/r_data/steady_tracking"
st_files <- list(all = c("data_summary_all.rds", "metadata_all.rds"),
                 fast_speeds = c("data_summary_fast_speeds.rds", "metadata_fast_speeds.rds"))

# aborted sims (one row each) and their combinations
read_crashed <- function(dir = st_data_dir, speeds = plot_speeds) {
  f <- file.path(dir, "crashed_sims.csv")
  if (!file.exists(f)) return(data.table(shift_speed = numeric(), reg = character(), psi = numeric(),
                                         rep = integer(), ID = integer(), crash_tick = integer()))
  keep_speeds(fread(f), speeds)
}

# per-tick summary data of every run. cols: columns to keep (NULL = all; the parameter
# columns, rep, ID and tick are always kept). runs: names in st_files.
read_summary <- function(cols = NULL, runs = names(st_files), keep_partial = FALSE, dir = st_data_dir,
                         speeds = plot_speeds) {
  keys <- c("shift_speed", "reg", "n", "psi", "rep", "ID", "tick")
  out <- list()
  for (r in runs) {
    f <- file.path(dir, st_files[[r]][1])
    if (!file.exists(f)) { message("read_summary: no ", f, ", skipped"); next }
    x <- keep_speeds(as.data.table(readRDS(f)), speeds)
    if (!is.null(cols)) x <- x[, intersect(c(keys, cols), names(x)), with = FALSE]
    x[, run := r]
    out[[r]] <- x
    rm(x); gc(verbose = FALSE)
  }
  d <- rbindlist(out, use.names = TRUE, fill = TRUE)
  if (!keep_partial) d <- drop_partial(d)
  d
}

read_meta <- function(runs = names(st_files), keep_partial = FALSE, dir = st_data_dir, speeds = plot_speeds) {
  m <- rbindlist(lapply(runs, function(r) {
    f <- file.path(dir, st_files[[r]][2])
    if (file.exists(f)) keep_speeds(as.data.table(readRDS(f)), speeds)[, run := r]
  }), use.names = TRUE, fill = TRUE)
  if (!keep_partial) m <- drop_partial(m)
  m
}

# drop the combinations (speed x reg x psi) that have any aborted rep
drop_partial <- function(x, dir = st_data_dir) {
  cr <- unique(read_crashed(dir)[, .(shift_speed, reg, psi)])
  if (!nrow(cr) || !"shift_speed" %in% names(x)) return(x)
  key <- function(sp, reg, psi) paste(round(sp, 2), reg, round(psi, 2))
  x[!key(shift_speed, reg, psi) %in% key(cr$shift_speed, cr$reg, cr$psi)]
}

# page width for a plot with one facet column per speed
speed_width <- function(n_speeds, per = 1.6, extra = 3) max(12, extra + per * n_speeds)

# display names of the copying rules: use these in every figure (strips, legends, titles);
# the data keep ave / fit / ext
reg_names <- c(ave = "Copy-the-average", fit = "Copy-the-best", ext = "Copy-the-extreme")
reg_factor <- function(r) factor(reg_names[as.character(r)], levels = reg_names)

# Output root of the current analyses. Every plot script writes under plot_root; plots/ holds the
# earlier figures, which used the metadata's per-tick `failed` flag.
plot_root <- "analysis/steady_tracking/server/plots_fit500"

# Tracking outcome per sim (classify_tracking.R): failed = mean fitness over the last 500 ticks
# < 10% of the fitness at the optimum, or aborted by SLiM. Use this instead of metadata `failed`.
class_file <- file.path(plot_root, "tracking_class", "tracking_class_per_sim.csv")
read_class <- function(ds = "steady_tracking", speeds = plot_speeds) {
  if (!file.exists(class_file)) stop("no ", class_file, ": run classify_tracking.R first")
  keep_speeds(fread(class_file), speeds)[dataset == ds, .(ID, failed = outcome == "failed", class, fit_500, lag_500)]
}

# Speed range of every plot: the loaders above (read_summary, read_meta, read_crashed, read_class)
# keep only these speeds by default (speeds = NULL: all), and scripts that read their own files
# call keep_speeds() on them. classify_tracking.R classifies every speed.
plot_speeds <- c(1.1, 3.5)
keep_speeds <- function(x, speeds = plot_speeds) {
  if (is.null(speeds) || !"shift_speed" %in% names(x)) return(x)
  x[shift_speed >= speeds[1] - 1e-9 & shift_speed <= speeds[2] + 1e-9]
}

# rule display name -> data code (ave / fit / ext), e.g. for file names
reg_code <- function(x) names(reg_names)[match(as.character(x), reg_names)]

# ---- colours: the paper conventions (as in analysis/abrupt_shift) ----
# psi: MetBrewer Hokusai2 over psi 0-1 by 0.1; the stepped colour bar shows all 11 boxes and labels
# only the plotted psi. Map colour to numeric psi (psi_num() also takes a factor such as psi_f or
# "psi 0.2") and set group explicitly: a continuous colour does not split lines into groups.
# Group size n (reversed MetBrewer Greek) does not apply here: every steady_tracking sim has n = 8.
psi_grid <- seq(0, 1, by = 0.1)
psi_pal  <- as.character(MetBrewer::met.brewer("Hokusai2", length(psi_grid)))
psi_num  <- function(x) as.numeric(sub("^psi ", "", as.character(x)))
psi_scale <- function(show = c(0, 0.2, 0.5, 0.8), aesthetics = c("colour", "fill"), guide = TRUE) list(
  ggplot2::scale_colour_gradientn(
    expression(psi), colours = psi_pal, values = psi_grid, limits = range(psi_grid),
    breaks = show, labels = format(show, drop0trailing = TRUE), aesthetics = aesthetics,
    guide = if (guide) ggplot2::guide_colourbar(nbin = length(psi_grid), display = "rectangles",
                                                direction = "horizontal", title.position = "left",
                                                barwidth = grid::unit(4, "cm"),
                                                barheight = grid::unit(0.3, "cm"),
                                                position = "top") else "none"))

# rule colours, keyed by display name (rule = reg_factor(reg))
rule_cols <- setNames(c("#2a6fdb", "#e08a00", "#7a3fb0"), reg_names)
