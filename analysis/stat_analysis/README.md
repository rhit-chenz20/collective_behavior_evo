# stat_analysis

Statistical comparison of the social rules (`ave`, `fit`, `ext`) across psi and group size n.
Split from `analysis/stat_analysis.Rmd` (sections 1–10). Run every script from the project
root, e.g. `source("analysis/stat_analysis/06_genetic_variance.R")`.

Outputs go to `analysis/figure/<date>/` (`date` is set in `00_setup.R`). Sections 1–5 prefix
their files with the trait (`g_` or `z_`); sections 6–10 write to their own subfolder.

## Setup and data loading

| File | What it does |
|---|---|
| `00_setup.R` | Parameters (`date`, `trait`, `n_boot`, `models`, ...), packages, data file paths, and helpers shared by several sections (`boot_diff_ratio`, `classify_equiv`, `ratio_heatmap`, `boot_mean_ci`, `fmt`, ...). Every script sources it. |
| `00_load_data.R` | Sections 1–9: reads `cons_data_df` (only if not already loaded) and builds `burn_df`, the burn-in trajectories (ticks < 10000), one shift size per run. |
| `00_load_data_px.R` | Section 10: reads `px_data_df` (constant env) and, unless `px_load_fluc <- FALSE`, `fluc_data_df` (fluctuating env). |

Data files are set in `00_setup.R`. They are made by
`analysis/load_data_from_storage/load_data.R`; see `load_examples.sh` there for commands.

| Object | File (under `analysis/rds/`) | Needs |
|---|---|---|
| `cons_data_df` | `cons_env_data_with_adjusted_shift_size/data_summary_all_nosfs` | type `data`, all sims; columns `tick, psi, n, reg, rep, shift_size, z_mean, z_sd, g_mean, g_sd, fitness_mean` (`-t data --tag all_nosfs`, ~19 GB) |
| `px_data_df` | `cons_env_data_with_adjusted_shift_size/data_summary_fit_zg_all` | type `data`, all sims, `fitness_mean, z_mean, g_mean` |
| `fluc_data_df` | `fluc_env_constant_psi/data_summary_fitness_all` | type `data`, all sims, `fitness_mean` at **every tick** (no `--gen-values` thinning) |

## Sections

Sections 1–9 use the constant-environment burn-in (ticks 9950–9999, stabilizing selection).
A run is one psi × n × reg × rep and is the unit of analysis.
`trait` in `00_setup.R` picks phenotype (`z`) or genotype (`g`) for sections 1–5.

| File | Analysis and plots | Needs (besides setup) | Main outputs |
|---|---|---|---|
| `01_burnin_diagnostics.R` | Trajectory of the trait mean and SD over the burn-in, one line per run, by psi × n, coloured by model. | `cons_data_df` | `<trait>_burnin_<sd/mean col>.pdf` |
| `02_per_run_summary.R` | Collapses each run to single values: SD and mean averaged over its time points (`mean_sd_all`, `mean_mean_all`). Checks every run has 50 time points. | `cons_data_df` | `run_df`, `<trait>_per_run_summary.csv` |
| `03_model_comparison.R` | Model comparison within each psi × n: group stats, ANOVA, pairwise Welch t-tests (Holm-adjusted), bootstrap difference and SD ratio with 90%/95% CIs. | `run_df` (runs 02 if missing) | `<trait>_model_group_stats.csv`, `_model_anova.csv`, `_model_pairwise.csv` |
| `04_labeled_summary.R` | Readable text summary of section 3. | section 3 results (runs 03 if missing) | `<trait>_model_comparison_summary.txt` |
| `05_effect_size_report.R` | SD ratios ave/fit and ext/fit as psi × n heatmaps, classified Equivalent / Different / inconclusive against a ±5% margin (90% CI = TOST), plus a text report. | section 3 results (runs 03 if missing) | `<trait>_sd_ratio_equivalence.csv`, `_sd_ratio_heatmap_*.pdf`, `_sd_ratio_report.txt` |
| `06_genetic_variance.R` | Genetic SD (`g_sd`, whatever `trait` is): vs psi per n and model; selected-n versions with several legend styles (6.0.1); relative to the neutral ceiling sqrt(2 N_e V_m) (6.0.2); ratio to psi = 0 (6.1); trend tests in psi (Spearman, Jonckheere–Terpstra, slope) (6.2); relaxed-selection mechanism check, log V_g vs −2 log(1 − psi) (6.3); neutral ceiling at psi = 1 (6.4); group-size effect (6.5); text report (6.6). | `cons_data_df` | `genetic_variance/`: `gv_per_run.csv`, `g_sd_vs_psi*.pdf`, `g_sd_ratio_vs_psi0.csv`, `g_sd_trend_tests.csv`, `mechanism_*.{pdf,csv}`, `neutral_ceiling_psi1.csv`, `group_size_tests.csv`, `genetic_variance_report.txt` |
| `07_phenotypic_variance.R` | Phenotypic SD (`z_sd`) and its link to `g_sd`: theoretical variance factor Var(P)/V_g for each rule (7.0); ratio to psi = 0 (7.1); z_sd and z_sd/g_sd vs n with theory (7.2); z_sd/g_sd decomposition (7.3); predicted vs observed phenotypic variance (7.4); text report (7.5). | `cons_data_df` | `phenotypic_variance/`: `theory_variance_factor.csv`, `z_sd_*`, `zg_ratio_*`, `pred_vs_obs_z_var.pdf`, `prediction_*.csv`, `phenotypic_variance_report.txt` |
| `08_run_heterogeneity.R` | Run-to-run CV of per-run z_sd and g_sd per psi × n; flags cells with CV > 0.2 and plots their run trajectories. | `pv_run`, `pv_df` (runs 07 if missing) | `heterogeneity/`: `cv_of_per_run_sd.csv`, `cv_heatmap_*.pdf`, `flagged_trajectories_*.pdf`, `heterogeneity_report.txt` |
| `09_fitness_load.R` | Mean relative fitness and genetic load (1 − fitness) per psi × n, load relative to psi = 0, and observed vs Gaussian-predicted load from z_mean and z_sd. | `cons_data_df` (needs `fitness_mean`) | `fitness_load/`: `fitness_per_run.csv`, `fitness_load_by_cell.csv`, `fitness_heatmap_*.pdf`, `load_*.{pdf,csv}`, `fitness_load_report.txt` |

### Section 10: psi = x vs psi = 0

For each model × n and each metric, compares the mean over runs at every psi = x > 0 with
psi = 0. A cell's colour is the fraction of x values better than psi = 0 (blue = all
better, red = all worse); the label gives how many. All psi, including 1, are compared (`px_psi_exclude`). Runs that do not reach a "Gens to" threshold by generation 11000 count as slowest (Inf), so a psi = x with any such run is worse than psi = 0; such psi are listed in `psi_unreached` in the CSV.

Outputs go to `analysis/figure/<date>/psi_vs_psi0/`, in subfolders (set by `px_out()` in `00_setup.R`): `constant_only/` (constant-only tables), `results/` (the full comparison CSV that every figure is drawn from, plus its settings and the peak / worst phase CSV), `overview/` (combined table and per-metric fluctuating plots), `paper/` (paper figures).

| File | Analysis and plots | Needs | Main outputs (`psi_vs_psi0/<subfolder>/`) |
|---|---|---|---|
| `10_psi_vs_psi0.R` | Constant environment only. Metrics: mean fitness over the burn-in (stabilizing); fitness at generation 10000, generations until z_mean / g_mean reach 50% and 90% of the shift, and "Phenotype lag behind optimum" over the first 25 and 100 generations (`px_lag_T`): the sum over the first T generations after the shift of (new optimum − z_mean) / shift_size, in generations, lower = faster, not biased by generation-to-generation noise like the crossing time (directional, per shift size). The median and fluctuating versions below include the same lag columns. | `px_data_df` | `constant_only/psi_vs_psi0_table.csv`, `.pdf` |
| `10_psi_vs_psi0_median.R` | As `10_psi_vs_psi0.R`, but the "generations to reach 50% / 90% of the shift" metrics compare the **median** over runs instead of the mean (fitness metrics still use the mean), because the time to reach varies a lot between runs (especially `ext`). With `px_unreached_as_inf = TRUE` (default), runs that do not reach the threshold by generation 11000 count as slowest, so a psi is worse than psi = 0 if its median run did not reach it; `FALSE` leaves out a psi with any unreached run, as in the mean version. Median-based metrics are marked `[median]` in the plot. | `px_data_df` | `constant_only/psi_vs_psi0_median_table.csv` (`statistic`, `value_psi_*` = compared value, `n_runs_unreached`), `.pdf` |
| `10_psi_vs_psi0_fluc.R` | The same, plus fluctuating-environment panels, one per (fluc_opt, fluc_interval). Fitness is summarised per half cycle (flip to flip; flips assumed at 50 + k × interval), using only whole half cycles and an even number of them: cycle average; peak = half-cycle max and mean of its last `px_edge_k` (5) ticks; worst = half-cycle min and mean of its first 5 ticks. Each is averaged over half cycles per run. Writes the CSVs, then runs the plot script. | `px_data_df`, `fluc_data_df` | `results/psi_vs_psi0_table_with_fluc.csv` (every setting, metric and n), `_settings.csv`, `psi_vs_psi0_fluc_peak_worst_phase.csv` (where max/min fall in the half cycle) |
| `10_psi_vs_psi0_fluc_plot.R` | Plots from those CSVs only (no data, nothing recomputed): combined table with the fluctuating settings chosen by `px_fluc_opt_main` / `px_fluc_interval_main`, and one overview per fluctuating metric over all settings. | the two CSVs above | `overview/psi_vs_psi0_table_with_fluc.pdf`, `overview/psi_vs_psi0_fluc_<metric>.pdf` |
| `10_psi_vs_psi0_figure.R` | Paper figures from the same CSVs (nothing recomputed). The directional shift size is matched to the fluctuating `fluc_opt` (nearest `fluc_opt` with a warning until matched runs exist), and only n values present in both environments are shown (`px_common_n`). **Main:** Stabilizing, Directional at `px_main_shift` (18.8), one Fluctuating panel per interval in `px_main_interval` (40, 80) at `fluc_opt` = that shift (or `px_main_fluc_opt`); metrics chosen with `px_main_metrics` (NULL = all); the 5-tick peak / worst versions are dropped (`px_drop_metrics`); text size `px_text_scale` (3), tile size `px_col_width` / `px_row_height`; only the `px_reach_keep` (50%) "Gens to X% of shift" threshold is shown in both figures. **Extra main-style figures** (`px_pair_figs`, file name -> fluctuating (fluc_opt, interval) pairs): e.g. `psi_vs_psi0_main_opt_interval_pairs.pdf` with (20, 20) and (80, 80). **Supplement:** one file per shift size: Stabilizing, Directional at that shift, its matched fluctuation at every interval, all metrics. | the two CSVs above | `paper/psi_vs_psi0_main.pdf`, `paper/psi_vs_psi0_main_opt_interval_pairs.pdf`, `paper/psi_vs_psi0_supp_shift_<S>.pdf` |
| `10_psi_vs_psi0_figure_n8.R` | Compact paper figure for one group size (`px_n`, 8): models (ave, fit, ext) as columns, one row per selection × metric (Stabilizing, Directional at `px_shift` (18.8), then each fluctuating setting). One figure per entry of `px_n8_figs` (file name -> fluctuating (fluc_opt, interval) rows): opt 40 at intervals 40 / 80, and (20, 20) / (80, 80). Same threshold, dropped metrics and text-size settings as `10_psi_vs_psi0_figure.R`. | the two CSVs above | `paper/psi_vs_psi0_n8_main.pdf`, `paper/psi_vs_psi0_n8_opt_interval_pairs.pdf` |

## Notes

- Results from one section kept in the R session (`run_df`, `pairwise_df`, `pv_run`) are
  reused by later ones. After changing `trait`, `rm(run_df, pairwise_df)` or start a new session.
- Data objects are read only once per session; `rm()` them to re-read after replacing a file.
