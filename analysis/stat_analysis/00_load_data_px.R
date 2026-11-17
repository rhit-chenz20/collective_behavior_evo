# Load the data for section 10 (psi = x vs psi = 0). Sourced by 10_psi_vs_psi0.R and
# 10_psi_vs_psi0_fluc.R after 00_setup.R. Each file is read only if its object is not
# already loaded. These are separate from cons_data_df (sections 1-9), which is a
# different file. rm(px_data_df) / rm(fluc_data_df) to reload after changing a file.
#   px_data_df   constant env, all ticks: fitness_mean, z_mean, g_mean   (px_cons_file)
#   fluc_data_df fluctuating env, every tick: fitness_mean               (px_fluc_file)
# Set px_load_fluc <- FALSE before sourcing to skip the fluctuating data.

if (!exists("px_load_fluc")) px_load_fluc <- TRUE

if (!exists("px_data_df")) px_data_df <- read_rds_data(px_cons_file)
if (px_load_fluc && !exists("fluc_data_df")) fluc_data_df <- read_rds_data(px_fluc_file)
