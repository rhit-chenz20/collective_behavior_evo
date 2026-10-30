# Shared settings for submit_slurm.sh and the slurm_*.sh job scripts.
# Values already set in the environment (e.g. DIR_NAME=test_run via --export) take precedence.

PROJ=${PROJ:-/home/$USER/collective_behavior_evo}   # code
DIR_NAME=${DIR_NAME:-cons_env_moving_opt}

# all outputs go to permanent Messer lab storage, not /home (200 GB quota).
# It is an NFS mount that every job must create itself (BioHPC docs): call mount_lab_storage
LAB_NODE=cbsubscb10
OUT=${OUT:-/fs/${LAB_NODE}/storage/$USER/${DIR_NAME}}
mount_lab_storage() { /programs/bin/labutils/mount_server ${LAB_NODE} /storage; }
# burnin pop files (input), also on lab storage: cons_n_<n>_psi_<psi>_reg_<reg>_<rep>.pop
# (one per n/psi/reg/rep, >= 60 reps each; every shift speed starts from the same burnin)
POP_DIR=${POP_DIR:-/fs/${LAB_NODE}/storage/$USER/cons_env_pop}

# parameter grid: one array task per (shift speed, n, psi, reg) combination, running all
# replicates in a row. Task index order: speed > n > psi > reg (reg varies fastest).
# The optimum moves by shift_speed every tick from tick 10000 (stops at 11000 or when
# mean fitness drops below half the fitness at the optimum, see base_script.slim).
SPEEDS=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)
NS=(2 4 5 8 10 20 25 40 50 100)
PSIS=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
REGS=(ave ext fit)
N_REPS_FULL=60          # fixes the sim IDs; do not change
REPS=${REPS:-$N_REPS_FULL}   # replicates per task (fewer only for tests)
N_COMBOS=$(( ${#SPEEDS[@]} * ${#NS[@]} * ${#PSIS[@]} * ${#REGS[@]} ))   # 3300

# nodes and time limit for every sbatch call. In the fluc_env_constant_psi run (Oct 2026) no
# task finished on cbsubscb02-12: polars crashes there with "Illegal instruction" (the CPUs
# lack AVX2 etc.) and 02-07 were also too slow. Sims took 380-600 s each in the
# cons_env_data_with_adjusted_shift_size run (60 reps ~6.5-10 h); sims that stop early are shorter.
EXCLUDE_NODES=${EXCLUDE_NODES:-cbsubscb[02-04,06-12]}
TIME_LIMIT=${TIME_LIMIT:-24:00:00}

# tsv output folders converted to parquet; log/ holds .txt and is only archived
TSV_TYPES=(phenotype genotype data vars fitness allele copied_val group_tag metadata)
ARCHIVE_TYPES=("${TSV_TYPES[@]}" log)

PY=/home/$USER/miniconda3/envs/simdata/bin/python
ZSTD=/home/$USER/miniconda3/envs/simdata/bin/zstd

# only used by slurm_convert.sh / slurm_archive.sh, the old flow for loose tsv folders
# already on disk; 1 = delete them once the parquet and the archive are verified
DELETE_RAW=${DELETE_RAW:-1}
