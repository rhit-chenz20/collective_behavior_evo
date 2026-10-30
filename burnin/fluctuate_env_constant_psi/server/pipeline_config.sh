# Shared settings for submit_slurm.sh and slurm_job.sh.
# Values already set in the environment (e.g. DIR_NAME=test_run via --export) take precedence.

PROJ=${PROJ:-/home/$USER/collective_behavior_evo}   # code
DIR_NAME=${DIR_NAME:-fluc_env_adjusted_shift_size}

# all outputs go to permanent Messer lab storage, not /home (200 GB quota).
# It is an NFS mount that every job must create itself (BioHPC docs): call mount_lab_storage
LAB_NODE=cbsubscb10
OUT=${OUT:-/fs/${LAB_NODE}/storage/$USER/${DIR_NAME}}
mount_lab_storage() { /programs/bin/labutils/mount_server ${LAB_NODE} /storage; }


# parameter grid. Task index order: gap > sz > n > psi > reg (reg varies fastest), then the
# chunk of replicates. Sims start from scratch (no burnin file): ~31 min each in the test run.
GAPS=(20 50 80 110)   # fluc_interval: ticks between optimum flips
SZS=(9.4 28.2 47.0 65.8 84.6 103.4)     # 1x 3x 5x 7x 9x 11x the base shift size 9.4
NS=(2 5 8 10 20 50)
PSIS=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
REGS=(ave fit ext)
N_REPS_FULL=20               # fixes the sim IDs; do not change
REPS=${REPS:-$N_REPS_FULL}   # replicates per combination (fewer only for tests)
REPS_PER_TASK=${REPS_PER_TASK:-10}   # 10 sims took 4.5-10.4 h on average per node (Oct 2026), under TIME_LIMIT
N_CHUNKS=$(( (REPS + REPS_PER_TASK - 1) / REPS_PER_TASK ))
N_COMBOS=$(( ${#GAPS[@]} * ${#SZS[@]} * ${#NS[@]} * ${#PSIS[@]} * ${#REGS[@]} ))   # 4752
N_TASKS=$(( N_COMBOS * N_CHUNKS ))                                               # 9504
# tasks per sbatch call; 4950 is the largest array this cluster is known to accept
MAX_ARRAY=${MAX_ARRAY:-4950}

# nodes and time limit for every sbatch call (submit_slurm.sh, rerun_failed.sh). In the
# fluc_env_constant_psi run (Oct 2026) no task finished on cbsubscb02-12: polars crashes there
# with "Illegal instruction" (the CPUs lack AVX2 etc.) and 02-07 also hit the time limit.
# On cbsubscb13-20 a task (10 sims) took 4.5-10.4 h on average per node and up to ~12 h,
# so 12 h was too short; 24 h leaves room. The test run of this grid (2026-10-26, gap 20,
# sz 9.4, n 2) took ~32 min/sim on cbsubscb16, in line with the old grid there (~37 min).
# Memory (--mem=4G in slurm_job.sh): max RSS was 334 MB over all 11,882 old tasks, but SLURM
# samples every 30 s and can miss the short parquet conversion peak, so 4G is kept.
EXCLUDE_NODES=${EXCLUDE_NODES:-cbsubscb[02-04,06-12]}
TIME_LIMIT=${TIME_LIMIT:-24:00:00}

# tsv output folders converted to parquet (the tsv themselves are not kept); log/ holds
# SLiM's .txt output and is only archived.
# vars, copied_val and fitness are not recorded: they can be recomputed from
# genotype + group_tag + phenotype (+ the optimum, from data: z_opt + z_mean)
TSV_TYPES=(phenotype genotype data allele group_tag metadata)
ARCHIVE_TYPES=("${TSV_TYPES[@]}" log)

PY=/home/$USER/miniconda3/envs/simdata/bin/python
ZSTD=/home/$USER/miniconda3/envs/simdata/bin/zstd
