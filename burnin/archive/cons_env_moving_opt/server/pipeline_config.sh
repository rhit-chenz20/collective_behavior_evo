# Shared settings for submit_slurm.sh and the slurm_*.sh job scripts.
# Source it with the run mode:   source pipeline_config.sh <test|full>
# Nothing here is read from the environment: every setting depends only on the mode, so a
# variable left over in your shell (DIR_NAME, REPS, OUT, ...) cannot change a run.

if [[ ${1:-} != test && ${1:-} != full ]]; then
    echo "pipeline_config.sh: mode must be 'test' or 'full' (got '${1:-}')" >&2
    return 1
fi
MODE=$1
ME=$(id -un)

PROJ=/home/${ME}/collective_behavior_evo   # code
DIR_NAME=cons_env_moving_opt
[[ $MODE == test ]] && DIR_NAME=${DIR_NAME}_test

# all outputs, including every log, go to permanent Messer lab storage, not /home (200 GB quota).
# It is an NFS mount that every job must create itself (BioHPC docs): call mount_lab_storage
LAB_NODE=cbsubscb10
OUT=/fs/${LAB_NODE}/storage/${ME}/${DIR_NAME}
mount_lab_storage() { /programs/bin/labutils/mount_server ${LAB_NODE} /storage; }
# burnin pop files (input), also on lab storage: cons_n_<n>_psi_<psi>_reg_<reg>_<rep>.pop
# (one per n/psi/reg/rep, >= 60 reps each; every shift speed starts from the same burnin)
POP_DIR=/fs/${LAB_NODE}/storage/${ME}/cons_env_pop

# SLURM opens a job's --output file before the job can mount lab storage, so it has to be on
# /home. The job sends everything to its own log instead, copies that to $OUT/logs and then
# deletes the /home file; one is only left here if the job died before it could do that.
LOG_DIR=/home/${ME}/slurm_logs
# frozen copy of this folder made at each submit; the jobs run that copy, so editing this
# folder while tasks are queued cannot change their code (also kept at $OUT/code)
CODE_ROOT=/home/${ME}/slurm_code

# parameter grid: one array task per (shift speed, n, psi, reg) combination, running all
# replicates in a row. Task index order: speed > n > psi > reg (reg varies fastest).
# The optimum moves by shift_speed every tick from tick 10000 (stops at 11000 or when
# mean fitness drops below half the fitness at the optimum, see base_script.slim).
SPEEDS=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)
NS=(2 4 5 8 10 20 25 40 50 100)
PSIS=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
REGS=(ave ext fit)
N_REPS_FULL=60          # fixes the sim IDs; do not change
N_COMBOS=$(( ${#SPEEDS[@]} * ${#NS[@]} * ${#PSIS[@]} * ${#REGS[@]} ))   # 3300
if [[ $MODE == full ]]; then
    REPS=$N_REPS_FULL
    N_TASKS=$N_COMBOS
else
    # tasks 0-2: speed 0.01 (the optimum moves slowest, so sims track it longest: the slowest
    # case), n 2, psi 0.0, one task per reg, so all three SLiM scripts run
    REPS=2
    N_TASKS=3
fi
# tasks per sbatch call; 4950 is the largest array this cluster is known to accept
# (fluctuate_env_constant_psi); larger runs are split into arrays of at most this many
MAX_ARRAY=4950

# task index -> parameters, the same in every script:
#   task_params <task>  sets sp_i n_i psi_i reg_i, sp n psi reg, combo, and n_sp n_n n_psi n_reg
task_params() {
    local r=$1
    n_sp=${#SPEEDS[@]}; n_n=${#NS[@]}; n_psi=${#PSIS[@]}; n_reg=${#REGS[@]}
    reg_i=$(( r % n_reg ));  r=$(( r / n_reg ))
    psi_i=$(( r % n_psi ));  r=$(( r / n_psi ))
    n_i=$(( r % n_n ));      r=$(( r / n_n ))
    sp_i=$r
    reg=${REGS[$reg_i]}; psi=${PSIS[$psi_i]}; n=${NS[$n_i]}; sp=${SPEEDS[$sp_i]}
    combo="sp_${sp}_reg_${reg}_n_${n}_psi_${psi}"
}

# nodes and time limit for every sbatch call. In the fluc_env_constant_psi run (Oct 2026) no
# task finished on cbsubscb02-12: polars crashes there with "Illegal instruction" (the CPUs
# lack AVX2 etc.) and 02-07 were also too slow.
# Time limit: sims to tick 11000 (1000 ticks after the shift) took 380-600 s each in the
# cons_env_data_with_adjusted_shift_size run. Sims here may run 2000 ticks after the shift, so
# assume up to ~2x that (~1200 s): 60 reps ~20 h at worst, hence 24 h for the full run.
# The test runs only 2 reps per task (~40 min at worst), so it gets 2 h: a short limit
# schedules faster and shows a hung sim sooner. The test records each sim's time
# (summary/<combo>/sim_seconds.tsv) and `submit_slurm.sh full` checks the 60-rep projection
# against the full limit before submitting.
EXCLUDE_NODES='cbsubscb[02-04,06-12]'
if [[ $MODE == full ]]; then TIME_LIMIT_H=24; else TIME_LIMIT_H=2; fi
TIME_LIMIT=${TIME_LIMIT_H}:00:00

# tsv output folders converted to parquet; log/ holds .txt and is only archived
TSV_TYPES=(phenotype genotype data vars fitness allele copied_val group_tag metadata)
ARCHIVE_TYPES=("${TSV_TYPES[@]}" log)

PY=/home/${ME}/miniconda3/envs/simdata/bin/python
ZSTD=/home/${ME}/miniconda3/envs/simdata/bin/zstd

# the files that decide what a run computes. `submit_slurm.sh full` is refused unless the
# test run used the same ones (the TIME_LIMIT/EXCLUDE_NODES lines of this file may differ)
code_fingerprint() {
    { md5sum base_script.slim make_model.py slurm_job.sh convert_parquet.py bin/slim5.0
      grep -vE '^(TIME_LIMIT|EXCLUDE_NODES)=|TIME_LIMIT_H=[0-9]' pipeline_config.sh | md5sum
    } | md5sum | cut -d' ' -f1
}

# only used by slurm_convert.sh / slurm_archive.sh, the old flow for loose tsv folders
# already on disk; 1 = delete them once the parquet and the archive are verified
DELETE_RAW=1
