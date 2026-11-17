# Shared settings for submit_slurm.sh and slurm_job.sh.
# Source it with the run name and the mode:   source pipeline_config.sh <run> <test|pilot|full>
# <run> names a file runs/<run>.sh next to this file that sets the grid, end time, reps and
# output folder of that run (see runs/README in that folder's header comments).
# Nothing here is read from the environment: every setting depends only on the run and the
# mode, so a variable left over in your shell (DIR_NAME, REPS, OUT, TIME_LIMIT, ...) cannot
# change a run.

RUN=${1:-}
MODE=${2:-}
CODE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
RUN_FILE=${CODE_DIR}/runs/${RUN}.sh
if [[ -z $RUN || ! -f $RUN_FILE ]]; then
    echo "pipeline_config.sh: unknown run '${RUN}'; runs: $(cd "${CODE_DIR}/runs" && ls *.sh | sed 's/\.sh$//' | tr '\n' ' ')" >&2
    return 1
fi
case $MODE in
    test|pilot|full) ;;
    *) echo "pipeline_config.sh: mode must be test, pilot or full (got '${MODE}')" >&2; return 1 ;;
esac
ME=$(id -un)

# this package: code/ (this folder) and pop/ (burnin pop files) next to it. A code snapshot
# made by submit_slurm.sh has the same layout, with pop/ a link to the real pop folder.
# -P: the real path, so POP_DIR stays on /home when entered through a symlink (e.g. /local/workdir)
PKG_DIR=$(cd "${CODE_DIR}/.." && pwd -P)

# ---- the run's own settings ----
# A run file sets:
#   RUN_DIR         output folder name on lab storage for full mode. Several runs may share one
#                   (e.g. fast_speeds adds speeds to steady_tracking) as long as their
#                   combinations differ; test and pilot always go to test_<run> / pilot_<run>
#   END_TIME        last tick (the optimum starts moving at tick 10000)
#   ID_OFFSET       added to every sim ID, so IDs never repeat across runs (0 = the original run)
#   N_REPS_FULL     replicates per combination in full mode (at most 60: there are 60 burnins)
#   REPS_PER_TASK   replicates one array task runs; a combination's reps are split over
#                   ceil(N_REPS_FULL / REPS_PER_TASK) tasks
#   TIME_LIMIT_H_TEST / _PILOT / _FULL   Slurm time limit (hours) per mode
#   and either the grid  SPEEDS NS PSIS REGS  (every combination)
#   or a list            COMBOS=("<speed> <n> <psi> <reg>" ...)
# Write speeds and psi exactly as in the pop / folder names: 1.1 not 1, 0.0 not 0.
COMBOS=()
source "$RUN_FILE"

case $MODE in
    test)  DIR_NAME=test_${RUN} ;;
    pilot) DIR_NAME=pilot_${RUN} ;;
    full)  DIR_NAME=${RUN_DIR} ;;
esac

# all outputs, including every log, go to permanent Messer lab storage, not /home (200 GB quota).
# It is an NFS mount that every job must create itself (BioHPC docs): call mount_lab_storage
LAB_NODE=cbsubscb10
OUT=/fs/${LAB_NODE}/storage/${ME}/${DIR_NAME}
mount_lab_storage() { /programs/bin/labutils/mount_server ${LAB_NODE} /storage; }
# burnin pop files (input): cons_n_<n>_psi_<psi>_reg_<reg>_<rep>.pop, one per n/psi/reg/rep;
# every shift speed starts from the same burnin. The package carries the 720 this grid needs
# (n = 8, psi 0.0/0.2/0.5/0.8, 3 regs, reps 1-60; also in /fs/cbsubscb10/storage/<user>/cons_env_pop)
POP_DIR=${PKG_DIR}/pop

# SLURM opens a job's --output file before the job can mount lab storage, so it has to be on
# /home. The job sends everything to its own log instead, copies that to $OUT/logs and then
# deletes the /home file; one is only left here if the job died before it could do that.
LOG_DIR=/home/${ME}/slurm_logs
# frozen copy of the package made at each submit (<CODE_ROOT>/<DIR_NAME>_<time>/{code,pop});
# the jobs run that copy, so editing or re-syncing code/ while tasks are queued cannot
# change their code (also kept at $OUT/code)
CODE_ROOT=/home/${ME}/slurm_code
JOB_NAME=${RUN}_${MODE}

# A folder can hold several runs. Each run's frozen code goes to its own folder: code/ for the
# run the folder is named after (and in test_ / pilot_ folders, which belong to one run), and
# code_<run>/ for a run added to another run's folder, so the record of the code that produced
# the folder's original data is never replaced.
if [[ $MODE != full || $RUN == "$RUN_DIR" ]]; then
    CODE_OUT=${OUT}/code
    FOLDER_OWNER=$RUN
else
    CODE_OUT=${OUT}/code_${RUN}
    FOLDER_OWNER=$RUN_DIR
fi
# each task's summary/<combo>/_RUN names the run that wrote it; a summary without one (written
# before runs existed) belongs to the folder's own run. A run may only (re)write its own tasks.
#   combo_owner <combo>  prints the owning run, or nothing if the task is not in $OUT yet
combo_owner() {
    local d=${OUT}/summary/$1
    if [[ -f $d/_RUN ]]; then cat "$d/_RUN"
    elif [[ -d $d || -f ${OUT}/raw/$1.tar.zst ]]; then echo "$FOLDER_OWNER"
    fi
}

# ---- combinations ----
# COMBO_LIST holds "<speed> <n> <psi> <reg>" per combination. From a grid it is built in the
# order speed > n > psi > reg (reg varies fastest), the order of the original run, so that run's
# sim IDs are reproduced exactly.
COMBO_LIST=()
if (( ${#COMBOS[@]} )); then
    COMBO_LIST=("${COMBOS[@]}")
else
    for _sp in "${SPEEDS[@]}"; do for _n in "${NS[@]}"; do for _psi in "${PSIS[@]}"; do for _reg in "${REGS[@]}"; do
        COMBO_LIST+=("$_sp $_n $_psi $_reg")
    done; done; done; done
fi
N_COMBOS=${#COMBO_LIST[@]}

# ---- replicates and tasks per mode ----
#   test:  the first 3 combinations (put one of each reg first in a COMBOS list) x 2 reps
#   pilot: every combination x 3 reps, one task each
#   full:  every combination x N_REPS_FULL reps, split into tasks of REPS_PER_TASK reps
case $MODE in
    test)  REPS=2;            RPT=2;              TIME_LIMIT_H=$TIME_LIMIT_H_TEST ;;
    pilot) REPS=3;            RPT=3;              TIME_LIMIT_H=$TIME_LIMIT_H_PILOT ;;
    full)  REPS=$N_REPS_FULL; RPT=$REPS_PER_TASK; TIME_LIMIT_H=$TIME_LIMIT_H_FULL ;;
esac
(( RPT > REPS )) && RPT=$REPS
N_CHUNKS=$(( (REPS + RPT - 1) / RPT ))
if [[ $MODE == test ]]; then
    N_TASKS=$(( N_COMBOS < 3 ? N_COMBOS : 3 ))
else
    N_TASKS=$(( N_COMBOS * N_CHUNKS ))
fi
TIME_LIMIT=${TIME_LIMIT_H}:00:00
# largest array this cluster is known to accept (fluctuate_env_constant_psi)
MAX_ARRAY=4950

# task index -> parameters, the same in every script. Tasks cycle through the combinations
# first, then the rep chunks (task = chunk * N_COMBOS + combo index), so a test's first tasks
# are different combinations.
#   task_params <task>  sets combo_i sp n psi reg, rep_first rep_last, combo (the name used for
#                       summary/, logs/ and raw/), part_prefix (parquet part file prefix)
task_params() {
    local t=$1
    combo_i=$(( t % N_COMBOS ))
    local chunk=$(( t / N_COMBOS ))
    read -r sp n psi reg <<< "${COMBO_LIST[$combo_i]}"
    rep_first=$(( chunk * RPT + 1 ))
    rep_last=$(( rep_first + RPT - 1 ))
    (( rep_last > REPS )) && rep_last=$REPS
    combo="sp_${sp}_reg_${reg}_n_${n}_psi_${psi}"
    part_prefix="psi_${psi}-"
    if (( N_CHUNKS > 1 )); then
        combo+="_reps_${rep_first}-${rep_last}"
        part_prefix="psi_${psi}-reps_${rep_first}-${rep_last}-"
    fi
}
# unique sim ID of combination combo_i, replicate rep: rep > combination order, plus 2 (as in
# the original run), plus the run's ID_OFFSET
sim_id() { echo $(( ID_OFFSET + (($2 - 1) * N_COMBOS + $1) + 2 )); }

# nodes to avoid. In the fluc_env_constant_psi run (Oct 2026) no task finished on
# cbsubscb02-12: polars crashes there with "Illegal instruction" (the CPUs lack AVX2 etc.) and
# 02-07 were also too slow. Sim times in the steady_tracking run (2000 ticks after the shift):
# mean 640 s, max 1740 s per sim, about the same at every speed; nodes differ by > 2x.
EXCLUDE_NODES='cbsubscb[02-04,06-12]'

# tsv output folders converted to parquet; log/ holds .txt and is only archived
TSV_TYPES=(phenotype genotype data allele group_tag metadata)
ARCHIVE_TYPES=("${TSV_TYPES[@]}" log)

SLIM=/home/${ME}/bin/slim5.0
PY=/home/${ME}/miniconda3/envs/simdata/bin/python
ZSTD=/home/${ME}/miniconda3/envs/simdata/bin/zstd

# the files that decide what a run computes. Pilot and full submits are refused unless the
# test run of the same run used the same ones. Comment lines of slurm_job.sh (incl. its
# #SBATCH defaults), of this file and of the run file, and the time-limit and EXCLUDE_NODES
# lines, may differ.
code_fingerprint() {
    { md5sum base_script.slim make_model.py convert_parquet.py "$SLIM"
      grep -v '^[[:space:]]*#' slurm_job.sh | md5sum
      grep -vE '^[[:space:]]*#|^EXCLUDE_NODES=|TIME_LIMIT' pipeline_config.sh | md5sum
      grep -vE '^[[:space:]]*#|TIME_LIMIT' "runs/${RUN}.sh" | md5sum
    } | md5sum | cut -d' ' -f1
}
