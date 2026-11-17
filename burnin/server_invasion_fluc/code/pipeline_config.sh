# Shared settings for submit_slurm.sh and slurm_job.sh (invasion under a fluctuating optimum).
# Two kinds of run (RUN_KIND in the run file):
#   burnin    psi = 0 populations from scratch to tick 9949 under each fluctuation regime, one
#             per (fluc_opt, fluc_interval, rep): <RUN_DIR>/pop/fluc_sz_<sz>_gap_<gap>_<rep>.pop
#             (+ .pop.L). No parquet; only the pop files, SLiM logs and summaries.
#   invasion  read those, put INVASION_IND social invaders in at the first flip at or after
#             10000, run to END_TIME or until they are lost or fixed; parquet as usual.
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

# this package: code/ (this folder). A code snapshot made by submit_slurm.sh has the same layout.
# -P: the real path, so it stays on /home when entered through a symlink (e.g. /local/workdir)
PKG_DIR=$(cd "${CODE_DIR}/.." && pwd -P)

# ---- the run's own settings ----
# A run file sets:
#   RUN_KIND        burnin or invasion
#   RUN_DIR         output folder name on lab storage for full mode. Several runs may share one
#                   as long as their combinations differ (inv10 / inv500 share invasion_fluc);
#                   test and pilot always go to test_<run> / pilot_<run>
#   END_TIME        last tick (invasion: a sim stops earlier once the invaders are lost or fixed;
#                   burnin: ignored, the burnin stops at 9949)
#   ID_OFFSET       added to every sim ID, so IDs never repeat across runs
#   N_REPS_FULL     replicates per combination in full mode
#   REPS_PER_TASK   replicates one array task runs; a combination's reps are split over
#                   ceil(N_REPS_FULL / REPS_PER_TASK) tasks
#   TIME_LIMIT_H_TEST / _PILOT / _FULL   Slurm time limit (hours) per mode
#   the grid  SZS GAPS NS PSIS REGS  (every combination; burnin: NS=(8) PSIS=(0.0) REGS=(ave),
#             as the residents do not copy and the rule does not matter)
#   or a list COMBOS=("<sz> <gap> <n> <psi> <reg>" ...)
#   invasion only:
#     INVASION_IND  number of social invaders (one per run)
#     BURNIN_RUN    the burnin run whose pop files the residents come from (its full-mode
#                   folder, also in test and pilot), and N_POP_REPS its N_REPS_FULL: rep r
#                   starts from burnin ((r - 1) % N_POP_REPS) + 1, with its own random seed
# Write sz and psi exactly as in the folder names: 47.0 not 47, 0.0 not 0.
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
# resident burnin pop files: a burnin run writes them to $OUT/pop; an invasion run reads them
# from the full-mode folder of its BURNIN_RUN (on lab storage, mounted by every job)
case ${RUN_KIND:-} in
    burnin)   POP_DIR=${OUT}/pop ;;
    invasion) POP_DIR=$(source "${CODE_DIR}/runs/${BURNIN_RUN}.sh" && echo "/fs/${LAB_NODE}/storage/${ME}/${RUN_DIR}/pop") ;;
    *) echo "pipeline_config.sh: RUN_KIND must be burnin or invasion in ${RUN_FILE}" >&2; return 1 ;;
esac
pop_name() { echo "fluc_sz_$1_gap_$2_$3.pop"; }   # sz gap rep
pop_file() { echo "${POP_DIR}/$(pop_name $1 $2 $(( ($3 - 1) % ${N_POP_REPS:-1} + 1 )))"; }   # sz gap rep (invasion)

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
# COMBO_LIST holds "<sz> <gap> <n> <psi> <reg>" per combination. From a grid it is built in the
# order sz > gap > n > psi > reg (reg varies fastest).
COMBO_LIST=()
if (( ${#COMBOS[@]} )); then
    COMBO_LIST=("${COMBOS[@]}")
else
    for _sz in "${SZS[@]}"; do for _gap in "${GAPS[@]}"; do for _n in "${NS[@]}"; do for _psi in "${PSIS[@]}"; do for _reg in "${REGS[@]}"; do
        COMBO_LIST+=("$_sz $_gap $_n $_psi $_reg")
    done; done; done; done; done
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
#   task_params <task>  sets combo_i sz gap n psi reg, rep_first rep_last, combo (the name used for
#                       summary/, logs/ and raw/), part_prefix (parquet part file prefix)
task_params() {
    local t=$1
    combo_i=$(( t % N_COMBOS ))
    local chunk=$(( t / N_COMBOS ))
    read -r sz gap n psi reg <<< "${COMBO_LIST[$combo_i]}"
    rep_first=$(( chunk * RPT + 1 ))
    rep_last=$(( rep_first + RPT - 1 ))
    (( rep_last > REPS )) && rep_last=$REPS
    if [[ $RUN_KIND == burnin ]]; then
        combo="burnin_sz_${sz}_gap_${gap}"
    else
        combo="sz_${sz}_gap_${gap}_reg_${reg}_n_${n}_psi_${psi}_inv_${INVASION_IND}"
    fi
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
# 02-07 were also too slow. Sim times in fluc_env_adjusted_shift_size (from scratch to tick
# 11000): ~31 min per sim on cbsubscb16; 10 sims took 4.5-10.4 h per task depending on the node.
EXCLUDE_NODES='cbsubscb[02-04,06-12]'

# tsv output folders converted to parquet; log/ holds .txt and is only archived
# (a burnin run converts nothing: its output is the pop files, plus SLiM logs)
if [[ $RUN_KIND == burnin ]]; then
    TSV_TYPES=()
else
    TSV_TYPES=(phenotype genotype data allele group_tag social metadata)
fi
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
