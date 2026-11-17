#!/bin/bash
# Run this from the cluster login node (cbsulogin); this folder must be under /home.
# Submits one array task per (shift speed, n, psi, reg) combination (slurm_job.sh). Each task
# runs all its replicates, converts them to parquet and archives the raw files in local
# /workdir, then copies everything to $OUT on Messer lab storage:
#   $OUT/parquet/<type>/shift_speed=../reg=../n=../part-psi_<psi>-0.parquet
#   $OUT/raw/<combo>.tar.zst              raw tsv + SLiM logs
#   $OUT/summary/<combo>/                 conversion summaries, failed_sims.txt,
#                                         sim_seconds.tsv, _DONE when all ok
#   $OUT/logs/<combo>.<jobid>.log         the task's own log (+ .convert.log)
#   $OUT/logs/submit_<time>.log           this script's output
#   $OUT/code/                            the exact code the tasks ran
# Settings are in pipeline_config.sh; nothing is taken from your shell environment.
#
#   bash submit_slurm.sh test   3 tasks (speed 0.01, n 2, psi 0.0, each reg) x 2 reps
#                               into <lab storage>/cons_env_moving_opt_test
#   bash submit_slurm.sh full   all 3300 tasks x 60 reps into <lab storage>/cons_env_moving_opt;
#                               refused until the test run has finished with the same code

set -euo pipefail

if (( $# != 1 )) || [[ $1 != test && $1 != full ]]; then
    echo "usage: bash submit_slurm.sh <test|full>" >&2
    exit 2
fi
mode=$1
here=$(dirname "$(readlink -f "$0")")
cd "$here"
source ./pipeline_config.sh "$mode"

# a nonzero exit can mean "already mounted"; the writable check below is what counts
mount_lab_storage > /dev/null || true
stamp=$(date +%Y%m%d_%H%M%S)
mkdir -p "$OUT/logs"
exec > >(tee -a "$OUT/logs/submit_${stamp}.log") 2>&1
echo "submit_slurm.sh $mode  $(date)  from $here"

# ---- checks before submitting anything ----
fail=0
check() {  # check "description" command...
    local msg=$1; shift
    local out
    if out=$("$@" 2>&1); then echo "  ok    $msg"
    else echo "  FAIL  $msg"; echo "$out" | tail -5 | sed 's/^/          /'; fail=1; fi
}
echo "Checking setup:"
check "code snapshots on a /home path (cluster nodes can see it): $CODE_ROOT" test "${CODE_ROOT#/home/}" != "$CODE_ROOT"
check "lab storage writable: $OUT" test -w "$OUT"
check "SLURM account bscb10" bash -c "sacctmgr -n show assoc user=$ME format=account | grep -qw bscb10"
check "SLiM runs: bin/slim5.0" bin/slim5.0 -v
check "pop files in ${POP_DIR}" test -f ${POP_DIR}/cons_n_2_psi_0.0_reg_ave_1.pop
# only check polars is installed: importing it fails on the login node, whose CPU lacks
# instructions polars is built for (AVX2 etc.); the compute nodes that run the jobs have them
check "python env (polars installed) at $PY" $PY -c "import importlib.util, sys; sys.exit(importlib.util.find_spec('polars') is None)"
check "zstd at $ZSTD" $ZSTD --version
check "rsync" which rsync
for f in base_script.slim make_model.py slurm_job.sh convert_parquet.py; do
    check "file $f" test -f $f
done
if (( fail )); then
    echo "Setup checks failed; nothing submitted."
    exit 1
fi

# every task needs a pop file for each of its REPS replicates
n_missing=0
for (( t=0; t<N_TASKS; t++ )); do
    task_params $t
    for (( rep=1; rep<=REPS; rep++ )); do
        [[ -f ${POP_DIR}/cons_n_${n}_psi_${psi}_reg_${reg}_${rep}.pop ]] || (( ++n_missing ))
    done
done
echo "  pop files missing for these tasks: ${n_missing}"

fingerprint=$(code_fingerprint)

# ---- full run: only after a finished test run of the same code ----
if [[ $mode == full ]]; then
    read -r test_out test_tasks < <(source ./pipeline_config.sh test && echo "$OUT $N_TASKS")
    echo "Checking the test run in ${test_out}:"
    test_ok=1
    if [[ $(cat "$test_out/code/code.md5" 2>/dev/null) != "$fingerprint" ]]; then
        echo "  FAIL  the code changed since the test run (or there is none): run 'bash submit_slurm.sh test' first"
        test_ok=0
    fi
    for (( t=0; t<test_tasks; t++ )); do
        task_params $t
        if [[ -f $test_out/summary/${combo}/_DONE ]]; then echo "  ok    ${combo}"
        else echo "  FAIL  ${combo} has no _DONE (still running, failed, or not run)"; test_ok=0; fi
    done
    if (( ! test_ok )); then
        echo "Test run not finished with the current code; nothing submitted."
        exit 1
    fi
    # project the slowest task from the test's slowest sim
    max_s=$(cat "$test_out"/summary/*/sim_seconds.tsv | awk -F'\t' '$1!="rep" && $4>m {m=$4} END{print m+0}')
    proj_s=$(( max_s * REPS + 1800 ))   # + 30 min for parquet, archive, copy
    echo "  slowest test sim: ${max_s} s -> ${REPS} reps ~ $(( proj_s / 60 )) min; limit ${TIME_LIMIT}"
    if (( proj_s > TIME_LIMIT_H * 3600 )); then
        echo "  FAIL  that is over TIME_LIMIT; raise TIME_LIMIT_H in pipeline_config.sh (no re-test needed)"
        exit 1
    fi
    if compgen -G "$OUT/summary/*/_DONE" > /dev/null; then
        echo "WARNING: ${OUT} already has finished tasks; a full submit reruns and overwrites all of them."
    fi
    read -r -p "Submit ${N_TASKS} tasks x ${REPS} reps to ${OUT}? Type yes: " answer < /dev/tty
    [[ $answer == yes ]] || { echo "Nothing submitted."; exit 1; }
fi

# ---- frozen code snapshot: the tasks run this copy, never this folder ----
snap=${CODE_ROOT}/${DIR_NAME}_${stamp}
mkdir -p "$CODE_ROOT" "$LOG_DIR"
cp -rL "$here" "$snap"
cd "$snap"
$PY make_model.py
echo "$fingerprint" > code.md5
rsync -a --delete "$snap/" "$OUT/code/"
echo "code:    $snap  (copy in $OUT/code)"

# --export=NONE: the jobs get no variables from this shell; mode and offset are arguments
for (( off=0; off<N_TASKS; off+=MAX_ARRAY )); do
    count=$(( N_TASKS - off < MAX_ARRAY ? N_TASKS - off : MAX_ARRAY ))
    id=$(sbatch --parsable --export=NONE --time=$TIME_LIMIT --exclude="$EXCLUDE_NODES" \
        --array=0-$(( count - 1 )) slurm_job.sh "$mode" "$off")
    echo "sims:    ${id}  (tasks ${off}-$(( off + count - 1 )))"
done
echo "${N_TASKS} tasks x ${REPS} reps, --time=$TIME_LIMIT --exclude=$EXCLUDE_NODES"
echo "output:  ${OUT}"
echo "finished tasks: ls ${OUT}/summary/*/_DONE | wc -l   (expect ${N_TASKS})"
if [[ $mode == test ]]; then
    echo "when all ${N_TASKS} are done: bash submit_slurm.sh full"
fi
