#!/bin/bash
# Run this from the cluster login node (cbsulogin); this package must be under /home.
# Invasion of social learners into a psi = 0 population under a moving optimum.
# Submits one array job for a run defined in runs/<run>.sh (its combinations, number of
# invaders, invasion tick, end time, reps and output folder): one task per combination, or per chunk of its reps when the run splits
# them (REPS_PER_TASK). Each task runs its replicates, converts them to parquet and archives
# the raw files in local /workdir, then copies everything to $OUT on Messer lab storage:
#   $OUT/parquet/<type>/shift_speed=../invasion_ind=../reg=../n=../part-psi_<psi>-[reps_<a>-<b>-]0.parquet
#   $OUT/raw/<combo>.tar.zst              raw tsv + SLiM logs
#   $OUT/summary/<combo>/                 conversion summaries, failed_sims.txt,
#                                         (<combo> ends in _reps_<a>-<b> when reps are split)
#                                         sim_seconds.tsv, _DONE when all ok
#   $OUT/logs/<combo>.<jobid>.log         the task's own log (+ .convert.log)
#   $OUT/logs/submit_<time>.log           this script's output
#   $OUT/code/ or $OUT/code_<run>/        the exact code the tasks ran (code_<run>/ when the
#                                         run adds tasks to another run's folder)
#   $OUT/summary/<combo>/_RUN             the run a task belongs to; a run never overwrites
#                                         another run's tasks in a shared folder
# Settings are in pipeline_config.sh and runs/<run>.sh; nothing is taken from your shell
# environment. Runs: moving_inv10, moving_inv500 (both into invasion_moving).
#
#   bash submit_slurm.sh <run> test     the run's first 3 combinations x 2 reps
#   bash submit_slurm.sh <run> pilot    every combination x 3 reps
#   bash submit_slurm.sh <run> full     every combination x N_REPS_FULL reps
#   bash submit_slurm.sh <run> <test|pilot|full> --status
#                                 how many of the run's tasks are done; submits nothing
#   bash submit_slurm.sh <run> <pilot|full> --unfinished
#                                 only the tasks without a _DONE that are not still queued
#                                 or running; refused if the code changed since that run
# pilot and full are refused until the run's test has finished with the same code.

set -euo pipefail

# Slurm emails this address once per submitted array: when the whole array has ended (END) and
# if any task failed (FAIL); the subject gives the job id, name and state. Kept here rather
# than in pipeline_config.sh so changing it does not change the code fingerprint.
MAIL_USER=zc524@cornell.edu

usage="usage: bash submit_slurm.sh <run> <test|pilot|full> [--unfinished | --status]   (runs: $(cd "$(dirname "$(readlink -f "$0")")/runs" && ls *.sh | sed 's/\.sh$//' | tr '\n' ' '))"
run=${1:-}
mode=${2:-}
case $mode in test|pilot|full) ;; *) echo "$usage" >&2; exit 2 ;; esac
unfinished=0 status=0
if (( $# == 3 )) && [[ $3 == --unfinished && $mode != test ]]; then
    unfinished=1
elif (( $# == 3 )) && [[ $3 == --status ]]; then
    status=1
elif (( $# != 2 )); then
    echo "$usage  (--unfinished works with pilot or full only)" >&2; exit 2
fi

code_dir=$(dirname "$(readlink -f "$0")")
cd "$code_dir"
source ./pipeline_config.sh "$run" "$mode"

# a nonzero exit can mean "already mounted"; the writable check below is what counts
mount_lab_storage > /dev/null || true

if (( status )); then
    n_done=0 n_other=0 pending=()
    for (( t=0; t<N_TASKS; t++ )); do
        task_params $t
        owner=$(combo_owner $combo)
        if [[ -n $owner && $owner != "$run" ]]; then (( ++n_other )); pending+=("$combo (belongs to run $owner)")
        elif [[ -f ${OUT}/summary/${combo}/_DONE ]]; then (( ++n_done ))
        else pending+=("$combo"); fi
    done
    echo "${run} ${mode}: ${n_done} of ${N_TASKS} tasks done in ${OUT}"
    (( ${#pending[@]} )) && printf '  not done: %s\n' "${pending[@]}" | head -20
    (( ${#pending[@]} > 20 )) && echo "  ... $(( ${#pending[@]} - 20 )) more"
    squeue -u "$ME" -h -n "$JOB_NAME" -o '%i %T %M %N' 2>/dev/null | sed 's/^/  queued: /' | head -10
    exit 0
fi
stamp=$(date +%Y%m%d_%H%M%S)
mkdir -p "$OUT/logs"
exec > >(tee -a "$OUT/logs/submit_${stamp}.log") 2>&1
echo "submit_slurm.sh $*  $(date)  from $code_dir"
echo "run ${run} (${RUN_FILE}): ${N_COMBOS} combinations, ${INVASION_IND} invaders at tick ${INVASION_TICK}, to tick ${END_TIME} at most, ${mode}: ${REPS} reps, ${RPT} per task -> ${N_TASKS} tasks, IDs from $(sim_id 0 1)"

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
check "pop folder on a /home path: $POP_DIR" test "${POP_DIR#/home/}" != "$POP_DIR"
check "lab storage writable: $OUT" test -w "$OUT"
check "SLURM account bscb10" bash -c "sacctmgr -n show assoc user=$ME format=account | grep -qw bscb10"
check "SLiM 5 runs: $SLIM" bash -c "$SLIM -v | grep -q 'version 5'"
# only check polars is installed: importing it fails on the login node, whose CPU lacks
# instructions polars is built for (AVX2 etc.); the compute nodes that run the jobs have them
check "python env (polars installed) at $PY" $PY -c "import importlib.util, sys; sys.exit(importlib.util.find_spec('polars') is None)"
check "zstd at $ZSTD" $ZSTD --version
check "rsync" which rsync
check "at most $MAX_ARRAY tasks in one array (${N_TASKS})" test $N_TASKS -le $MAX_ARRAY
for f in base_script.slim make_model.py slurm_job.sh convert_parquet.py; do
    check "file $f" test -f $f
done
check "END_TIME after the shift (10000): $END_TIME" test "$END_TIME" -gt 10000
check "replicates >= 1: $REPS" test "$REPS" -ge 1
check "1 <= INVASION_IND < 1000 (N): ${INVASION_IND:-}" test "${INVASION_IND:-0}" -ge 1 -a "${INVASION_IND:-0}" -lt 1000
check "10000 <= INVASION_TICK < END_TIME: ${INVASION_TICK:-}" test "${INVASION_TICK:-0}" -ge 10000 -a "${INVASION_TICK:-0}" -lt "$END_TIME"
check "ID_OFFSET set: ${ID_OFFSET:-}" test -n "${ID_OFFSET:-}"
bad_combos=$(for c in "${COMBO_LIST[@]}"; do read -r a b c2 d <<< "$c"
    [[ $a == *.* && $c2 == *.* && $b =~ ^[0-9]+$ && $d =~ ^(ave|ext|fit)$ ]] || echo "'$c'"; done)
check "combinations are '<speed> <n> <psi> <reg>' with decimal speed and psi (1.1, 0.0)" test -z "$bad_combos"
[[ -n $bad_combos ]] && echo "          bad: $bad_combos"
dup_combos=$(printf '%s\n' "${COMBO_LIST[@]}" | sort | uniq -d)
check "no combination listed twice" test -z "$dup_combos"
# a folder can hold several runs: refuse if any of this run's tasks already belong to another
others=()
for (( t=0; t<N_TASKS; t++ )); do
    task_params $t
    owner=$(combo_owner $combo)
    [[ -n $owner && $owner != "$run" ]] && others+=("$combo ($owner)")
done
check "none of the ${N_TASKS} tasks belongs to another run in ${OUT}" test ${#others[@]} -eq 0
(( ${#others[@]} )) && printf '          %s\n' "${others[@]:0:10}"
if (( fail )); then
    echo "Setup checks failed; nothing submitted."
    exit 1
fi

# every task needs a pop file for each of its replicates
n_missing=0 n_needed=0
for (( t=0; t<N_TASKS; t++ )); do
    task_params $t
    for (( rep=rep_first; rep<=rep_last; rep++ )); do
        (( ++n_needed ))
        [[ -f $(pop_file $n $reg $rep) ]] || (( ++n_missing ))
    done
done
echo "  pop files in ${POP_DIR}: $(( n_needed - n_missing )) of ${n_needed} needed"
if (( n_missing )); then
    echo "  WARNING: ${n_missing} pop files missing; those sims will be listed as failed"
fi

fingerprint=$(code_fingerprint)

# ---- pilot / full: only after a finished test of the same code ----
if [[ $mode != test ]]; then
    # the test's output folder and its task names, from the test settings of the same run
    test_info=$(source ./pipeline_config.sh "$run" test && echo "$OUT $CODE_OUT" &&
                for (( t=0; t<N_TASKS; t++ )); do task_params $t; echo "$combo"; done)
    read -r test_out test_code < <(head -1 <<< "$test_info")
    mapfile -t test_combos < <(tail -n +2 <<< "$test_info")
    echo "Checking the test run in ${test_out}:"
    test_ok=1
    if [[ $(cat "$test_code/code.md5" 2>/dev/null) != "$fingerprint" ]]; then
        echo "  FAIL  the code changed since the test run (or there is none): run 'bash submit_slurm.sh test' first"
        test_ok=0
    fi
    # a _DONE only counts if it was written after the test's code was submitted (code.md5):
    # an older one is left from an earlier test and is replaced only when its task finishes
    for combo in "${test_combos[@]}"; do
        done_f=$test_out/summary/${combo}/_DONE
        if [[ -f $done_f && $done_f -nt $test_code/code.md5 ]]; then echo "  ok    ${combo}"
        elif [[ -f $done_f ]]; then echo "  FAIL  ${combo}: its _DONE is from an earlier test (the latest is still running or failed)"; test_ok=0
        else echo "  FAIL  ${combo} has no _DONE (still running, failed, or not run)"; test_ok=0; fi
    done
    if (( ! test_ok )); then
        echo "Test run not finished with the current code; nothing submitted."
        exit 1
    fi
    # project the slowest task from the test's slowest sim
    max_s=$(cat "$test_out"/summary/*/sim_seconds.tsv | awk -F'\t' '$1!="rep" && $4>m {m=$4} END{print m+0}')
    proj_s=$(( max_s * RPT + 1800 ))   # + 30 min for parquet, archive, copy
    echo "  slowest test sim: ${max_s} s -> ${RPT} reps per task ~ $(( proj_s / 60 )) min; limit ${TIME_LIMIT}"
    if (( proj_s > TIME_LIMIT_H * 3600 )); then
        echo "  FAIL  that is over the ${mode} time limit; raise TIME_LIMIT_H_${mode^^} or lower REPS_PER_TASK in runs/${run}.sh (a time-limit change needs no re-test)"
        exit 1
    fi
fi

# ---- tasks to submit ----
if (( unfinished )); then
    if [[ $(cat "$CODE_OUT/code.md5" 2>/dev/null) != "$fingerprint" ]]; then
        echo "The code differs from the ${mode} run already in ${OUT}; rerunning would mix results of two versions. Nothing submitted."
        exit 1
    fi
    # tasks of this mode still pending or running (array indices = task numbers)
    declare -A queued=()
    while read -r k; do [[ -n $k ]] && queued[$k]=1; done \
        < <(squeue -u "$ME" -h -r -n "$JOB_NAME" -o %K)
    todo=() n_queued=0
    for (( t=0; t<N_TASKS; t++ )); do
        task_params $t
        [[ -f ${OUT}/summary/${combo}/_DONE ]] && continue
        if [[ -n ${queued[$t]:-} ]]; then (( ++n_queued )); continue; fi
        todo+=($t)
    done
    echo "  ${#todo[@]} tasks to rerun; ${n_queued} without _DONE are still queued or running (left alone)"
    if (( ${#todo[@]} == 0 )); then
        echo "Nothing to submit."
        exit 0
    fi
    array=$(IFS=,; echo "${todo[*]}")
    n_submit=${#todo[@]}
else
    n_own_done=0
    for (( t=0; t<N_TASKS; t++ )); do task_params $t; [[ -f ${OUT}/summary/${combo}/_DONE ]] && (( ++n_own_done )); done
    if (( n_own_done )); then
        echo "WARNING: ${n_own_done} of this run's tasks in ${OUT} are already finished; this submit reruns and overwrites them."
        [[ $mode != test ]] && echo "         (bash submit_slurm.sh ${run} ${mode} --unfinished submits only the tasks without a _DONE)"
    fi
    if [[ $mode == test ]]; then
        if squeue -u "$ME" -h -n "$JOB_NAME" | grep -q .; then
            echo "A test job (${JOB_NAME}) is still queued or running; cancel it or wait. Nothing submitted."
            exit 1
        fi
    fi
    array=0-$(( N_TASKS - 1 ))
    n_submit=$N_TASKS
fi

if [[ $mode != test ]]; then
    read -r -p "Submit ${n_submit} ${run} ${mode} tasks (${RPT} reps each, to tick ${END_TIME}) to ${OUT}? Type yes: " answer < /dev/tty
    [[ $answer == yes ]] || { echo "Nothing submitted."; exit 1; }
fi

# ---- frozen code snapshot: the tasks run this copy, never code/ itself ----
snap=${CODE_ROOT}/${DIR_NAME}_${stamp}
mkdir -p "$snap" "$LOG_DIR"
cp -r "$code_dir" "$snap/code"
ln -s "$POP_DIR" "$snap/pop"
cd "$snap/code"
$PY make_model.py
echo "$fingerprint" > code.md5
mkdir -p "$CODE_OUT"
rsync -a --delete "$snap/code/" "$CODE_OUT/"
echo "code:    $snap/code  (copy in $CODE_OUT)"

# --export=NONE: the jobs get no variables from this shell; the mode is their argument
job_id=$(sbatch --parsable --export=NONE --job-name="$JOB_NAME" --time=$TIME_LIMIT \
    --exclude="$EXCLUDE_NODES" --mail-type=END,FAIL --mail-user="$MAIL_USER" \
    --array=$array slurm_job.sh "$run" "$mode")

echo "sims:    ${job_id}  (${n_submit} tasks x up to ${RPT} reps, to tick ${END_TIME}, --time=$TIME_LIMIT --exclude=$EXCLUDE_NODES)"
echo "output:  ${OUT}"
echo "email:   ${MAIL_USER} when the array ends (and if a task fails)"
echo "progress: bash submit_slurm.sh ${run} ${mode} --status"
[[ $mode == test ]] && echo "when all ${N_TASKS} are done: bash submit_slurm.sh ${run} pilot   (or ${run} full)"
exit 0
