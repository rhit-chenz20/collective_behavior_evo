#!/bin/bash
# Run the full grid on the Mac (the cluster version is in ../server). Keep this file working with macOS bash 3.2.
# One worker = one (shift speed, n, psi, reg) combination, running its REPS replicates and
# then converting and archiving them (run_combo.sh); up to N_JOBS workers run at once.
# Settings are in pipeline_config.sh. Results go to $OUT, as on the cluster:
#   $OUT/parquet/<type>/shift_speed=../reg=../n=../part-psi_<psi>-0.parquet
#   $OUT/raw/<combo>.tar.zst              raw tsv + SLiM logs
#   $OUT/summary/<combo>/                 conversion summaries, failed_sims.txt, _DONE when all ok
#   $OUT/task_logs/task_<index>.log       a worker's own log
#   $OUT/failed_combos.txt                combos without _DONE at the end of the run
#
# Rerunning skips every combo that has a _DONE, so a stopped run or failed combos only need
# the same command again. Posts to ntfy.sh/$NTFY_TOPIC when the run ends or is stopped, and
# at the first failed combo.
#
#   bash run_sim.sh --test    2 combos x 2 reps into <data>/test_run
#   nohup bash run_sim.sh > run_sim.log 2>&1 &    the full run; survives closing the terminal
# Needs on the Mac: SLiM 5 (slim), python3 with polars, zstd, and the pop files in $POP_DIR.

set -uo pipefail
cd "$(dirname "$0")"

if [[ ${1:-} == "--test" ]]; then
    # exported, so the workers get them too
    export DIR_NAME=test_run REPS=2
    N_SIM_TASKS=2
fi
source pipeline_config.sh

notify() {
    curl -s -o /dev/null --max-time 30 -d "$1" ntfy.sh/${NTFY_TOPIC} || echo "WARNING: could not send: $1" >&2
}

# ---- one worker (started by xargs below) ----
if [[ ${1:-} == "--task" ]]; then
    if ! bash run_combo.sh $2 > ${OUT}/task_logs/task_$2.log 2>&1; then
        # mkdir is atomic, so only the first failed worker sends this
        mkdir ${OUT}/.failure_notified 2> /dev/null &&
            notify "first failure: task $2, see ${OUT}/task_logs/task_$2.log; run continues"
    fi
    exit 0
fi

# ---- checks before starting anything ----
fail=0
check() {  # check "description" command...
    local msg=$1; shift
    local out
    if out=$("$@" 2>&1); then echo "  ok    $msg"
    else echo "  FAIL  $msg"; echo "$out" | tail -5 | sed 's/^/          /'; fail=1; fi
}
mkdir -p ${OUT}/task_logs
echo "Checking setup:"
check "no spaces in paths: $PWD $OUT" test "$PWD$OUT" = "${PWD// /}${OUT// /}"
check "output writable: $OUT" test -w "$OUT"
check "SLiM 5 runs: $SLIM" bash -c "$SLIM -v | grep -q 'version 5'"
check "pop files in ${POP_DIR}" test -f ${POP_DIR}/cons_n_2_psi_0.0_reg_ave_1.pop
check "python with polars: $PY" $PY -c "import polars"
check "zstd: $ZSTD" $ZSTD --version
check "rsync" which rsync
check "curl" which curl
for f in base_script.slim make_model.py run_combo.sh convert_parquet.py; do
    check "file $f" test -f $f
done
if (( fail )); then
    echo "Setup checks failed; nothing started."
    exit 1
fi
# every (n, psi, reg) needs a pop file for each of its REPS replicates
n_missing=0
for n in "${NS[@]}"; do for psi in "${PSIS[@]}"; do for reg in "${REGS[@]}"; do
    for (( rep=1; rep<=REPS; rep++ )); do
        [[ -f ${POP_DIR}/cons_n_${n}_psi_${psi}_reg_${reg}_${rep}.pop ]] || (( ++n_missing ))
    done
done; done; done
if (( n_missing )); then
    echo "  WARNING: ${n_missing} pop files missing for reps 1-${REPS}; those sims will be listed as failed"
fi

# ---- the run ----
$PY make_model.py || exit 1
# for record keeping
cp ave.slim ${OUT}/moving_ave.slim
cp ext.slim ${OUT}/moving_ext.slim
cp fit.slim ${OUT}/moving_fit.slim

n_tasks=${N_SIM_TASKS:-$N_COMBOS}
list_unfinished() {  # "<task index> <combo>" per combo without a _DONE, same order as run_combo.sh
    local i=0 sp n psi reg combo
    for sp in "${SPEEDS[@]}"; do for n in "${NS[@]}"; do for psi in "${PSIS[@]}"; do for reg in "${REGS[@]}"; do
        (( i < n_tasks )) || return 0
        combo="sp_${sp}_reg_${reg}_n_${n}_psi_${psi}"
        [[ -f ${OUT}/summary/${combo}/_DONE ]] || echo "$i $combo"
        (( ++i ))
    done; done; done; done
}

# leftovers from a run that was cut off (e.g. power loss) would mix into the rerun
rm -rf ${WORK_ROOT} ${OUT}/.failure_notified

todo=$(list_unfinished | cut -d' ' -f1)
n_todo=$(echo "$todo" | grep -c .)
echo "$(date) ${n_todo} of ${n_tasks} combos to run, ${REPS} reps each, ${N_JOBS} at a time"
echo "output:  ${OUT}"

# TERM or Ctrl-C: stop the workers too (kill 0 reaches the whole process group). Each
# worker removes its scratch folder, and the rerun continues from the finished combos.
trap 'trap - INT TERM; notify "run stopped before finishing, rerun run_sim.sh to continue"; kill 0' INT TERM

# keep the Mac from idle sleeping while this script runs
command -v caffeinate > /dev/null && { caffeinate -i -w $$ & }

if (( n_todo )); then
    echo "$todo" | xargs -P ${N_JOBS} -n 1 bash run_sim.sh --task &
    wait $!
fi

# ---- final report ----
n_failed=$(list_unfinished | tee ${OUT}/failed_combos.txt | wc -l | tr -d ' ')
echo "$(date) finished, ${n_failed} of ${n_tasks} combos failed"
rm -rf ${WORK_ROOT}
if (( n_failed )); then
    notify "job done, but ${n_failed} of ${n_tasks} combos failed, see ${OUT}/failed_combos.txt"
    exit 1
fi
notify "job done"
