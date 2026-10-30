#!/bin/bash
# Run this from the cluster login node (cbsulogin), inside this directory under /home.
# Submits one array job, one task per (shift speed, n, psi, reg) combination (slurm_job.sh).
# Each task runs all its replicates, converts them to parquet and archives the raw files in
# local /workdir, then copies the results to $OUT on Messer lab storage:
#   $OUT/parquet/<type>/shift_speed=../reg=../n=../part-psi_<psi>-0.parquet
#   $OUT/raw/<combo>.tar.zst              raw tsv + SLiM logs
#   $OUT/summary/<combo>/                 conversion summaries, failed_sims.txt, _DONE when all ok
# Settings are in pipeline_config.sh, including the nodes to avoid (EXCLUDE_NODES) and the
# time limit (TIME_LIMIT), which are passed to sbatch.
#
#   bash submit_slurm.sh --test   2 combos x 2 reps into <lab storage>/test_run
#   bash submit_slurm.sh          the full run

set -euo pipefail

if [[ ${1:-} == "--test" ]]; then
    # exported, so sbatch passes them on to every job
    export DIR_NAME=test_run REPS=2
    N_SIM_TASKS=2
fi
source pipeline_config.sh
mount_lab_storage > /dev/null

# ---- checks before submitting anything ----
fail=0
check() {  # check "description" command...
    local msg=$1; shift
    local out
    if out=$("$@" 2>&1); then echo "  ok    $msg"
    else echo "  FAIL  $msg"; echo "$out" | tail -5 | sed 's/^/          /'; fail=1; fi
}
echo "Checking setup:"
check "on a /home path (cluster nodes can see it): $PWD" test "${PWD#/home/}" != "$PWD"
check "lab storage writable: $(dirname $OUT)" test -w "$(dirname $OUT)"
check "SLURM account bscb10" bash -c "sacctmgr -n show assoc user=$USER format=account | grep -qw bscb10"
check "SLiM runs: $HOME/bin/slim5.0" $HOME/bin/slim5.0 -v
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
echo "  pop files: $(ls ${POP_DIR} | wc -l) (expect 19844)"
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
if [[ -d ${OUT}/parquet || -d ${OUT}/raw ]]; then
    echo "WARNING: ${OUT} already has parquet/ or raw/; combos that run again overwrite their own files."
fi

mkdir -p ${OUT} /home/$USER/slurm_logs

$PY make_model.py
# for record keeping
cp ave.slim ${OUT}/moving_ave.slim
cp ext.slim ${OUT}/moving_ext.slim
cp fit.slim ${OUT}/moving_fit.slim

n_tasks=${N_SIM_TASKS:-$N_COMBOS}
sim_id=$(sbatch --parsable --time=$TIME_LIMIT --exclude="$EXCLUDE_NODES" --array=0-$(( n_tasks - 1 )) slurm_job.sh)

echo "sims:    ${sim_id}  (${n_tasks} tasks x ${REPS} reps, --time=$TIME_LIMIT --exclude=$EXCLUDE_NODES)"
echo "output:  ${OUT}"
echo "finished combos: ls ${OUT}/summary/*/_DONE | wc -l   (expect ${n_tasks})"
