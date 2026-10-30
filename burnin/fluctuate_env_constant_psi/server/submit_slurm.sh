#!/bin/bash
# Run this from the cluster login node (cbsulogin), inside this directory under /home.
# Submits the fluctuating-environment sims (slurm_job.sh): one task per (gap, sz, n, psi, reg)
# combination and chunk of REPS_PER_TASK replicates, split over several arrays of at most
# MAX_ARRAY tasks. Each task runs its replicates, converts them to parquet and archives the
# SLiM logs in local /workdir (raw tsv are not kept), then copies the results to $OUT on
# Messer lab storage:
#   $OUT/parquet/<type>/fluc_opt=../fluc_interval=../reg=../n=../part-*.parquet
#   $OUT/logs/<chunk>_logs.tar.zst         SLiM logs
#   $OUT/summary/<chunk>/                  conversion summaries, failed_sims.txt, _DONE when all ok
# Settings are in pipeline_config.sh, including the nodes to avoid (EXCLUDE_NODES) and the
# time limit (TIME_LIMIT), which are passed to every sbatch call.
#
#   bash submit_slurm.sh --test   2 tasks x 2 reps into <lab storage>/test_run
#   bash submit_slurm.sh          the full run

set -euo pipefail

if [[ ${1:-} == "--test" ]]; then
    # exported, so sbatch passes them on to every job
    export DIR_NAME=test_run REPS=2 REPS_PER_TASK=2
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
check "SLiM runs: bin/slim5.0" bin/slim5.0 -v
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
if [[ -d ${OUT}/parquet || -d ${OUT}/logs ]]; then
    echo "WARNING: ${OUT} already has parquet/ or logs/; tasks that run again overwrite their own files."
fi

mkdir -p ${OUT} /home/$USER/slurm_logs

$PY make_model.py
# for record keeping
cp ave.slim ${OUT}/fluc_ave.slim
cp ext.slim ${OUT}/fluc_ext.slim
cp fit.slim ${OUT}/fluc_fit.slim

n_tasks=${N_SIM_TASKS:-$N_TASKS}
for (( off=0; off<n_tasks; off+=MAX_ARRAY )); do
    count=$(( n_tasks - off < MAX_ARRAY ? n_tasks - off : MAX_ARRAY ))
    id=$(sbatch --parsable --time=$TIME_LIMIT --exclude="$EXCLUDE_NODES" \
        --array=0-$(( count - 1 )) --export=ALL,TASK_OFFSET=${off} slurm_job.sh)
    echo "sims:    ${id}  (tasks ${off}-$(( off + count - 1 )))"
done
echo "${n_tasks} tasks x ${REPS_PER_TASK} reps, --time=$TIME_LIMIT --exclude=$EXCLUDE_NODES, output: ${OUT}"
echo "finished tasks: ls ${OUT}/summary/*/_DONE | wc -l   (expect ${n_tasks})"
