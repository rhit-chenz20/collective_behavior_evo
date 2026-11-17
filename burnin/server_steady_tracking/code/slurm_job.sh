#!/bin/bash -l
#SBATCH --account=bscb10
#SBATCH --partition=regular,long7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --job-name=steady_tracking
#SBATCH --output=/home/%u/slurm_logs/steady_tracking_%A_%a.out

# One array task = one (shift speed, n, psi, reg) combination, or one chunk of its replicates
# when the run splits them (REPS_PER_TASK); task order in pipeline_config.sh.
#   sbatch [...] slurm_job.sh <run> <test|pilot|full>
# Submit only via submit_slurm.sh, which runs it from a frozen code snapshot with
# --export=NONE: every setting comes from pipeline_config.sh, runs/<run>.sh and the mode.
# Every sim runs to the run's END_TIME (passed to SLiM as end_time); failing the fitness check
# and a steady lag are only recorded in the metadata; see base_script.slim. --time and
# --job-name above are defaults: submit_slurm.sh passes the mode's TIME_LIMIT, EXCLUDE_NODES
# and job name.
# Everything happens in local /workdir. Inputs (the code snapshot, SLiM and the combo's pop
# files) are copied once at the start:
#   1. run the task's replicates one after another   -> $WD/raw/<type>/*.tsv, $WD/raw/log/*.txt
#   2. convert each type to parquet                    -> $WD/parquet/<type>/shift_speed=../reg=../n=../part-<part_prefix>0.parquet
#      (part_prefix: psi_<psi>- , or psi_<psi>-reps_<first>-<last>- when reps are split)
#   3. archive the raw files and check the archive     -> $WD/<combo>.tar.zst
#   4. copy parquet, archive and summaries to $OUT     -> $OUT/{parquet,raw,summary}
# $OUT/summary/<combo>/_DONE is written last, only if every step succeeded.
# All output of this script goes to $OUT/logs/<combo>.<jobid>.log, copied at exit for any
# reason; the SLURM file in /home/<user>/slurm_logs is then deleted (it is only kept if
# the job died before it could copy its log, e.g. a broken mount).

set -u
run=${1:?usage: sbatch slurm_job.sh <run> <test|pilot|full>}
mode=${2:?usage: sbatch slurm_job.sh <run> <test|pilot|full>}
source ${SLURM_SUBMIT_DIR}/pipeline_config.sh "$run" "$mode" || exit 1

# ---- this task's parameters ----
task=$SLURM_ARRAY_TASK_ID
task_params $task
slurm_log=${LOG_DIR}/steady_tracking_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out

# local scratch for this job; removed when the job ends for any reason (done, error,
# scancel or time limit: SLURM sends TERM first, which the trap turns into a normal exit)
WD=/workdir/${ME}/${SLURM_JOB_ID}
mkdir -p $WD/raw $WD/pop $WD/logs || exit 1
job_log=$WD/logs/${combo}.${SLURM_JOB_ID}.log
echo "task ${task} (${combo}): log goes to ${OUT}/logs/$(basename $job_log)"
exec > $job_log 2>&1

save_logs() {
    local status=$?
    echo "$(date '+%F %T') job exit status ${status}"
    [[ -f $WD/convert.log ]] && cp $WD/convert.log $WD/logs/${combo}.${SLURM_JOB_ID}.convert.log
    if mkdir -p $OUT/logs && rsync -a $WD/logs/ $OUT/logs/; then
        rm -f $slurm_log
    else
        # lab storage unreachable: keep the logs on /home instead
        mkdir -p $LOG_DIR && cp $WD/logs/* $LOG_DIR/
    fi
    cd /; rm -rf $WD
}
trap save_logs EXIT
trap 'echo "$(date "+%F %T") received TERM: time limit (${TIME_LIMIT}) or scancel"; exit 143' TERM

echo "task ${task} (array ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}, run ${run}, mode ${mode}): ${combo}, reps ${rep_first}-${rep_last}, to tick ${END_TIME}, on $(hostname), $(date '+%F %T')"
mount_lab_storage

# polars (parquet conversion) needs CPU instructions some nodes lack (see EXCLUDE_NODES);
# fail now rather than after hours of sims
if ! $PY -c "import polars" > /dev/null 2>&1; then
    echo "ERROR: polars cannot run on $(hostname); add it to EXCLUDE_NODES in pipeline_config.sh and rerun this task" >&2
    exit 1
fi

# read inputs once, up front: the code snapshot (the scripts), the SLiM binary and only this
# combo's pop files (one per rep) from $POP_DIR; nothing is read from them during the sims.
# The pop file depends only on n, psi, reg and rep, so all shift speeds of a rep start from
# the same burnin
cp -r ${SLURM_SUBMIT_DIR} $WD/code
cp -L $SLIM $WD/slim
SLIM=$WD/slim
for (( rep=rep_first; rep<=rep_last; rep++ )); do
    f=${POP_DIR}/cons_n_${n}_psi_${psi}_reg_${reg}_${rep}.pop
    [[ -f $f ]] && cp $f $WD/pop/
done
echo "copied inputs: $(ls $WD/pop | wc -l) pop files, $(du -sh $WD/code | cut -f1) code"

# never overwrite a task of another run sharing this folder (submit_slurm.sh checks this too)
owner=$(combo_owner $combo)
if [[ -n $owner && $owner != "$run" ]]; then
    echo "ERROR: ${OUT}/summary/${combo} belongs to run '${owner}', not '${run}'; nothing run" >&2
    exit 1
fi

for d in "${ARCHIVE_TYPES[@]}"; do mkdir -p $WD/raw/$d; done
cd $WD/raw   # SLiM writes its outputs relative to here

# ---- 1. simulations ----
failed_sims=()
printf 'rep\tID\tstatus\tseconds\n' > $WD/sim_seconds.tsv
for (( rep=rep_first; rep<=rep_last; rep++ )); do
    # unique sim ID (sim_id in pipeline_config.sh): rep > combination, plus the run's ID_OFFSET
    i=$(sim_id $combo_i $rep)
    tag="n_${n}_psi_${psi}_sp_${sp}_reg_${reg}_${rep}_${i}"
    POP_FILE=$WD/pop/cons_n_${n}_psi_${psi}_reg_${reg}_${rep}.pop
    echo "$(date +%T) rep=${rep} ID=${i}"
    if [[ ! -f "$POP_FILE" ]]; then
        echo "  missing pop file: ${POP_DIR}/$(basename $POP_FILE)"
        failed_sims+=($tag)
        continue
    fi

    t0=$SECONDS
    $SLIM -d psi=$psi \
        -d ID="${i}" \
        -d rep="${rep}" \
        -d n=$n \
        -d shift_speed=$sp \
        -d end_time=$END_TIME \
        -d "metadata_fn='metadata/${tag}.tsv'" \
        -d "phenotype_fn='phenotype/${tag}.tsv'" \
        -d "genotype_fn='genotype/${tag}.tsv'" \
        -d "fre_output='data/${tag}.tsv'" \
        -d "burnin_fn='${POP_FILE}'" \
        -d "group_tag_fn='group_tag/${tag}.tsv'" \
        -d "allele_fn='allele/${tag}.tsv'" \
        $WD/code/${reg}.slim &> log/${tag}.txt
    status=$?
    secs=$(( SECONDS - t0 ))
    printf '%s\t%s\t%s\t%s\n' $rep $i $status $secs >> $WD/sim_seconds.tsv
    echo "  exit status ${status}, ${secs} s"
    if (( status != 0 )); then
        # drop any partial outputs so only complete sims reach the parquet; keep its log
        failed_sims+=($tag)
        rm -f */${tag}.tsv
    fi
done
cd $WD
max_s=$(awk -F'\t' 'NR>1 && $4>m {m=$4} END{print m+0}' $WD/sim_seconds.tsv)
echo "$(date +%T) sims done, ${#failed_sims[@]} failed; slowest sim ${max_s} s (x ${REPS_PER_TASK} reps per full task = $(( max_s * REPS_PER_TASK / 60 )) min)"

# ---- 2. convert to parquet (one part file per type; one python process for all types) ----
$PY $WD/code/convert_parquet.py --layout moving --raw-dir $WD/raw --type "${TSV_TYPES[@]}" --out-dir $WD/parquet \
    --workers 1 --sims-per-part $(( rep_last - rep_first + 1 )) --part-prefix "$part_prefix" \
    --summary-dir $WD/summary/${combo} > $WD/convert.log 2>&1
# a type succeeded iff its summary was written
failed_types=()
for t in "${TSV_TYPES[@]}"; do
    [[ -f $WD/summary/${combo}/${t}.json ]] || failed_types+=($t)
done
if (( ${#failed_types[@]} )); then
    echo "  convert FAILED for ${failed_types[*]}:"; grep ERROR $WD/convert.log
fi
echo "$(date +%T) convert done, ${#failed_types[@]} types failed"

# ---- 3. archive the raw files and check the archive holds every file ----
arc=$WD/${combo}.tar.zst
n_files=$(find $WD/raw -type f | wc -l)
tar -cf - -C $WD/raw . | $ZSTD -T${SLURM_CPUS_PER_TASK:-1} -10 -q -o $arc
n_archived=$($ZSTD -dc $arc | tar -tf - | grep -vc '/$')
if (( n_archived != n_files )); then
    echo "ERROR: archive has ${n_archived} files, raw has ${n_files}" >&2
    archive_ok=0
else
    archive_ok=1
    echo "$(date +%T) archive ok: ${n_files} files, $(du -h $arc | cut -f1)"
fi

# ---- 4. copy to lab storage (write only) ----
# rsync writes each file under a temporary name and renames it when complete, so a
# killed job never leaves a truncated file under its final name
mkdir -p $WD/summary/${combo}
echo "$run" > $WD/summary/${combo}/_RUN
printf '%s\n' "${failed_sims[@]}" > $WD/summary/${combo}/failed_sims.txt
cp $WD/sim_seconds.tsv $WD/summary/${combo}/
all_ok=1
(( ${#failed_sims[@]} || ${#failed_types[@]} || ! archive_ok )) && all_ok=0

mount_lab_storage
copy_ok=1
mkdir -p $OUT/raw $OUT/summary/${combo} || copy_ok=0
# a rerun first removes this combo's outputs from an earlier run (exact paths, nothing is
# listed or read), so storage never mixes old files with the new run's results
rm -f $OUT/summary/${combo}/* $OUT/raw/${combo}.tar.zst
for t in "${TSV_TYPES[@]}"; do
    rm -f $OUT/parquet/$t/shift_speed=${sp}/reg=${reg}/n=${n}/part-${part_prefix}0.parquet
done
if (( archive_ok )); then
    rsync -a $arc $OUT/raw/ || copy_ok=0
fi
[[ -d $WD/parquet ]] && { rsync -a $WD/parquet $OUT/ || copy_ok=0; }
rsync -a $WD/summary/ $OUT/summary/ || copy_ok=0
if (( copy_ok && all_ok )); then
    date > $WD/_DONE
    rsync -a $WD/_DONE $OUT/summary/${combo}/ || copy_ok=0
fi

if (( ! copy_ok )); then
    echo "ERROR: copy to ${OUT} failed; rerun this task" >&2
    exit 1
fi

if (( ! all_ok )); then
    echo "ERROR: failed sims: ${failed_sims[*]:-none}; failed conversions: ${failed_types[*]:-none}; archive ok: ${archive_ok}" >&2
    exit 1
fi
echo "$(date +%T) all done: ${OUT}/summary/${combo}/_DONE"
