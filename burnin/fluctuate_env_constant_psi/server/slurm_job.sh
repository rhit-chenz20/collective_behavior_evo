#!/bin/bash -l
#SBATCH --account=bscb10
#SBATCH --partition=regular,long7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --exclude=cbsubscb[02-04,06-12]
#SBATCH --job-name=fluc_collect
#SBATCH --output=/home/%u/slurm_logs/fluc_collect_%A_%a.out
#SBATCH --error=/home/%u/slurm_logs/fluc_collect_%A_%a.err

# One array task = one (gap, sz, n, psi, reg) combination x one chunk of REPS_PER_TASK
# replicates (order in pipeline_config.sh). Sims start from scratch (no burnin file); 10 reps
# took 4.5-10.4 h on average per node and up to ~12 h (Oct 2026). The --time and --exclude
# above are defaults: submit_slurm.sh and rerun_failed.sh pass TIME_LIMIT and EXCLUDE_NODES
# from pipeline_config.sh.
# Everything happens in local /workdir. Inputs (this folder incl. SLiM) are copied from /home
# once at the start; lab storage is only written to, never read:
#   1. run the replicates one after another   -> $WD/raw/<type>/*.tsv, $WD/raw/log/*.txt
#   2. convert each type to parquet            -> $WD/parquet/<type>/fluc_opt=../fluc_interval=../reg=../n=../part-psi_<psi>-reps_<first>-0.parquet
#   3. archive the SLiM logs and check it      -> $WD/<chunk>_logs.tar.zst  (raw tsv are not kept)
#   4. copy everything to $OUT                 -> $OUT/{parquet,logs,summary/<chunk>}
# $OUT/summary/<chunk>/_DONE is written last, only if every step succeeded.
# Submit via submit_slurm.sh from this folder; settings come from pipeline_config.sh.

set -u
source ${SLURM_SUBMIT_DIR}/pipeline_config.sh

# ---- this task's parameters (submit_slurm.sh splits the tasks over several arrays) ----
n_gap=${#GAPS[@]}; n_sz=${#SZS[@]}; n_n=${#NS[@]}; n_psi=${#PSIS[@]}; n_reg=${#REGS[@]}
task=$(( ${TASK_OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
chunk=$(( task % N_CHUNKS ))
r=$(( task / N_CHUNKS ))
reg_i=$(( r % n_reg ));  r=$(( r / n_reg ))
psi_i=$(( r % n_psi ));  r=$(( r / n_psi ))
n_i=$(( r % n_n ));      r=$(( r / n_n ))
sz_i=$(( r % n_sz ));    r=$(( r / n_sz ))
gap_i=$r
reg=${REGS[$reg_i]}; psi=${PSIS[$psi_i]}; n=${NS[$n_i]}; sz=${SZS[$sz_i]}; gap=${GAPS[$gap_i]}
rep_first=$(( chunk * REPS_PER_TASK + 1 ))
rep_last=$(( rep_first + REPS_PER_TASK - 1 ))
(( rep_last > REPS )) && rep_last=$REPS
combo="gap_${gap}_sz_${sz}_reg_${reg}_n_${n}_psi_${psi}"
part="${combo}_reps_${rep_first}-${rep_last}"
echo "task ${task}: ${combo}, reps ${rep_first}-${rep_last}, on $(hostname)"

# stop before the sims if this node cannot run the parquet conversion: on nodes whose CPUs
# lack AVX2 etc., importing polars dies with "Illegal instruction" (see EXCLUDE_NODES)
if ! $PY -c "import polars" > /dev/null 2>&1; then
    echo "ERROR: polars cannot run on $(hostname); add it to EXCLUDE_NODES in pipeline_config.sh and rerun this task" >&2
    exit 1
fi

# local scratch for this job; removed when the job ends for any reason (done, error,
# scancel or time limit: SLURM sends TERM first, which the trap turns into a normal exit)
WD=/workdir/$USER/$SLURM_JOB_ID
mkdir -p $WD/raw
trap 'cd /; rm -rf $WD' EXIT
trap 'exit 143' TERM
# read inputs from /home once, up front: this folder (scripts, and bin/slim5.0, a symlink
# that -L copies as the real binary)
cp -rL ${SLURM_SUBMIT_DIR} $WD/code
SLIM=$WD/code/bin/slim5.0
echo "copied inputs: $(du -sh $WD/code | cut -f1) code"

for d in "${ARCHIVE_TYPES[@]}"; do mkdir -p $WD/raw/$d; done
cd $WD/raw   # SLiM writes its outputs relative to here

# ---- 1. simulations ----
failed_sims=()
for (( rep=rep_first; rep<=rep_last; rep++ )); do
    # unique sim ID: index over rep > gap > sz > n > psi > reg, plus 1
    idx=$(( ((((((rep - 1) * n_gap + gap_i) * n_sz + sz_i) * n_n + n_i) * n_psi + psi_i) * n_reg) + reg_i ))
    i=$(( idx + 1 ))
    tag="n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}_${i}"
    echo "$(date +%T) rep=${rep} ID=${i}"

    $SLIM -d psi=$psi \
        -d ID="${i}" \
        -d rep="${rep}" \
        -d n=$n \
        -d fluc_opt=$sz \
        -d fluc_interval=$gap \
        -d "metadata_fn='metadata/${tag}.tsv'" \
        -d "phenotype_fn='phenotype/${tag}.tsv'" \
        -d "genotype_fn='genotype/${tag}.tsv'" \
        -d "fre_output='data/${tag}.tsv'" \
        -d "allele_fn='allele/${tag}'" \
        -d "group_tag_fn='group_tag/${tag}.tsv'" \
        $WD/code/${reg}.slim &> log/${tag}.txt
    status=$?
    # SLiM also exits 0 when it stops early (e.g. "max L has reached"), so check the outputs
    missing=()
    for t in "${TSV_TYPES[@]}"; do [[ -s $t/${tag}.tsv ]] || missing+=($t); done
    echo "  exit status ${status}${missing:+, missing: ${missing[*]}}"
    if (( status != 0 || ${#missing[@]} )); then
        # drop any partial outputs so only complete sims reach the parquet; keep its log
        failed_sims+=($tag)
        rm -f */${tag}.tsv
    fi
done
cd $WD
echo "$(date +%T) sims done, ${#failed_sims[@]} failed"

# ---- 2. convert to parquet (one part file per type; one python process for all types) ----
$PY $WD/code/convert_parquet.py --layout fluc --raw-dir $WD/raw --type "${TSV_TYPES[@]}" \
    --out-dir $WD/parquet --workers 1 --sims-per-part $REPS_PER_TASK \
    --part-prefix "psi_${psi}-reps_${rep_first}-" \
    --summary-dir $WD/summary/${part} > $WD/convert.log 2>&1
# a type succeeded iff its summary was written
failed_types=()
for t in "${TSV_TYPES[@]}"; do
    [[ -f $WD/summary/${part}/${t}.json ]] || failed_types+=($t)
done
if (( ${#failed_types[@]} )); then
    echo "  convert FAILED for ${failed_types[*]}:"; grep ERROR $WD/convert.log
fi
echo "$(date +%T) convert done, ${#failed_types[@]} types failed"

# ---- 3. archive the SLiM logs and check the archive holds every log ----
# the raw tsv are not kept: the parquet holds the data, and they go with $WD at the end
arc=$WD/${part}_logs.tar.zst
n_files=$(find $WD/raw/log -type f | wc -l)
tar -cf - -C $WD/raw log | $ZSTD -T${SLURM_CPUS_PER_TASK:-1} -10 -q -o $arc
n_archived=$($ZSTD -dc $arc | tar -tf - | grep -vc '/$')
if (( n_archived != n_files )); then
    echo "ERROR: archive has ${n_archived} files, log/ has ${n_files}" >&2
    archive_ok=0
else
    archive_ok=1
    echo "$(date +%T) log archive ok: ${n_files} files, $(du -h $arc | cut -f1)"
fi

# ---- 4. copy to lab storage (write only) ----
# rsync writes each file under a temporary name and renames it when complete, so a
# killed job never leaves a truncated file under its final name
mkdir -p $WD/summary/${part}
printf '%s\n' "${failed_sims[@]}" > $WD/summary/${part}/failed_sims.txt
all_ok=1
(( ${#failed_sims[@]} || ${#failed_types[@]} || ! archive_ok )) && all_ok=0

mount_lab_storage
copy_ok=1
mkdir -p $OUT/logs $OUT/summary/${part} || copy_ok=0
# a rerun first removes this chunk's outputs from an earlier run (exact paths, nothing is
# listed or read), so storage never mixes old files with the new run's results
rm -f $OUT/summary/${part}/* $OUT/logs/${part}_logs.tar.zst
for t in "${TSV_TYPES[@]}"; do
    rm -f $OUT/parquet/$t/fluc_opt=${sz}/fluc_interval=${gap}/reg=${reg}/n=${n}/part-psi_${psi}-reps_${rep_first}-0.parquet
done
if (( archive_ok )); then
    rsync -a $arc $OUT/logs/ || copy_ok=0
fi
[[ -d $WD/parquet ]] && { rsync -a $WD/parquet $OUT/ || copy_ok=0; }
rsync -a $WD/summary/ $OUT/summary/ || copy_ok=0
if (( copy_ok && all_ok )); then
    date > $WD/_DONE
    rsync -a $WD/_DONE $OUT/summary/${part}/ || copy_ok=0
fi

if (( ! copy_ok )); then
    echo "ERROR: copy to ${OUT} failed; rerun this task" >&2
    exit 1
fi

if (( ! all_ok )); then
    echo "ERROR: failed sims: ${failed_sims[*]:-none}; failed conversions: ${failed_types[*]:-none}; archive ok: ${archive_ok}" >&2
    exit 1
fi
echo "$(date +%T) all done: ${OUT}/summary/${part}/_DONE"
