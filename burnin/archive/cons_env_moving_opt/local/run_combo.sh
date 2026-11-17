#!/bin/bash
# bash run_combo.sh <task index>, started by run_sim.sh from this folder; keep it working with
# macOS bash 3.2. One task = one (shift speed, n, psi, reg) combination (index order in
# pipeline_config.sh). Everything happens in $WD, which is removed at the end:
#   1. run the REPS replicates one after another      -> $WD/raw/<type>/*.tsv, $WD/raw/log/*.txt
#   2. convert each type to parquet                    -> $WD/parquet/<type>/shift_speed=../reg=../n=../part-psi_<psi>-0.parquet
#   3. archive the raw files and check the archive     -> $WD/<combo>.tar.zst
#   4. copy parquet, archive and summaries to $OUT     -> $OUT/{parquet,raw,summary}
# $OUT/summary/<combo>/_DONE is written last, only if every step succeeded.

set -u
TASK=$1
CODE_DIR=$PWD
source pipeline_config.sh

# ---- this task's parameters ----
n_sp=${#SPEEDS[@]}; n_n=${#NS[@]}; n_psi=${#PSIS[@]}; n_reg=${#REGS[@]}
r=$TASK
reg_i=$(( r % n_reg ));  r=$(( r / n_reg ))
psi_i=$(( r % n_psi ));  r=$(( r / n_psi ))
n_i=$(( r % n_n ));      r=$(( r / n_n ))
sp_i=$r
reg=${REGS[$reg_i]}; psi=${PSIS[$psi_i]}; n=${NS[$n_i]}; sp=${SPEEDS[$sp_i]}
combo="sp_${sp}_reg_${reg}_n_${n}_psi_${psi}"
echo "task ${TASK}: ${combo}, ${REPS} reps, on $(hostname)"

# scratch for this task; removed when it ends for any reason (TERM included)
WD=$WORK_ROOT/task_$TASK
mkdir -p $WD/raw $WD/pop
trap 'cd /; rm -rf $WD' EXIT
trap 'exit 143' TERM

# work from a copy, so editing this folder during a long run does not touch running tasks
cp -r ${CODE_DIR} $WD/code
for (( rep=1; rep<=REPS; rep++ )); do
    f=${POP_DIR}/cons_n_${n}_psi_${psi}_reg_${reg}_${rep}.pop
    [[ -f $f ]] && cp $f $WD/pop/
done
echo "copied inputs: $(ls $WD/pop | wc -l) pop files, $(du -sh $WD/code | cut -f1) code"

for d in "${ARCHIVE_TYPES[@]}"; do mkdir -p $WD/raw/$d; done
cd $WD/raw   # SLiM writes its outputs relative to here

# ---- 1. simulations ----
failed_sims=()
for (( rep=1; rep<=REPS; rep++ )); do
    # unique sim ID: index over rep > speed > n > psi > reg, plus 2
    idx=$(( (((((rep - 1) * n_sp + sp_i) * n_n + n_i) * n_psi + psi_i) * n_reg) + reg_i ))
    i=$(( idx + 2 ))
    tag="n_${n}_psi_${psi}_sp_${sp}_reg_${reg}_${rep}_${i}"
    POP_FILE=$WD/pop/cons_n_${n}_psi_${psi}_reg_${reg}_${rep}.pop
    echo "$(date +%T) rep=${rep} ID=${i}"
    if [[ ! -f "$POP_FILE" ]]; then
        echo "  missing pop file: ${POP_DIR}/$(basename $POP_FILE)"
        failed_sims+=($tag)
        continue
    fi

    $SLIM -d psi=$psi \
        -d ID="${i}" \
        -d rep="${rep}" \
        -d n=$n \
        -d shift_speed=$sp \
        -d "metadata_fn='metadata/${tag}.tsv'" \
        -d "phenotype_fn='phenotype/${tag}.tsv'" \
        -d "genotype_fn='genotype/${tag}.tsv'" \
        -d "fre_output='data/${tag}.tsv'" \
        -d "vars_fn='vars/${tag}.tsv'" \
        -d "burnin_fn='${POP_FILE}'" \
        -d "fitness_fn='fitness/${tag}.tsv'" \
        -d "group_tag_fn='group_tag/${tag}.tsv'" \
        -d "allele_fn='allele/${tag}.tsv'" \
        -d "copied_val_fn='copied_val/${tag}.tsv'" \
        $WD/code/${reg}.slim &> log/${tag}.txt
    status=$?
    echo "  exit status ${status}"
    if (( status != 0 )); then
        # drop any partial outputs so only complete sims reach the parquet; keep its log
        failed_sims+=($tag)
        rm -f */${tag}.tsv
    fi
done
cd $WD
echo "$(date +%T) sims done, ${#failed_sims[@]} failed"

# ---- 2. convert to parquet (one part file per type; one python process for all types) ----
$PY $WD/code/convert_parquet.py --layout moving --raw-dir $WD/raw --type "${TSV_TYPES[@]}" --out-dir $WD/parquet \
    --workers 1 --sims-per-part $REPS --part-prefix "psi_${psi}-" \
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
tar -cf - -C $WD/raw . | $ZSTD -T1 -10 -q -o $arc
n_archived=$($ZSTD -dc $arc | tar -tf - | grep -vc '/$')
if (( n_archived != n_files )); then
    echo "ERROR: archive has ${n_archived} files, raw has ${n_files}" >&2
    archive_ok=0
else
    archive_ok=1
    echo "$(date +%T) archive ok: ${n_files} files, $(du -h $arc | cut -f1)"
fi

# ---- 4. copy to $OUT ----
# rsync writes each file under a temporary name and renames it when complete, so a
# killed job never leaves a truncated file under its final name
mkdir -p $WD/summary/${combo}
# (the :+ form keeps an empty array from being an unbound variable in bash 3.2)
printf '%s\n' ${failed_sims[@]+"${failed_sims[@]}"} > $WD/summary/${combo}/failed_sims.txt
all_ok=1
(( ${#failed_sims[@]} || ${#failed_types[@]} || ! archive_ok )) && all_ok=0

copy_ok=1
mkdir -p $OUT/raw $OUT/summary/${combo} || copy_ok=0
# a rerun first removes this combo's outputs from an earlier run (exact paths), so $OUT never mixes old files with the new run's results
rm -f $OUT/summary/${combo}/* $OUT/raw/${combo}.tar.zst
for t in "${TSV_TYPES[@]}"; do
    rm -f $OUT/parquet/$t/shift_speed=${sp}/reg=${reg}/n=${n}/part-psi_${psi}-0.parquet
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
