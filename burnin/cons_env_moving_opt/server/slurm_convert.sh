#!/bin/bash -l
#SBATCH --account=bscb10
#SBATCH --partition=regular,long7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --job-name=cons_convert
#SBATCH --output=/home/%u/slurm_logs/cons_convert_%A_%a.out
#SBATCH --error=/home/%u/slurm_logs/cons_convert_%A_%a.err

# One array task per output type (index into TSV_TYPES): tsv -> partitioned parquet.
# Builds the dataset on local /workdir, then copies it to $OUT/parquet/<type>.

set -euo pipefail
source ${SLURM_SUBMIT_DIR}/pipeline_config.sh
mount_lab_storage

type=${TSV_TYPES[$SLURM_ARRAY_TASK_ID]}
WD=/workdir/$USER/$SLURM_JOB_ID
mkdir -p $WD ${OUT}/parquet
trap 'rm -rf $WD' EXIT

# one worker per core; each holds ~100 sims of one type in memory (well under 2 GB)
$PY ${SLURM_SUBMIT_DIR}/convert_parquet.py \
    --raw-dir ${OUT} \
    --type ${type} \
    --out-dir ${WD}/parquet \
    --workers ${SLURM_CPUS_PER_TASK}

rm -rf ${OUT}/parquet/${type}
cp -r ${WD}/parquet/${type} ${OUT}/parquet/
echo "copied to ${OUT}/parquet/${type}"
