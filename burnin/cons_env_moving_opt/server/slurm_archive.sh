#!/bin/bash -l
#SBATCH --account=bscb10
#SBATCH --partition=regular,long7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=cons_archive
#SBATCH --output=/home/%u/slurm_logs/cons_archive_%A_%a.out
#SBATCH --error=/home/%u/slurm_logs/cons_archive_%A_%a.err

# One array task per folder (index into ARCHIVE_TYPES): pack $OUT/<type>/ into
# $ARCHIVE_DIR/<type>.tar.zst, verify it, and (if DELETE_RAW=1) remove the loose files.
# ARCHIVE_DIR defaults to $OUT/raw; set it to archive loose files left elsewhere (e.g. on /home)
# straight to lab storage:  sbatch --export=ALL,OUT=<loose dir>,ARCHIVE_DIR=<lab>/raw ...

set -euo pipefail
source ${SLURM_SUBMIT_DIR}/pipeline_config.sh
mount_lab_storage

type=${ARCHIVE_TYPES[$SLURM_ARRAY_TASK_ID]}
src=${OUT}/${type}
ARCHIVE_DIR=${ARCHIVE_DIR:-${OUT}/raw}
arc=${ARCHIVE_DIR}/${type}.tar.zst
mkdir -p ${ARCHIVE_DIR}

if [[ ! -d $src ]]; then
    if [[ -f $arc ]]; then
        echo "${src} not found; already archived in ${arc}"
        exit 0
    fi
    echo "ERROR: neither ${src} nor ${arc} exists (wrong OUT?)" >&2
    exit 1
fi
# never overwrite a finished archive: if a run was stopped during rm -rf, the folder is
# incomplete and re-archiving it would replace the good archive with a partial one
if [[ -f $arc ]]; then
    echo "ERROR: ${arc} already exists; check it and remove it by hand to re-archive ${src}" >&2
    exit 1
fi
# a .part file is either left by a cancelled run (zstd will not overwrite it, so remove it)
# or still being written by another running job (leave it alone)
if [[ -e ${arc}.part ]]; then
    if [[ -n $(find ${arc}.part -mmin -10) ]]; then
        echo "ERROR: ${arc}.part was written in the last 10 min; another job is probably archiving ${type}" >&2
        exit 1
    fi
    rm -f ${arc}.part
fi

n_files=$(find $src -type f | wc -l)
echo "archiving ${n_files} files from ${src}"
tar -cf - -C ${OUT} ${type} | $ZSTD -T${SLURM_CPUS_PER_TASK} -10 -q -o ${arc}.part
mv ${arc}.part ${arc}

# check the archive decompresses and holds every file
n_archived=$($ZSTD -dc ${arc} | tar -tf - | grep -vc '/$')
if (( n_archived != n_files )); then
    echo "ERROR: archive has ${n_archived} files, folder has ${n_files}; keeping ${src}" >&2
    exit 1
fi
echo "archive ok: ${n_archived} files, $(du -h ${arc} | cut -f1)"

if (( DELETE_RAW != 1 )); then
    echo "DELETE_RAW=${DELETE_RAW}; keeping ${src}"
    exit 0
fi
# tsv folders are only deleted once their parquet conversion has succeeded
if [[ " ${TSV_TYPES[*]} " == *" ${type} "* && ! -f ${OUT}/parquet/${type}/_SUCCESS.json ]]; then
    echo "ERROR: no ${OUT}/parquet/${type}/_SUCCESS.json; keeping ${src}" >&2
    exit 1
fi
rm -rf ${src}
echo "removed ${src}"
