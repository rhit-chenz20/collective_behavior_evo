#!/bin/bash

alpha=$1
n=$2

mkdir -p "/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/pop"
mkdir -p "/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log"

echo "rec, SLURM_JOB_ID=$SLURM_JOB_ID, TASK_ID=$SLURM_ARRAY_TASK_ID"
# Example simulation command (e.g., Python, R, binary)
bin/slim5.0 -d alpha=$alpha -d n=$n -d ID="${SLURM_JOB_ID}" -d "folder='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/pop'" n_group_bi.slim &> /workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log/rec_alpha_${alpha}_n_${n}log_${SLURM_JOB_ID}.txt 

# Copy results back to permanent storage if needed
