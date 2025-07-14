#!/bin/bash

psi=$1
echo "rec, SLURM_JOB_ID=$SLURM_JOB_ID, TASK_ID=$SLURM_ARRAY_TASK_ID"
# Example simulation command (e.g., Python, R, binary)
bin/slim5.0 -d psi11=$psi -d ID="${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}" -d "folder='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'" social_new.slim &> /workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/rec_log_${SLURM_JOB_ID}.txt 

# Copy results back to permanent storage if needed
