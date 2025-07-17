#!/bin/bash
mkdir -p "/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/pop"
mkdir -p "/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log"

psi=$1
echo "rec, SLURM_JOB_ID=$SLURM_JOB_ID, TASK_ID=$SLURM_ARRAY_TASK_ID"
# Example simulation command (e.g., Python, R, binary)
bin/slim5.0 -d psi11=$psi -d ID="${SLURM_JOB_ID}" -d "folder='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/pop'" non0_burnin.slim &> /workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log/rec_psi_${psi}_log_${SLURM_JOB_ID}.txt 

# Copy results back to permanent storage if needed
