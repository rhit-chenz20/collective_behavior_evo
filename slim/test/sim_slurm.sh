#!/bin/bash
echo "Test mode: SLURM_JOB_ID=$SLURM_JOB_ID, TASK_ID=$SLURM_ARRAY_TASK_ID"
bin/slim5.0 -d ID=${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}-$1 -d "folder='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'" test.slim &> /workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log_${SLURM_JOB_ID}.txt 
exit 0