#!/bin/bash

echo "SLURM_JOB_ID=$SLURM_JOB_ID, TASK_ID=$SLURM_ARRAY_TASK_ID"

bin/slim5.0 -d ID=$SLURM_JOB_ID -d "folder='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'" newmodel.slim 

# Copy results back to permanent storage if needed
