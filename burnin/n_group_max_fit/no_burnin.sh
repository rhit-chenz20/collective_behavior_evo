#!/bin/bash

psi=$1
n=$2
echo "no burnin max fit, SLURM_JOB_ID=$SLURM_JOB_ID, TASK_ID=$SLURM_ARRAY_TASK_ID"

bin/slim5.0 -d psi=$psi \
            -d ID="${SLURM_JOB_ID}" \
            -d n=$n \
            -d "phenotype_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/phenotype/n_${n}/psi_${psi}_${SLURM_JOB_ID}.tsv'" \
            -d "genotype_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/genotype/n_${n}/psi_${psi}_${SLURM_JOB_ID}.tsv'" \
            -d "fre_output='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/data/n_${n}/psi_${psi}_${SLURM_JOB_ID}.tsv'" \
            -d "ind_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/ind/n_${n}/psi_${psi}_${SLURM_JOB_ID}.tsv'" \
            -d "mut_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/mut/n_${n}/psi_${psi}_${SLURM_JOB_ID}.tsv'" \
            post_burnin_n_group.slim &> /workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log/n_${n}/psi_${psi}_${SLURM_JOB_ID}.txt &

wait