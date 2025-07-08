#!/bin/bash

psi=$1
burnin_id=$2
simid=$3
echo "rec, SLURM_JOB_ID=$SLURM_JOB_ID, TASK_ID=$SLURM_ARRAY_TASK_ID"

bin/slim5.0 -d psi11=$psi \
            -d ID="${SLURM_JOB_ID}" \
            -d "burnin_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/pop/rec_N_1000_psi_${psi}_pop_${burnin_id}.pop'" \
            -d "phenotype_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/phenotype/rec_N_1000_psi_${psi}_pop_${burnin_id}_${simid}.tsv'" \
            -d "genotype_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/genotype/rec_N_1000_psi_${psi}_pop_${burnin_id}_${simid}.tsv'" \
            -d "fre_output='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/data/rec_N_1000_psi_${psi}_pop_${burnin_id}_${simid}.tsv'" \
            -d "mutation_prefix='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/mut_ind/rec_N_1000_psi_${psi}_pop_${burnin_id}_${simid}'" \
            after_burnin_rec.slim &> /workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log/rec_N_1000_psi_${psi}_pop_${burnin_id}_${simid}.txt &

# 2702577
# bin/slim5.0 -d psi11=$psi \
            # -d ID="${SLURM_JOB_ID}" \
#             -d "folder='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'"  \
#             -d "burnin_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/pop/rec_N_1000_psi_${psi}_pop_2702581.pop'" \
#             -d "phenotype_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/phenotype/rec_N_1000_psi_${psi}_pop_2702581_${SLURM_JOB_ID}.tsv'" \
#             -d "fre_output='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/data/rec_N_1000_psi_${psi}_pop_2702581_${SLURM_JOB_ID}.tsv'" \
#             after_burnin_rec.slim &> /workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log/rec_N_1000_psi_${psi}_pop_2702581_${SLURM_JOB_ID}.txt &

wait

