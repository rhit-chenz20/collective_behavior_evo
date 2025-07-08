#!/bin/bash

psi=$1
jobid=$2
echo "nonrec, SLURM_JOB_ID=$SLURM_JOB_ID, TASK_ID=$SLURM_ARRAY_TASK_ID"
# Example simulation command (e.g., Python, R, binary)
bin/slim5.0 -d psi11=$psi \
            -d ID="${SLURM_JOB_ID}" \
            -d "folder='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'"  \
            -d "burnin_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/pop/nrec_N_1000_psi_${psi}_pop_2702577.pop'" \
            -d "phenotype_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/phenotype/nrec_N_1000_psi_${psi}_pop_2702577_${SLURM_JOB_ID}.tsv'" \
            -d "fre_output='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/data/nrec_N_1000_psi_${psi}_pop_2702577_${SLURM_JOB_ID}.tsv'" \
            after_burnin_nonrec.slim &> /workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log/nrec_N_1000_psi_${psi}_pop_2702577_${SLURM_JOB_ID}.txt &

# bin/slim5.0 -d psi11=$psi \
#             -d ID="${SLURM_JOB_ID}" \
#             -d "folder='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}'"  \
#             -d "burnin_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/pop/nrec_N_1000_psi_${psi}_pop_2702581.pop'" \
#             -d "phenotype_fn='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/phenotype/nrec_N_1000_psi_${psi}_pop_2702581_${SLURM_JOB_ID}.tsv'" \
#             -d "fre_output='/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/data/nrec_N_1000_psi_${psi}_pop_2702581_${SLURM_JOB_ID}.tsv'" \
#             after_burnin_nonrec.slim &> /workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}/log/nrec_N_1000_psi_${psi}_pop_2702581_${SLURM_JOB_ID}.txt &

wait
