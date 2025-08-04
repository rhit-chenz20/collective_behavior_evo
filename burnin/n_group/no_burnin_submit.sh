#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=bscb10
#SBATCH --job-name=noburnin_ngroup
#SBATCH --array=1-20
#SBATCH --output="/home/zc524/slurm-outputs/%x-%j-%a.out"
#SBATCH --time=24:00:00      
#SBATCH --cpus-per-task=7   
#SBATCH --mem=14G              
#SBATCH --partition=short
#SBATCH --mail-user=zc524@cornell.edu
#SBATCH --mail-type=FAIL,END

echo "Workstation is ${HOSTNAME}, partition is ${SLURM_JOB_PARTITION}."

# Create working directory and the destination folder for results.
WORKDIR=/workdir/$USER/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}
DATAHOME=/fs/cbsubscb10/storage/zc524/col_beh
RESULTSHOME=${DATAHOME}/social/${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}

# Create relevant directory structure
mkdir -p ${WORKDIR}
cd ${WORKDIR}
echo "Working directory is ${WORKDIR}."

# Mount storage
/programs/bin/labutils/mount_server cbsubscb10 /storage

echo "Copying analysis scripts."
cp -r ~/collective_behavior_evo/slim/burnin/n_group/* .
unzip pop.zip
echo "Linking SLiM executable."
mkdir bin
ln -s ~/bin/slim5.0 bin/slim5.0

echo "Running simulations."

for n in 2 3 5 7 9; do
    mkdir -p phenotype/n_${n}
    mkdir -p genotype/n_${n}
    mkdir -p data/n_${n}
    mkdir -p ind/n_${n}
    mkdir -p mut/n_${n}
    mkdir -p log/n_${n}
done


for n in 2 3 5 7 9; do
    echo "Running simulation with n=${n}"
    for psi in 0.0 0.3 0.6 0.9 -0.3 -0.6 -0.9; do
        echo "Submitting job for psi=${psi}, n=${n}."
        bash no_burnin.sh ${psi} ${n} &
    done
    wait
done


echo "Moving results into storage."
mkdir -p ${RESULTSHOME}/data
# mkdir -p ${RESULTSHOME}/pop
mkdir -p ${RESULTSHOME}/phenotype
mkdir -p ${RESULTSHOME}/genotype
mkdir -p ${RESULTSHOME}/log
mkdir -p ${RESULTSHOME}/ind
mkdir -p ${RESULTSHOME}/mut

zip -r data_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip data
# zip -r pop_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip pop
zip -r phenotype_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip phenotype
zip -r genotype_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip genotype
zip -r log_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip log
zip -r ind_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip ind
zip -r mut_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip mut

mv data_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip $RESULTSHOME/data/
# mv pop_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip $RESULTSHOME/pop/
mv phenotype_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip $RESULTSHOME/phenotype/
mv genotype_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip $RESULTSHOME/genotype/
mv log_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip $RESULTSHOME/log/
mv ind_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip $RESULTSHOME/ind/
mv mut_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip $RESULTSHOME/mut/

# mv data/* $RESULTSHOME/data/
# # mv pop/* $RESULTSHOME/pop/
# mv phenotype/* $RESULTSHOME/phenotype/
# mv genotype/* $RESULTSHOME/genotype/
# mv log/* $RESULTSHOME/log/
# mv ind/* $RESULTSHOME/ind/

echo "Cleaning up working directory..."
rm -r $WORKDIR
