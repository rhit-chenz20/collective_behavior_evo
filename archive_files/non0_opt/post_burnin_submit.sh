#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=bscb10
#SBATCH --job-name=non0_opt_post
#SBATCH --array=1-3
#SBATCH --output="/home/zc524/slurm-outputs/%x-%j-%a.out"
#SBATCH --time=4:00:00      
#SBATCH --cpus-per-task=3  
#SBATCH --mem=24G              
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
cp -r ~/collective_behavior_evo/slim/burnin/non0_opt/* .
cp -r ~/collective_behavior_evo/slim/burnin/split_ind_data.py .

mkdir -p log
mkdir -p data
mkdir -p phenotype
mkdir -p genotype
mkdir -p ind
mkdir -p mut

unzip pop.zip
echo "Linking SLiM executable."
mkdir bin
ln -s ~/bin/slim5.0 bin/slim5.0

echo "Running simulations."

for psi in 0.5 0.7 0.9; do
# for psi in 0.0 0.1 -0.1; do
    echo "Running simulation with psi=${psi}"

    for jobid in 2729516; do
        for i in 1; do
            echo "Submitting job for psi=${psi}, jobid=${jobid}, simid=${i}."
            # burnin phase 
            bash post_burnin.sh ${psi} ${jobid} ${i} &
        done
    done
done
wait

for psi in 0.0 0.1 0.3; do
# for psi in 0.0 0.1 -0.1; do
    echo "Running simulation with psi=${psi}"

    for jobid in 2729516; do
        for i in 1; do
            echo "Submitting job for psi=${psi}, jobid=${jobid}, simid=${i}."
            # burnin phase 
            bash post_burnin.sh ${psi} ${jobid} ${i} &
        done
    done
done
wait

for psi in -0.7 -0.9; do
# for psi in 0.0 0.1 -0.1; do
    echo "Running simulation with psi=${psi}"

    for jobid in 2729516; do
        for i in 1; do
            echo "Submitting job for psi=${psi}, jobid=${jobid}, simid=${i}."
            # burnin phase 
            bash post_burnin.sh ${psi} ${jobid} ${i} &
        done
    done
done
wait

for psi in -0.1 -0.3 -0.5 ; do
# for psi in 0.0 0.1 -0.1; do
    echo "Running simulation with psi=${psi}"

    for jobid in 2729516; do
        for i in 1; do
            echo "Submitting job for psi=${psi}, jobid=${jobid}, simid=${i}."
            # burnin phase 
            bash post_burnin.sh ${psi} ${jobid} ${i} &
        done
    done
done
wait

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
