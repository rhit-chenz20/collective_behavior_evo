#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=bscb10
#SBATCH --job-name=postburnin
#SBATCH --array=1-3
#SBATCH --output="/home/zc524/slurm-outputs/%x-%j-%a.out"
#SBATCH --time=4:00:00      
#SBATCH --cpus-per-task=12    
#SBATCH --mem=2000              
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
cp -r ~/collective_behavior_evo/slim/burnin/* .
echo "Linking SLiM executable."
mkdir bin
ln -s ~/bin/slim5.0 bin/slim5.0

echo "Running simulations."

for psi in 0.0 0.1 0.3 0.5 0.7 0.9 -0.1 -0.3 -0.5 -0.7 -0.9; do
    echo "Running simulation with psi=${psi}"

    # # burnin phase
    # bash sim_slurm.sh ${psi} &
    # bash sim_nonrec_slurm.sh ${psi} &

    for jobid in {2702577..2702583}; do
        for i in {1..2}; do
            echo "Submitting job for psi=${psi}, jobid=${jobid}, simid=${i}."
            # burnin phase 
            bash burnin_rec.sh ${psi} ${jobid} ${i} &
            # bash burnin_nrec.sh ${psi} ${jobid} ${i} &
        done
    done
    # bash post_burnin_nrec.sh ${psi} &
    echo "Submitted jobs for psi=${psi}."
    wait
    echo "All jobs for psi=${psi} completed."
done

echo "Moving results into storage."
mkdir -p ${RESULTSHOME}/data
mkdir -p ${RESULTSHOME}/pop
mkdir -p ${RESULTSHOME}/phenotype
mkdir -p ${RESULTSHOME}/genotype
mkdir -p ${RESULTSHOME}/log
mkdir -p ${RESULTSHOME}/mut_ind
mv data/*.tsv $RESULTSHOME/data/
mv log/*.txt $RESULTSHOME/log/
mv pop/*.pop $RESULTSHOME/pop/
mv phenotype/*.tsv $RESULTSHOME/phenotype/
mv genotype/*.tsv $RESULTSHOME/genotype/
mv mut_ind/* $RESULTSHOME/mut_ind/

echo "Cleaning up working directory..."
rm -r $WORKDIR
