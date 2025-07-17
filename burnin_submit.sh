#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=bscb10
#SBATCH --job-name=n_group_burnin
#SBATCH --array=1-3
#SBATCH --output="/home/zc524/slurm-outputs/%x-%j-%a.out"
#SBATCH --time=24:00:00      
#SBATCH --cpus-per-task=22  
#SBATCH --mem=22G              
#SBATCH --partition=short,long7
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

for i in {1..2}; do
    echo "Running burnin simulations for simid=${i}."

    for psi in 0.1 0.3 0.5 0.7 0.9 0.0 -0.1 -0.3 -0.5 -0.7 -0.9; do
        echo "Running simulation with psi=${psi}"
        for n in 3 11;do
            bash n_group_bi.sh ${psi} ${n} &
        done
        # # burnin phase
        # bash burnin_rec.sh ${psi} &
        # bash burnin_relax_sel.sh ${psi} &
        # bash non0_burnin.sh ${psi} &
        # bash sim_nonrec_slurm.sh ${psi} &
    done
    wait
done


echo "Moving results into storage."
mkdir -p ${RESULTSHOME}/pop
mkdir -p ${RESULTSHOME}/log

zip -r pop_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip pop
zip -r log_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip log

mv log_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip $RESULTSHOME/log/
mv pop_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.zip $RESULTSHOME/pop/

# mv log/*.txt $RESULTSHOME/log/
# mv pop/*.pop $RESULTSHOME/pop/

echo "Cleaning up working directory..."
rm -r $WORKDIR
