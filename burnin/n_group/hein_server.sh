#!/bin/bash

ln -s ~/bin/slim5.0 bin/slim5.0

mkdir -p ../../data/n_group_base_data/phenotype
mkdir -p ../../data/n_group_base_data/genotype
mkdir -p ../../data/n_group_base_data/data
mkdir -p ../../data/n_group_base_data/ind
mkdir -p ../../data/n_group_base_data/mut
mkdir -p ../../data/n_group_base_data/log

for rep in {1..20}; do
    for n in 9 12; do
        echo "Running simulation with n=${n}"
        for psi in 0.0 0.3 0.6 0.9 -0.3 -0.6 -0.9; do
            echo "Submitting job for psi=${psi}, n=${n}."
            # bash no_burnin.sh ${psi} ${n} &

            bin/slim5.0 -d psi=$psi \
                -d ID="${rep}" \
                -d n=$n \
                -d "phenotype_fn='../../data/n_group_base_data/phenotype/n_${n}_psi_${psi}_${rep}.tsv'" \
                -d "genotype_fn='../../data/n_group_base_data/genotype/n_${n}_psi_${psi}_${rep}.tsv'" \
                -d "fre_output='../../data/n_group_base_data/data/n_${n}_psi_${psi}_${rep}.tsv'" \
                -d "ind_fn='../../data/n_group_base_data/ind/n_${n}_psi_${psi}_${rep}.tsv'" \
                -d "mut_fn='../../data/n_group_base_data/mut/n_${n}_psi_${psi}_${rep}.tsv'" \
                post_burnin_n_group.slim &> ../../data/n_group_base_data/log/n_${n}_psi_${psi}_${rep}.txt &
        done
    done
    wait
done
