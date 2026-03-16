#!/bin/bash

dir_name="e_psi_with_bg_psi"

mkdir -p ../../data/${dir_name}
mkdir -p ../../data/${dir_name}/phenotype
mkdir -p ../../data/${dir_name}/genotype
mkdir -p ../../data/${dir_name}/data
mkdir -p ../../data/${dir_name}/ind
mkdir -p ../../data/${dir_name}/mut
mkdir -p ../../data/${dir_name}/log
mkdir -p ../../data/${dir_name}/pop
mkdir -p ../../data/${dir_name}/vars
mkdir -p ../../data/${dir_name}/psis
mkdir -p ../../data/${dir_name}/group_psis
mkdir -p ../../data/${dir_name}/fitness

python make_model.py 

c=0
for rep in {1..5}; do  
   for sz in 40 80 120 300; do   
        for n in 2 8 20 50; do
            for psi in 0.0 0.2 0.5 0.8; do
                for reg in ave fit inv ext; do
                    if reg == "inv" && (n == 50 || n == 20); then
                        continue
                    fi
                    echo "Submitting job, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}, reg=${reg}."
                
                    bin/slim5.0 -d start_psi=$psi \
                    -d ID="${rep}" \
                    -d n=$n \
                    -d shift_size=$sz \
                    -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                    -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                    -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                    -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                    -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                    -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                    -d "pop_output='../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.pop'" \
                    -d "psi_fn='../../data/${dir_name}/psis/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                    -d "gpsi_fn='../../data/${dir_name}/group_psis/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                    -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                    ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.txt &
                    ((c++))
                    if (( c > 9)); 
                    then
                        wait
                        c=0
                    fi
                done
            done
        done
    done
done
