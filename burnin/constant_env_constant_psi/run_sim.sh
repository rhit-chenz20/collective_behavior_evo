#!/bin/bash

dir_name="constant_env_ext"

mkdir -p ../../data/${dir_name}
mkdir -p ../../data/${dir_name}/phenotype
mkdir -p ../../data/${dir_name}/genotype
mkdir -p ../../data/${dir_name}/data
mkdir -p ../../data/${dir_name}/ind
mkdir -p ../../data/${dir_name}/mut
mkdir -p ../../data/${dir_name}/log
mkdir -p ../../data/${dir_name}/pop
mkdir -p ../../data/${dir_name}/vars
mkdir -p ../../data/${dir_name}/fitness
mkdir -p ../../data/${dir_name}/allele

python make_model.py 

reset=false # reset the simulation even if the .pop file exists

c=0
for rep in {1..5}; do  
   for sz in 20; do   
        for n in 2 5 8 20; do
            for psi in 0.1 0.15 0.2 0.25 0.3; do
                for reg in ext; do
                    echo "Submitting job, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}, reg=${reg}."
                
                    POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.pop"

                    if [[  "$reset" == false && -f "$POP_FILE" ]]; then
                        # ---- read the pop file if the pop file EXISTS and reset is false ----
                        echo "Found pop file: $POP_FILE"

                        bin/slim5.0 -d psi=$psi \
                        -d ID="${rep}" \
                        -d n=$n \
                        -d shift_size=$sz \
                        -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                        -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                        -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "burnin_fn='$POP_FILE'" \
                        -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "allele_fn='../../data/${dir_name}/allele/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                        ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.txt &

                    else
                        # ---- run simulation from the start if the file DOES NOT exist or reset is true ----
                        echo "Missing pop file: $POP_FILE"

                        bin/slim5.0 -d psi=$psi \
                        -d ID="${rep}" \
                        -d n=$n \
                        -d shift_size=$sz \
                        -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                        -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                        -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                    -d "pop_output='$POP_FILE'" \
                        -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "allele_fn='../../data/${dir_name}/allele/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                        ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.txt &

                    fi
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
