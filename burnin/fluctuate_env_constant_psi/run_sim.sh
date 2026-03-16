#!/bin/bash

# dir_name="fluctuate_env_constant_psi"

# mkdir -p ../../data/${dir_name}
# mkdir -p ../../data/${dir_name}/phenotype
# mkdir -p ../../data/${dir_name}/genotype
# mkdir -p ../../data/${dir_name}/data
# mkdir -p ../../data/${dir_name}/ind
# mkdir -p ../../data/${dir_name}/mut
# mkdir -p ../../data/${dir_name}/log
# mkdir -p ../../data/${dir_name}/pop
# mkdir -p ../../data/${dir_name}/vars
# mkdir -p ../../data/${dir_name}/fitness

reset=true # reset the simulation even if the .pop file exists

python make_model.py 

c=0
for rep in {6..50}; do  
    for gap in 50 150 100; do
        dir_name="fluctuate_env_gap_${gap}"
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

        for sz in 20 40 60; do   
                for n in 2 8 20 50; do
                    for psi in 0.0 0.5 0.7 0.9 ; do
                        for reg in ave fit ext; do
                            echo "Submitting job, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}, reg=${reg}, gap=${gap}."

                            POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.pop"

                            if [[  "$reset" == false && -f "$POP_FILE" ]]; then
                                # ---- read the pop file if the pop file EXISTS and reset is false ----
                                echo "Found pop file: $POP_FILE"

                                bin/slim5.0 -d psi=$psi \
                                -d ID="${rep}" \
                                -d n=$n \
                                -d fluc_opt=$sz \
                                -d fluc_interval=$gap \
                                -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                                -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                                -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                -d "burnin_fn='$POP_FILE'" \
                                -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.txt &
                            else
                                # ---- run simulation from the start if the file DOES NOT exist or reset is true ----
                                echo "Missing pop file: $POP_FILE"

                                bin/slim5.0 -d psi=$psi \
                                -d ID="${rep}" \
                                -d n=$n \
                                -d fluc_opt=$sz \
                                -d fluc_interval=$gap \
                                -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                                -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                                -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                -d "pop_output='$POP_FILE'" \
                                -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                                ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.txt &
                                # alternative code here
                            fi
                            
                            
                            ((c++))
                            if (( c > 29)); 
                            then
                                wait
                                c=0
                            fi
                        done
                    done
                done
            done
    done
done
