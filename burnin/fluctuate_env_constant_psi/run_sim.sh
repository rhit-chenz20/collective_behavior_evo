#!/bin/bash

dir_name="pop_files"

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
# mkdir -p ../../data/${dir_name}/allele
# mkdir -p ../../data/${dir_name}/copied_val
# mkdir -p ../../data/${dir_name}/group_tag
mkdir -p ../../data/${dir_name}/pop/fluc_pop

python make_model.py 
cp ave.slim ../../data/${dir_name}/fluc_ave.slim
cp ext.slim ../../data/${dir_name}/fluc_ext.slim
cp fit.slim ../../data/${dir_name}/fluc_fit.slim

reset=true # reset the simulation even if the .pop file exists
stopIfPopfileNotFound=false # stop the simulation if the .pop file is not found. If false, it will run the simulation from the start if the .pop file is not found.
c=0

for rep in {1..50}; do  
    for gap in 20 40 60 80 100 120 ; do
        for sz in 20 40 60 80 100; do   
                # for n in 2 4 5 8 10 20 25 40 50 100; do
                for n in 2 5 8 10 20 50 ; do
                    for psi in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
                        for reg in ave fit ext; do
                            echo "Submitting job, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}, reg=${reg}, gap=${gap}."

                            POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.pop"
                            if [[  "$reset" == false && -f "$POP_FILE" ]]; then
                                # ---- read the pop file if the pop file EXISTS and reset is false ----
                                echo "Found pop file: $POP_FILE skipping simulation"
                            else
                                # ---- run simulation from the start if the file DOES NOT exist or reset is true ----
                                echo "Running simulation for pop file: $POP_FILE"
                                OUT_POP_FILE="../../data/${dir_name}/pop/fluc_pop/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.pop"
                            
                                bin/slim5.0 -d psi=$psi \
                                -d ID="${rep}" \
                                -d n=$n \
                                -d fluc_opt=$sz \
                                -d fluc_interval=$gap \
                                -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}'" \
                                -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}'" \
                                -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                -d "pop_output='$OUT_POP_FILE'" \
                                -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                -d "allele_fn='../../data/${dir_name}/allele/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}'" \
                                -d "group_tag_fn='../../data/${dir_name}/group_tag/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                -d "copied_val_fn='../../data/${dir_name}/copied_val/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.txt &
                                # alternative code here
                            fi
                            
                            ((c++))
                            if (( c > 15 )); 
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
