#!/bin/bash

dir_name="fluc_env"

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
mkdir -p ../../data/${dir_name}/copied_val
mkdir -p ../../data/${dir_name}/group_tag

python make_model.py 
cp ave.slim ../../data/${dir_name}/
cp ext.slim ../../data/${dir_name}/
cp fit.slim ../../data/${dir_name}/

reset=false # reset the simulation even if the .pop file exists
stopIfPopfileNotFound=false # stop the simulation if the .pop file is not found. If false, it will run the simulation from the start if the .pop file is not found.
c=0
for rep in {1..20}; do  
    for gap in 90; do
        for sz in 10 20 40 60 70 80 90 100; do   
                for n in 8; do
                    for psi in 0.0 0.2 0.5 0.8; do
                        for reg in ave fit ext; do
                            echo "Submitting job, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}, reg=${reg}, gap=${gap}."

                            POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.pop"
                             OUT_POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.pop"
                            if [[  "$reset" == false && -f "$POP_FILE" ]]; then
                                # ---- read the pop file if the pop file EXISTS and reset is false ----
                                echo "Found pop file: $POP_FILE"

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
                                -d "burnin_fn='$POP_FILE'" \
                                -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                -d "allele_fn='../../data/${dir_name}/allele/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}'" \
                                -d "group_tag_fn='../../data/${dir_name}/group_tag/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                -d "copied_val_fn='../../data/${dir_name}/copied_val/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.tsv'" \
                                ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.txt &
                            else
                                # ---- run simulation from the start if the file DOES NOT exist or reset is true ----
                                echo "Missing pop file: $POP_FILE"
                                if [[ "$stopIfPopfileNotFound" == true ]]; then
                                    echo "Stopping simulation because the pop file was not found and stopIfPopfileNotFound is set to true."
                                    exit 1
                                fi

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
                            if (( c > 19)); 
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
