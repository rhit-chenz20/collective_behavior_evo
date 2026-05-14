#!/bin/bash

dir_name="various_group_size"

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
mkdir -p ../../data/${dir_name}/group_tags

python make_model.py 

reset=false # reset the simulation even if the .pop file exists
stopIfPopfileNotFound=false # stop the simulation if the .pop file is not found. If false, it will run the simulation from the start if the .pop file is not found.
c=0
# burn in using sz=20
for rep in {1..5}; do  
   for sz in 40; do   
        for n in 5 8 20; do
            for psi in 0.0 0.2 0.5 0.8; do
                for reg in ave fit; do
                    echo "Submitting job, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}, reg=${reg}."
                
                    POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_${psi}_reg_${reg}_${rep}.pop"
                    # POP_FILE="../../data/constant_env_batch/pop/n_${n}_psi_${psi}_reg_${reg}_${rep}.pop"
                    # POP_FILE1="../../data/constant_env_batch/pop/n_${n}_psi_${psi}_sz_20_reg_${reg}_${rep}.pop"

                    # OUT_POP_FILE=POP_FILE
                    OUT_POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_${psi}_reg_${reg}_${rep}.pop"
                    if [[ "$reset" == false && ( -f "$POP_FILE" || -f "$POP_FILE1" ) ]]; then
                        # ---- read the pop file if the pop file EXISTS and reset is false ----

                        if [[ -f "$POP_FILE" ]]; then
                            POP_FILE="$POP_FILE"
                        elif [[ -f "$POP_FILE1" ]]; then
                            POP_FILE="$POP_FILE1"
                        fi
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
                        -d "copied_val_fn='../../data/${dir_name}/copied_val/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "group_tag_fn='../../data/${dir_name}/group_tags/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.txt &

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
                        -d shift_size=$sz \
                        -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                        -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                        -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "pop_output='$OUT_POP_FILE'" \
                        -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "allele_fn='../../data/${dir_name}/allele/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
                        -d "copied_val_fn='../../data/${dir_name}/copied_val/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        -d "group_tag_fn='../../data/${dir_name}/group_tags/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
                        ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.txt &

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


# # reusing burn in above to do other shift size
# for rep in {1..10}; do  
#    for sz in 40 60; do   
#         for n in 2 5 8 20; do
#             for psi in 0.1 0.2 0.3 0.8; do
#                 for reg in ext ave fit; do
#                     echo "Submitting job, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}, reg=${reg}."
                
#                     POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_${psi}_reg_${reg}_${rep}.pop"
#                     # POP_FILE="../../data/constant_env_ext/pop/n_${n}_psi_${psi}_sz_20_reg_${reg}_${rep}.pop"

#                     if [[  "$reset" == false && -f "$POP_FILE" ]]; then
#                         # ---- read the pop file if the pop file EXISTS and reset is false ----
#                         echo "Found pop file: $POP_FILE"

#                         bin/slim5.0 -d psi=$psi \
#                         -d ID="${rep}" \
#                         -d n=$n \
#                         -d shift_size=$sz \
#                         -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
#                         -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
#                         -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         -d "burnin_fn='$POP_FILE'" \
#                         -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         -d "allele_fn='../../data/${dir_name}/allele/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
#                         ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.txt &

#                     else
#                         # ---- run simulation from the start if the file DOES NOT exist or reset is true ----
#                         echo "Missing pop file for shift size $sz: $POP_FILE"

#                         # bin/slim5.0 -d psi=$psi \
#                         # -d ID="${rep}" \
#                         # -d n=$n \
#                         # -d shift_size=$sz \
#                         # -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         # -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         # -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         # -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
#                         # -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
#                         # -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         # -d "pop_output='$POP_FILE'" \
#                         # -d "fitness_fn='../../data/${dir_name}/fitness/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.tsv'" \
#                         # -d "allele_fn='../../data/${dir_name}/allele/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}'" \
#                         # ${reg}.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_${rep}.txt &

#                     fi
#                     ((c++))
#                     if (( c > 29)); 
#                     then
#                         wait
#                         c=0
#                     fi
#                 done
#             done
#         done
#     done
# done