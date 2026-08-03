#!/bin/bash

dir_name="cons_env_invasion"

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
stopIfPopfileNotFound=true # stop the simulation if the .pop file is not found. If false, it will run the simulation from the start if the .pop file is not found.
c=0
# burn in using sz=20
for rep in {1..20}; do  
   for sz in 40; do   
        for n in 2 5 8 20; do
            for psi in 0.2 0.5 0.8; do
                for reg in ext ave fit; do
                    for invasion_gen in 9000 9500 9900 10000 10020 10050 10100; do
                        for invasion_ind in 1 10 50 100; do
                            
                            echo "Submitting job, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}, reg=${reg}."

                            prefix="n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_inGen_${invasion_gen}_inInd_${invasion_ind}_${rep}"
                    
                            POP_FILE="../../data/cons_env/pop/n_${n}_psi_0.0_reg_${reg}_${rep}.pop"
                            # OUT_POP_FILE=POP_FILE
                            OUT_POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_0.0_reg_${reg}_${rep}.pop"

                            shift_offset=1000
                            if [[ "$reset" == false && ( -f "$POP_FILE" ) ]]; then
                                echo "Found pop file: $POP_FILE"
                                echo "output prefix: $prefix"

                                # if [[ (! -f "$OUT_POP_FILE" ) ]]; then
                                #     echo "pop file does not exist in current folder. copying $OUT_POP_FILE"
                                #     cp $POP_FILE $OUT_POP_FILE
                                # fi
                                
                                bin/slim5.0 -d psi=$psi \
                                -d ID="${rep}" \
                                -d n=$n \
                                -d shift_size=$sz \
                                -d invasion_gen=$(($invasion_gen + $shift_offset)) \
                                -d invasion_ind=$invasion_ind \
                                -d shift_offset=$shift_offset \
                                -d log_start_gen=$(($invasion_gen + $shift_offset - 30)) \
                                -d log_pop_gen=9000 \
                                -d "phenotype_fn='../../data/${dir_name}/phenotype/${prefix}.tsv'" \
                                -d "genotype_fn='../../data/${dir_name}/genotype/${prefix}.tsv'" \
                                -d "fre_output='../../data/${dir_name}/data/${prefix}.tsv'" \
                                -d "ind_fn='../../data/${dir_name}/ind/${prefix}'" \
                                -d "mut_fn='../../data/${dir_name}/mut/${prefix}'" \
                                -d "vars_fn='../../data/${dir_name}/vars/${prefix}.tsv'" \
                                -d "burnin_fn='$POP_FILE'" \
                                -d "fitness_fn='../../data/${dir_name}/fitness/${prefix}.tsv'" \
                                -d "allele_fn='../../data/${dir_name}/allele/${prefix}'" \
                                -d "group_tag_fn='../../data/${dir_name}/group_tag/${prefix}.tsv'" \
                                -d "copied_val_fn='../../data/${dir_name}/copied_val/${prefix}.tsv'" \
                                ${reg}.slim &> ../../data/${dir_name}/log/${prefix}.txt &

                                # mv $POP_FILE ../../data/${dir_name}/pop/
                            else
                                # ---- run simulation from the start if the file DOES NOT exist or reset is true ----
                                echo "Missing pop file: $POP_FILE or reset is true. Running simulation from the start."
                                if [[ "$stopIfPopfileNotFound" == true ]]; then
                                    echo "Stopping simulation because the pop file was not found and stopIfPopfileNotFound is set to true."
                                    continue
                                fi
                                bin/slim5.0 -d psi=$psi \
                                -d ID="${rep}" \
                                -d n=$n \
                                -d shift_size=$sz \
                                -d invasion_gen=$invasion_gen \
                                -d invasion_ind=$invasion_ind \
                                -d "phenotype_fn='../../data/${dir_name}/phenotype/${prefix}.tsv'" \
                                -d "genotype_fn='../../data/${dir_name}/genotype/${prefix}.tsv'" \
                                -d "fre_output='../../data/${dir_name}/data/${prefix}.tsv'" \
                                -d "ind_fn='../../data/${dir_name}/ind/${prefix}'" \
                                -d "mut_fn='../../data/${dir_name}/mut/${prefix}'" \
                                -d "vars_fn='../../data/${dir_name}/vars/${prefix}.tsv'" \
                                -d "pop_output='$OUT_POP_FILE'" \
                                -d "fitness_fn='../../data/${dir_name}/fitness/${prefix}.tsv'" \
                                -d "group_tag_fn='../../data/${dir_name}/group_tag/${prefix}.tsv'" \
                                -d "allele_fn='../../data/${dir_name}/allele/${prefix}'" \
                                -d "copied_val_fn='../../data/${dir_name}/copied_val/${prefix}.tsv'" \
                                ${reg}.slim &> ../../data/${dir_name}/log/${prefix}.txt &

                            fi
                            ((c++))
                            if (( c > 40)); 
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