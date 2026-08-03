#!/bin/bash

dir_name="fluc_env_invasion_ave"

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
for rep in {5..20}; do  
    # for gap in 10 20 30 50 60 70 80 90; do
    #     for sz in 10 20 40 60 70 80 90 100; do   
    for gap in 80 90; do
        for sz in 80 100; do  
            # for n in 2 5 8 20; do
            for n in 8; do
                for psi in 0.2 0.5 0.8; do
                    for reg in ave; do
                        for invasion_gen in 0; do
                            for invasion_ind in 500; do
                                echo "Submitting job, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}, reg=${reg}, gap=${gap}, in_Gen=${invasion_gen}, inInd=${invasion_ind}."

                                if (( invasion_gen >= gap )); then
                                    continue
                                fi
                                shift_offset=0
                                base=$((9900 + shift_offset))
                                current_gen_in_cycle=$(( (base - 50) % gap ))
                                # if [[ $invasion_gen == 0 ]]; then
                                    invasion_generation=$(( base + (gap - current_gen_in_cycle) ))
                                # else
                                #     # if (( current_gen_in_cycle < invasion_gen )); then
                                #     #     invasion_generation=$(( base + invasion_gen - current_gen_in_cycle ))
                                #     # else
                                #     #     invasion_generation=$(( base + gap - current_gen_in_cycle + invasion_gen ))
                                #     # fi
                                # fi
                                echo "$invasion_generation"

                                prefix="n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_inGen_${invasion_gen}_absInGen_${invasion_generation}_inInd_${invasion_ind}_${rep}"
                                POP_FILE="../../data/fluc_env/pop/n_${n}_psi_0.0_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.pop"
                                OUT_POP_FILE="../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_${reg}_gap_${gap}_${rep}.pop"

                                if [[  "$reset" == false && -f "$POP_FILE" ]]; then
                                    # ---- read the pop file if the pop file EXISTS and reset is false ----
                                    echo "Found pop file: $POP_FILE"
                                    echo "output prefix: $prefix"


                                    bin/slim5.0 -d psi=$psi \
                                    -d ID="${rep}" \
                                    -d n=$n \
                                    -d fluc_opt=$sz \
                                    -d fluc_interval=$gap \
                                    -d invasion_gen=$invasion_generation \
                                    -d invasion_ind=$invasion_ind \
                                    -d shift_offset=$shift_offset \
                                    -d log_start_gen=9901 \
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
                                    -d fluc_opt=$sz \
                                    -d fluc_interval=$gap \
                                    -d "phenotype_fn='../../data/${dir_name}/phenotype/${prefix}.tsv'" \
                                    -d "genotype_fn='../../data/${dir_name}/genotype/${prefix}.tsv'" \
                                    -d "fre_output='../../data/${dir_name}/data/${prefix}.tsv'" \
                                    -d "ind_fn='../../data/${dir_name}/ind/${prefix}'" \
                                    -d "mut_fn='../../data/${dir_name}/mut/${prefix}'" \
                                    -d "vars_fn='../../data/${dir_name}/vars/${prefix}.tsv'" \
                                    -d "pop_output='$OUT_POP_FILE'" \
                                    -d "fitness_fn='../../data/${dir_name}/fitness/${prefix}.tsv'" \
                                    -d "allele_fn='../../data/${dir_name}/allele/${prefix}'" \
                                    -d "group_tag_fn='../../data/${dir_name}/group_tag/${prefix}.tsv'" \
                                    -d "copied_val_fn='../../data/${dir_name}/copied_val/${prefix}.tsv'" \
                                    ${reg}.slim &> ../../data/${dir_name}/log/${prefix}.txt &
                                    # alternative code here
                                fi
                                
                                
                                ((c++))
                                if (( c > 15)); 
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
done
