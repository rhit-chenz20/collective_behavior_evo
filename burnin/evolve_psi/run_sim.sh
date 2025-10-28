#!/bin/bash

dir_name="evolve_psi"

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

# for rep in 1; do
#     for n in 2; do
#         for sz in 40; do
#             echo "Running simulation with n=${n}"
#             for psi in 0.1; do
#                 echo "Submitting job for fit, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}."

#                 bin/slim5.0 -d psi=$psi \
#                     -d beta=-100 \
#                     -d ID="${rep}" \
#                     -d n=$n \
#                     -d shift_size=$sz \
#                     -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "pop_output='../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.pop'" \
#                     -d "psi_fn='../../data/${dir_name}/psis/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     n_group_ave.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.txt &
#             done
#         done
#     done
# done
c=0
for rep in {1..10}; do  
   for sz in 40 80 120 300; do   
        for n in 2 8 20 50; do
            for psi in 0.0; do
                echo "Submitting job for fit, psi=${psi}, n=${n}, sz=${sz}, rep=${rep}."

                bin/slim5.0 -d psi=$psi \
                    -d beta=-100 \
                    -d ID="${rep}" \
                    -d n=$n \
                    -d shift_size=$sz \
                    -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "burnin_fn='../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.pop'" \
                    -d "psi_fn='../../data/${dir_name}/psis/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "gpsi_fn='../../data/${dir_name}/group_psis/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    n_group_maxfit.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.txt &
                ((c++))
# -d "pop_output='../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.pop'" \



                # bin/slim5.0 -d psi=$psi \
                #     -d beta=-100 \
                #     -d ID="${rep}" \
                #     -d n=$n \
                #     -d shift_size=$sz \
                #     -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
                #     -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
                #     -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
                #     -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
                #     -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
                #     -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
                #     -d "pop_output='../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.pop'" \
                #     -d "psi_fn='../../data/${dir_name}/psis/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
                #     -d "gpsi_fn='../../data/${dir_name}/group_psis/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
                #     n_group_ave.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.txt &
                # ((c++))
                # bin/slim5.0 -d psi=$psi \
                #     -d beta=-100 \
                #     -d ID="${rep}" \
                #     -d n=$n \
                #     -d shift_size=$sz \
                #     -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.tsv'" \
                #     -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.tsv'" \
                #     -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.tsv'" \
                #     -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.tsv'" \
                #     -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.tsv'" \
                #     -d "vars_fn='../../data/${dir_name}/vars/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.tsv'" \
                #     -d "pop_output='../../data/${dir_name}/pop/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.pop'" \
                #     -d "psi_fn='../../data/${dir_name}/psis/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.tsv'" \
                #     -d "gpsi_fn='../../data/${dir_name}/group_psis/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.tsv'" \
                #     n_group_extreme.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_ext_${rep}.txt &
                
                # ((c++))
                if (( c > 29)); 
                then
                    wait
                    c=0
                fi
            done
        done
    done
done
