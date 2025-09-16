#!/bin/bash







#TODO: add pop output to future scripts













ln -s ~/bin/slim5.0 bin/slim5.0

dir_name="n_group_translational"

mkdir -p ../../data/${dir_name}
mkdir -p ../../data/${dir_name}/phenotype
mkdir -p ../../data/${dir_name}/genotype
mkdir -p ../../data/${dir_name}/data
mkdir -p ../../data/${dir_name}/ind
mkdir -p ../../data/${dir_name}/mut
mkdir -p ../../data/${dir_name}/log
mkdir -p ../../data/${dir_name}/pop

# for rep in {9..10}; do
#     for n in 2 3 8; do
#         for sz in 0 20 60 100; do
#             echo "Running simulation with n=${n}"
#             for psi in 0.0 0.3 0.5 0.7 0.9; do
#                 echo "Submitting job for psi=${psi}, n=${n}, sz=${sz}, rep=${rep}."
#                 # bash no_burnin.sh ${psi} ${n} &

#                 bin/slim5.0 -d psi=$psi \
#                     -d ID="${rep}" \
#                     -d n=$n \
#                     -d shift_size=$sz \
#                     -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.tsv'" \
#                     n_group_new.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_ave_${rep}.txt &
#             done
#         done
#     done
#     wait
# done

for rep in {1..10}; do
    for n in 2 3 8; do
        for sz in 0 20 60 100; do
            echo "Running simulation with n=${n}"
            for psi in 0.0 0.3 0.5 0.7 0.9; do
                echo "Submitting job for psi=${psi}, n=${n}, sz=${sz}, rep=${rep}."
                # bash no_burnin.sh ${psi} ${n} &

                bin/slim5.0 -d psi=$psi \
                    -d ID="${rep}" \
                    -d n=$n \
                    -d shift_size=$sz \
                    -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "fre_output='../../data/${dir_name}/data/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "ind_fn='../../data/${dir_name}/ind/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "mut_fn='../../data/${dir_name}/mut/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.tsv'" \
                    -d "pop_output='../../data/${dir_name}/debug_n_8_psi_0.9.pop'" \
                    n_group_maxfit.slim &> ../../data/${dir_name}/log/n_${n}_psi_${psi}_sz_${sz}_reg_fit_${rep}.txt &
            done
        done
    done
    wait
done
