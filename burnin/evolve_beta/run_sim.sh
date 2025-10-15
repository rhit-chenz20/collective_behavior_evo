#!/bin/bash

# ln -s ../bin/slim5.0 bin/slim5.0

dir_name="evolve_beta"

mkdir -p ../../data/${dir_name}
mkdir -p ../../data/${dir_name}/phenotype
mkdir -p ../../data/${dir_name}/genotype
mkdir -p ../../data/${dir_name}/data
mkdir -p ../../data/${dir_name}/ind
mkdir -p ../../data/${dir_name}/mut
mkdir -p ../../data/${dir_name}/log
mkdir -p ../../data/${dir_name}/pop
mkdir -p ../../data/${dir_name}/vars
mkdir -p ../../data/${dir_name}/betas
mkdir -p ../../data/${dir_name}/group_betas

c=0
for rep in {1..10}; do  
   for sz in 40 80 120 300; do   
        for n in 2 8 20 50; do
            for beta in 0.0 0.05 0.5 1.0; do
                echo "Submitting job for beta=${beta}, n=${n}, sz=${sz}, rep=${rep}."

                bin/slim5.0 -d psi=-100 \
                    -d beta=$beta \
                    -d ID="${rep}" \
                    -d n=$n \
                    -d shift_size=$sz \
                    -d "phenotype_fn='../../data/${dir_name}/phenotype/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.tsv'" \
                    -d "genotype_fn='../../data/${dir_name}/genotype/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.tsv'" \
                    -d "fre_output='../../data/${dir_name}/data/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.tsv'" \
                    -d "ind_fn='../../data/${dir_name}/ind/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.tsv'" \
                    -d "mut_fn='../../data/${dir_name}/mut/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.tsv'" \
                    -d "vars_fn='../../data/${dir_name}/vars/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.tsv'" \
                    -d "pop_output='../../data/${dir_name}/pop/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.pop'" \
                    -d "beta_fn='../../data/${dir_name}/betas/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.tsv'" \
                    -d "gbeta_fn='../../data/${dir_name}/group_betas/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.tsv'" \
                    inverse.slim &> ../../data/${dir_name}/log/n_${n}_beta_${beta}_sz_${sz}_reg_inv_${rep}.txt &

                ((c++))
                if (( c > 10)); then
                    wait
                    c=0
                fi
            done
        done
    done
done

