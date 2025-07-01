
for i in {0..5}
do
    slim -d ID=$i -d "fn='tsv_output/social_new/data/sburnin-2702579/pop/rec_N_1000_psi_0.0_pop_2702577.pop'" -d psi11=0.0 -d "phenotype_fn='tsv_output/social_new/data/phenotype/rec_N_1000_psi_0.0_pop_270257_${i}.tsv'" after_burnin_rec.slim &

    slim -d ID=$i -d "fn='tsv_output/social_new/data/sburnin-2702579/pop/rec_N_1000_psi_0.0_pop_2702581.pop'" -d psi11=0.0 -d "phenotype_fn='tsv_output/social_new/data/phenotype/rec_N_1000_psi_0.0_pop_2702581_${i}.tsv'" after_burnin_rec.slim &

    slim -d ID=$i -d "fn='tsv_output/social_new/data/sburnin-2702579/pop/nrec_N_1000_psi_0.0_pop_2702577.pop'" -d psi11=0.0 -d "phenotype_fn='tsv_output/social_new/data/phenotype/nrec_N_1000_psi_0.0_pop_270257_${i}.tsv'" after_burnin_nonrec.slim &

    slim -d ID=$i -d "fn='tsv_output/social_new/data/sburnin-2702579/pop/nrec_N_1000_psi_0.0_pop_2702581.pop'" -d psi11=0.0 -d "phenotype_fn='tsv_output/social_new/data/phenotype/nrec_N_1000_psi_0.0_pop_2702581_${i}.tsv'" after_burnin_nonrec.slim &

    wait
done


# slim -d "fn='tsv_output/social_new/data/sburnin-2702580/pop/rec_N_1000_psi_0.0_pop_2702578.pop'" -d psi11=0.0 load_burnin.slim > tsv_output/social_new/data/phenotype/rec_N_1000_psi_0.0_pop_2702578.tsv

# slim -d "fn='tsv_output/social_new/data/sburnin-2702580/pop/rec_N_1000_psi_0.0_pop_2702582.pop'" -d psi11=0.0 load_burnin.slim > tsv_output/social_new/data/phenotype/rec_N_1000_psi_0.0_pop_2702582.tsv


# slim -d "fn='tsv_output/social_new/data/sburnin-2702581/pop/rec_N_1000_psi_0.0_pop_2702579.pop'" -d psi11=0.0 load_burnin.slim > tsv_output/social_new/data/phenotype/rec_N_1000_psi_0.0_pop_2702579.tsv

# slim -d "fn='tsv_output/social_new/data/sburnin-2702581/pop/rec_N_1000_psi_0.0_pop_2702583.pop'" -d psi11=0.0 load_burnin.slim > tsv_output/social_new/data/phenotype/rec_N_1000_psi_0.0_pop_2702583.tsv
