

# for i in $(seq 1 10);
# do
#     slim  -d ID=$i 1_trait.slim > log/log${i}.txt &
# done

CONFIGFILES=$(ls slim/simulation-parameters/*.yaml)
NUM_SIMS=5

for file in ${CONFIGFILES}
do
    echo ${file}
    snakemake --use-conda --snakefile simulation.smk --configfile ${file} --config slim=slim use_subdirectory=true 
done
