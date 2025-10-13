This folder has four R Markdown files and three data files:

R Markdown files: plotting_func.Rmd dist_plot.Rmd, geno_plot.Rmd, pheno_plot.Rmd

Data files: data_summary.rds, genotype_data.rds, phenotype_data.rds

Each R Markdown does one thing:

dist_plot.Rmd: plots distance to optimum overtime

geno_plot.Rmd: plots genotypic distributions over time / single generation

pheno_plot.Rmd: plots phenotypic distributions over time / single generation

What you need to do to generate plots

The only thing to change is the PARAMETER block inside each R Markdown.

Steps

1, Run plotting_func.Rmd first to load all plotting functions.

2, Open the markdown file and run those plotting parameters setting and data loading block

3, Change the PARAMETER blocks if there's a subset of data to be plotted

4, Run the RUN PLOTTING block. The ouput plot will be inside plots folder

Note:
1, it's possible the bin range doesn't cover the whole range of data and cause an error. Increase the num_bin will fix it.

2, genotype and phenotype data from gen 10000-10100 only records every 10 generation.

3, there's no phenotype data for inverse-variance model.

4, generation 9999 is the last generation before selection.