Workflow will need to encompass:

1. snakemake chrom dosage of LOP dihaploids. done for bins and chromosomes.
2. Clustering dosage variable states to get quasi-genotypes. Currently generates raw genotypes that are not filtered for outliers. I'm going to do this in a different repo.
3. Generates:

    Fig. 3B. Harder to do. May want to make a circos-specific repo.
    Fig. 3G. This would be part of a circos-specific repo as well.
    Coverage plots of Fig. S1. Make once SNP plots are made.
    Fig. S4. Done? Still need to do histogram/KDE part.
    
    
Figures 2A, 3A, S3, 4A-D are all made in the .Rmd notebook ```../../results/2019_1012_LOP_figures.Rmd```.
I am not going to automate these.
