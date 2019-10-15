# Experiment 2

## Summary

This Snakemake workflow takes short single end reads (gzip fastq format) as input and
produces the following:

 * 250kb bin coverage values (estimated copy number state) normalized  to tetraploid Alca Tarma.
   One .tsv per individual.
 * Called dosage states (quasi-genotypes) at each 250kb bin for each dihaploid, tsv format
 * List of outlier variants, i.e., dosage states present at ≤5% allele frequency in pop.
 * List of outlier duplication variants ≥ 750kb to be genotyped in experiments 3 and 4.

Steps:
1. Download reads from NCBI SRA (sra-tools 2.8.2)
2. Read trim (cutadapt 1.15)
3. Align reads to reference (bwa mem 0.7.12-r1039)
4. Remove PCR duplicates (Picard 2.18.27-0)
5. Merge libraries by biological sample (samtools 1.9)
6. Filter out ≤ MQ10 reads (samtools 1.9)
7. Count aligned reads in non-overlapping 250kb bins of the reference genome (bedtools 2.27.1)
8. Compute normalized coverage values for each dihaploid (custom R script available in ```scripts```)
9. Call dosage genotype states from clusters of normalized coverage values (custom R script available in ```scripts```)

Output:
 * 2019_0107_LOP_250k_dosage_genotypes_50percband.tsv
 * LOP_outlier_vars_all.tsv
 * LOP_outlier_vars_to_genotype.tsv"
 
To run, cd to this folder and run the following commands:

```
snakemake --cores 64 > smk.out 2> smk.err
```

Manuscript figures generated in this experiment are generated in ```results/2019_1012_LOP_figures.Rmd```