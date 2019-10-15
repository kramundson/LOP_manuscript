# Experiment 4

## Summary

This Snakemake workflow takes short paired-end reads (gzip fastq format) as input and
produces the following:

 * %HI allele in non-overlapping 1Mb bins for 167 dihaploids of LOP population
 * %HI allele in intervals identified as outlier dosage variants in Experiment 2

Steps:

1. Read trim (cutadapt 1.15)
2. Align reads to reference (bwa mem 0.7.12-r1039)
3. Remove PCR duplicates (Picard 2.18.27-0)
4. Merge libraries by biological sample (samtools 1.9) 
5. Filter out ≤ MQ10 reads (samtools 1.9)
6. Pileup file generation and parsing (custom Python scripts in ```scripts```, from Meric Lieberman)
7. Binned genotyping on non-overlapping 1Mb bins of the reference genome (custom Python scripts in ```scripts```, from Meric Lieberman)
8. Binned genotyping on dosage variant intervals (custom Python scripts in ```scripts```, from Meric Lieberman and Isabelle Henry)

Output:
 * 1Mb_bin_alleles_LOP.txt
 * indel_bin_alleles_LOP.txt
 
To run, cd to this folder and run the following commands:

```
cd ../4_low_cov_snp
snakemake --cores 24 > smk.out 2> smk.err
```

Manuscript figures from these datasets are made in ```results/2019_1012_LOP_figures.Rmd```
