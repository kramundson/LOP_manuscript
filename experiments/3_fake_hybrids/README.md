# Experiment 3

## Summary

This Snakemake workflow takes short paired-end reads (gzip fastq format) as input and
produces the following:

 * 100 simulated Alca Tarma x PL4 low-coverage hybrids
 * 100 simulated Alca Tarma x IVP101 low-coverage hybrids
 * 100 simulated matrilineal Alca Tarma low coverage samples
 * %HI allele in non-overlapping 1Mb bins for all 300 samples
 * %HI allele in intervals identified as outlier dosage variants in Experiment 2

Steps:
1. In-memory random sampling of parental read sets (custom Python script in ```scripts```)
2. Read trim (cutadapt 1.15)
3. Align reads to reference (bwa mem 0.7.12-r1039)
4. Remove PCR duplicates (Picard 2.18.27-0)
5. Merge libraries by biological sample (samtools 1.9) 
6. Filter out ≤ MQ10 reads (samtools 1.9)
7. Pileup file generation and parsing (custom Python scripts in ```scripts```, from Meric Lieberman)
8. Binned genotyping on non-overlapping 1Mb bins of the reference genome (custom Python scripts in ```scripts```, from Meric Lieberman)
9. Binned genotyping on dosage variant intervals (custom Python scripts in ```scripts```, from Meric Lieberman and Isabelle Henry)

Output:
 * 1Mb_bin_alleles_simhyb.txt
 * indel_bin_alleles_simhyb.txt
 
To run, cd to this folder and run the following commands:

```
cd ../3_fake_hybrids
snakemake -s 1_make_simhybs.snakes --cores 64 > 1_smk.out 2> 1_smk.err
snakemake -s 2_fastq_to_bam.snakes --cores 64 > 2_smk.out 2> 2_smk.err
```

Manuscript figures from these datasets are made in ```results/2019_1012_LOP_figures.Rmd```