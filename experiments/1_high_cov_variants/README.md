# Experiment 1

## Summary

This Snakemake workflow takes short paired end reads (gzip fastq format) as input and
produces the following:

 * Raw SNP and indel variants (VCF format)
 * Filtered SNP Datasets S1 and S2 used for analyses in manuscript (tsv format)
 * Filtered SNP Dataset S3 used for inference of SNP event run lengths (tsv format)

Steps:
1. Download reads from NCBI SRA (sra-tools 2.8.2)
2. Read trim (cutadapt 1.15)
3. Align reads to reference (bwa mem 0.7.12-r1039)
4. Remove PCR duplicates (Picard 2.18.27-0)
5. Filter out pairs with mates on different chromosomes (awk)
6. Soft clip one mate of an overlapping mate pair (bamUtil 1.0.14)
7. Merge libraries by biological sample (samtools 1.9)
8. Define indel realign targets and realign indels (GATK 3.8)
9. Call variants (freebayes 1.2.0). Variant calling is scattered across scaffolds.
10. Collect read depth values at all genomic positions (samtools 1.9)
11. Summarize median read depth in non-overlapping 10kb windows (custom Python and R scripts available in ```scripts/```

Output:
 * Dataset_S1.tsv
 * Dataset_S2.tsv
 * Dataset_S3.tsv
 * Fig. S7
 * Fig. S8
 * Fig. S9
 
To run, cd to this folder and enter the following:

```
# find scaffold breaks in reference genome
snakemake -s 1_init_files.snakes > 1_smk.out 2> 1_smk.err

# do all read processing, variant calling, and output generation using 64 cores.
# number of cores can be adjusted as necessary
snakemake -s 2_fastq_vcf.snakes --cores 64 > 2_smk.out 2> 2_smk.err
```