#!/Users/Kirk/miniconda3/envs/ximena/bin/Rscript
# make_simhyb_lesions.R
# Kirk Amundson

# USAGE: Rscript scripts/make_simhyb_lesions.R input_lesions_real output_lesions_simulated_hybs

library(tidyverse)

setwd("Desktop/Comai_Lab/github-repositories/squeaky_clean_LOP/experiments/3_fake_hybrids/")
args <- commandArgs(trailingOnly = F)

lesions <- read_tsv(args[1], col_names = T) %>%
  select(chrom, start, end, size)

simhybs <- dir(path = "data/reads/", pattern = "[0-9]{2,}_SRR6123032_", full.names = F) %>% 
  str_remove(".fastq.gz$") %>% 
  merge(lesions, .) %>% 
  rename(sample = y) %>% 
  arrange(sample, chrom, start, end) %>% 
  distinct() %>% 
  write_tsv(., args[2], col_names = T)