
# USAGE: scripts/hom_simhyb_alleles.R

library(tidyverse)

alleles <- read_tsv("het_hom_alleles_simhyb_pileup.txt", col_names = T)
variants = read_tsv("../1_high_cov_variants/Dataset_S1.tsv", col_names = T)
hom_alleles <- left_join(variants, alleles)
hom_alleles_with_AB <- left_join(hom_alleles, alleles)
write_tsv(hom_alleles_with_AB, "hom_alleles_simhyb_pileup.txt", col_names = T)