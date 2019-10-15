#!/usr/env Rscript

# This script takes in low pass dosage data for an entire population
# Using the distribution of relative coverage values in each bin,
# This script applies the mean shift algorithm on the 1-D distribution
# of relative coverage values for each bin to determine the number
# of clusters in each bin.

# Relative coverage values for each line are saved in their own file.
# Because of this, run this script in the same folder as all of the
# low coverage tables you wish to analyze.

# Libraries. meanShift masks dplyr::filter

library(MeanShift)
library(parallel)
library(tidyverse)

# Define function that will run msClustering on each bin
run_MSclustering <- function(x) {
  # clusters <- msClustering(x, kernel = "epanechnikovKernel", h = quantile(dist(t(x)), 0.48)) # changed from h = 0.5 to 0.48 on 27 December 2018, was underfitting some bins
  clusters <- msClustering(x, kernel = "epanechnikovKernel", h = quantile(dist(t(x)), 0.5)) # need to compare output of both
  # clusters <- msClustering(x, kernel = "epanechnikovKernel", h = quantile(dist(t(x)), 0.4))
  clusters$sample <- colnames(x)
  geno.df <- data.frame("sample" = clusters$sample,
                        "dosage.gt" = clusters$labels,
                        "cluster" = clusters$components[1, clusters$labels])
  return(geno.df)
}

# cd to working directory
# setwd("~/Desktop/Comai_Lab/github-repositories/kirk-LOP-all/experiments/2_low_cov_dosage/exploratory_analyses/") # local
# setwd("/share/comailab/kramundson/snakemake_chrom_dosage/data/mapQ-filter") # buffy

# Specify bad sequencing libraries/runs to withhold from analysis
bad_samples <- c("2x_LOP868_279-windowcov.tsv")

# Read in data
files <- dir(pattern = "2x_LOP.+-windowcov.tsv", path = "data/plots/", full.names = T)
print(files)

# Using vector of file names, read each file in, intially as a list with
# map. This will be flattened to a data frame later
all_dosage <- files %>% 
  map(read_delim, delim = " ", na = "NA")

# add sample names to each data frame in all_dosage
for (i in seq(1:length(files))) {all_dosage[[i]]$sample <- files[i]}

# flatten data frame. Also filter out non-chromosome scaffold sequences and bad samples
all_dosage <- do.call(rbind.data.frame, all_dosage) %>% 
  filter(chrom %in% c(paste0("chr0", 1:9), "chr10", "chr11", "chr12", NA)) %>% 
  filter(!(sample %in% bad_samples)) %>% 
  mutate(chrbin = floor(start / end[1])) %>%
  mutate(sample = str_replace(sample, "-windowcov.tsv", "")) %>% 
  filter(!(is.na(chrom))) # version-dependent, this line isn't necessary on my laptop, but not including it on the cluster results in crash

# change stuffer columns to NA
all_dosage[which(is.na(all_dosage$chrom)), ] <- NA

# initiate new list to store p number of 2 x q matrices
# where p is the number of nonoverlapping bins and q is the number of samples analysed
dosage.mtxs <- list()

# populate list with matrices. These are formatted so that MSclustering can work on them
print("populating list of matrices")
for (i in unique(all_dosage$bin)) {
  chrom <- unique(all_dosage$chrom[which(all_dosage$bin == i)])
  start <- unique(all_dosage$start[which(all_dosage$bin == i)])
  end <- unique(all_dosage$end[which(all_dosage$bin == i)])
  bin.pointer <- paste(chrom, start, end, sep = "_")
  sub <- subset(all_dosage, bin == i)$normcov # needs more quantile filter!
  mtx <- t(matrix(c(sub, rep(0, length(sub))), ncol = 2, nrow = length(sub)))
  colnames(mtx) <- unique(filter(all_dosage, bin == i)$sample)
  dosage.mtxs[[bin.pointer]] <- mtx
}

# This is the line that does the thing.
# Throw more cores at it if you dare.
# Caveat don't run it in RStudio, mclapply does not account for the CPU required to run the GUI
# test <- dosage.mtxs[1:2]
# ms.out <- lapply(test, function(x) run_MSclustering(x))
print("Running mean shift clustering")
ms.out <- mclapply(dosage.mtxs, function(x) run_MSclustering(x),
mc.cores = 12, mc.cleanup = T, mc.preschedule = T)

# transfer colnames from original matrices to new list
for (i in seq(1:length(ms.out))) {ms.out[[i]]$bin <- names(dosage.mtxs)[i]}

# flatten new list into data frame
ms.flattened <- do.call(rbind.data.frame, ms.out) %>%
  separate(bin, into = c("chrom", "start", "end"), sep = "_", convert = T)

# inspect the row names afterwards as a sanity check; they preserve chrom start end information

# write data frame to file
write_tsv(ms.flattened, "2019_0107_LOP_250k_dosage_genotypes_50percband.tsv", na = "NA",
          col_names = T)

# output kernel density plot of each bin with inferred peaks, experimental code chunk added 12/27/2018
for (i in unique(all_dosage$bin)) {
  print(i)
  p <- ggplot(filter(all_dosage, bin == i), aes(x = normcov)) +
    geom_density() +
    geom_vline(xintercept = unique(ms.out[[i]]$cluster))
  ggsave(paste0(i, "_bindist.pdf"), plot = p, device = "pdf")
}

print("Done")
# take us home snoop
# library(BRRR)
# skrrrahh("snoop")
