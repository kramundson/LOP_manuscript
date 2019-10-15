#!/Users/Kirk/miniconda3/envs/ximena/bin/Rscript
# Kirk Amundson
# dosage_outliers.R

library(tidyverse)
library(MeanShift)
library(plotly)
library(ggbeeswarm)
library(gridExtra)

gts <- read_tsv("2018_1018_LOP_250k_dosage_genotypes_50percband.tsv", col_names = T)

# based on consistent intermediate dosage profiles, I think these samples should be removed from the analysis altogether before
# considering what is commmon and what isn't in the dosage analysis.
likely_tetraploids <- c("097",  "269",  "280",  "466",  "509",  "1157", "1180", "1192", "1190")

gts %>% 
  group_by(chrom, start, end, cluster) %>% 
  summarize(cluster_count = n()) %>% 
  left_join(gts, .) %>% 
  mutate(sample = gsub("2x_LOP868_", "", .$sample)) %>% 
  filter(!sample %in% likely_tetraploids) %>% 
  filter(chrom == "chr06" & between(start, 3.175e7, 3.4e7)) %>% 
  filter(cluster_count >= 167 * 0.05) %>% 
  ggplot(., aes())

# add number of individuals in clusters
outliers <- gts %>% 
  group_by(chrom,start,end,cluster) %>% 
  summarize(cluster_count = n()) %>% 
  left_join(gts, .) %>% 
  mutate(rare = ifelse(cluster_count / length(unique(.$sample)) <= 0.05, TRUE, FALSE)) %>% 
  mutate(sample = gsub("2x_LOP868_", "", .$sample)) %>% 
  filter(!sample %in% likely_tetraploids) %>%
  filter(!(sample == "238" & chrom == "chr01")) %>% 
  filter(!(sample %in% c("259", "292", "452", "460", "1157", "1190") & chrom == "chr02")) %>% 
  filter(!(sample == "373" & chrom == "chr03")) %>%
  filter(!(sample %in% c("028","433") & chrom == "chr04")) %>% 
  filter(!(sample %in% c("293", "363") & chrom == "chr05")) %>%
  filter(!(sample == "172" & chrom == "chr06")) %>%
  filter(!(sample %in% c("065", "272") & chrom == "chr07")) %>%
  filter(!(sample %in% c("108", "488") & chrom == "chr08")) %>%
  filter(!(sample %in% c("397") & chrom == "chr10")) %>%
  filter(!(sample %in% c("128") & chrom == "chr12")) %>%
  filter(rare) %>% 
  dplyr::select(-rare, -cluster_count) %>% 
  mutate(Bin = paste(chrom, as.integer(start), as.integer(end), sep = "_")) %>% 
  arrange(sample, chrom, start, end)

# merge adjacent outlier bins
merged_vars <- data.frame("sample" = rep(NA,nrow(outliers)),
                          "dosage.gt" = rep(NA,nrow(outliers)),
                          "cluster" = rep(NA,nrow(outliers)),
                          "chrom" = rep(NA,nrow(outliers)),
                          "start" = rep(NA,nrow(outliers)),
                          "end" = rep(NA,nrow(outliers)),
                          "Bin" = rep(NA,nrow(outliers)))

clusters <- rep(NA, nrow(outliers))
m <- "off"
row_count <- 1

for (i in 1:(nrow(outliers)-1)) {
  # print(i)
  row1 <- outliers[i,]
  row2 <- outliers[i+1,]
  
  # Is a merge in progress?
  if (m == "off") {
    
    if (row1$sample == row2$sample && row1$chrom == row2$chrom && row1$end == row2$start) {
      m <- "on"
      merged_vars[row_count,] <- row1
      clusters[row_count] <- row1$cluster
    }
    
    # no merge coming, no merge to start
    else {
      merged_vars[row_count,] <- row1
      clusters[row_count] <- row1$cluster
      row_count <- row_count + 1
    }
  }
  
  # merge is in progress
  else {
    if (row1$sample == row2$sample && row1$chrom == row2$chrom && row1$end == row2$start) {
      clusters[row_count] <- paste(clusters[row_count], row1$cluster, sep = "_")
    }
    
    else {
      m <- "off"
      merged_vars[row_count,6] <- row1$end
      # merged_vars[row_count,7] <- row1$end
      clusters[row_count] <- paste(clusters[row_count], row1$cluster, sep = "_")
      row_count <- row_count + 1
    }
  }
}

merged_vars %>% 
  cbind(., clusters) %>% 
  filter(!is.na(sample)) %>%
  mutate(size = end-start) %>% 
  arrange(chrom, start, end, sample) %>%
  dplyr::select(chrom, start, end, size, sample, clusters) %>% 
  mutate(sample = gsub("^", "2x_LOP868_", .$sample)) %>% 
  write_tsv(., "LOP_outlier_vars_all.tsv")

merged_vars %>% 
  filter(!is.na(sample)) %>%
  filter(!sample %in% likely_tetraploids) %>% 
  mutate(size = end-start) %>% 
  filter(size >= 7.5e5) %>% 
  arrange(chrom, start, end, sample) %>% 
  dplyr::select(chrom, start, end, size, sample) %>% 
  mutate(sample = gsub("^", "2x_LOP868_", .$sample)) %>% 
  write_tsv(., "LOP_outlier_vars_to_genoype.tsv")
