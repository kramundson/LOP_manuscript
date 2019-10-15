#!/share/comailab/kramundson/miniconda3/envs/ximena/bin/Rscript
# Kirk Amundson
# expt1_datasets_figures.R

# USAGE: Rscript expt1_datasets_figures.R

library(tidyverse)
library(egg)
library(stringr) # have to load this in separately from tidyverse for reasons unknown
library(MeanShift)

args <- commandArgs(trailingOnly = T)

filtered <- dir(path = "data/calls/", pattern = "filtered", full.names = T)
vcf_all_filt <- filtered %>%
  map(read_tsv, col_names = T) %>%
  do.call(rbind.data.frame, .)

# # write out supplemental dataset S1, is homozygous SNP for experiment 3 at whole chromosome resolution
vcf_all_filt %>%
  filter(CHROM != "ChrUn") %>%
  filter(PL4_GT != "0/1") %>%
  filter(IVP101_GT != "0/1") %>%
  mutate(HI = ifelse(nonhi_gt == "1/1", REF, ALT)) %>%
  mutate(AlcaTarma = ifelse(nonhi_gt == "1/1", ALT, REF)) %>%
  select(CHROM, POS, REF, HI, AlcaTarma) %>%
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>%
  write_tsv(., "Dataset_S1.tsv", col_names = T) # was 2019_0926_LOP_hom_SNP.tsv

# write out supplemental dataset S2, is homozygous and heterozygous SNP for experiment 3
vcf_all_filt %>%
  filter(CHROM != "ChrUn") %>%
  mutate(HI = ifelse(nonhi_gt == "1/1", REF, ALT)) %>%
  mutate(AlcaTarma = ifelse(nonhi_gt == "1/1", ALT, REF)) %>%
  select(CHROM, POS, REF, HI, AlcaTarma) %>%
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>%
  write_tsv(., "Dataset_S2.tsv", col_names = T)
  # write_tsv(., "2019_0926_LOP_SNP.tsv", col_names = T)

# write out supplemental dataset S3, is combined filtered calls for each dihap
write_tsv(vcf_all_filt, "Dataset_S3.tsv", col_names = T)

# also write out BED3 for histogram counts, for Fig. 3
vcf_all_filt %>%
  filter(!(CHROM %in% c("chr00", "ChrUn"))) %>%
  select(CHROM, POS) %>%
  mutate(start = POS-1) %>%  # vcf is 1-based, changed to 0-based BED3 format here
  mutate(end = POS) %>%
  select(CHROM, start, end) %>%
  write_tsv(., "LOP_SNP.bed", col_names = F)

# # summarize rate of introgression-positive calls in each dihaploid
# # LOP868.004
# summary4 <- vcf_all_filt %>%
#   filter(!is.na(LOP868_004_GT)) %>%
#   filter(LOP868_004_DP >= 10) %>%
#   group_by(LOP868_004_GT) %>%
#   tally
# print(summary4)
#
# # LOP868.004, restricted
# summary4r <- vcf_all_filt %>%
#   filter(!is.na(LOP868_004_GT)) %>%
#   filter(LOP868_004_DP >= 10) %>%
#   filter(LOP868_004_ID >= 5) %>%
#   filter(LOP868_004_ID / LOP868_004_DP >= 0.1) %>%
#   group_by(LOP868_004_GT) %>%
#   tally
# print(summary4r)
#
# # LOP868.064
# summary64 <- vcf_all_filt %>%
#   filter(!is.na(LOP868_064_GT)) %>%
#   filter(LOP868_064_DP >= 10) %>%
#   group_by(LOP868_064_GT) %>%
#   tally
# print(summary64)
#
# # LOP868.064, restricted
# summary64r <- vcf_all_filt %>%
#   filter(!is.na(LOP868_064_GT)) %>%
#   filter(LOP868_064_DP >= 10) %>%
#   filter(LOP868_064_ID >= 5) %>%
#   filter(LOP868_064_ID / LOP868_064_DP >= 0.1) %>%
#   group_by(LOP868_064_GT) %>%
#   tally
# print(summary64r)
#
# # LOP868.305
# summary305 <- vcf_all_filt %>%
#   filter(!is.na(LOP868_305_GT)) %>%
#   filter(LOP868_305_DP >= 10) %>%
#   group_by(LOP868_305_GT) %>%
#   tally
# print(summary305)
#
# # LOP868.305, restricted
# summary305r <- vcf_all_filt %>%
#   filter(!is.na(LOP868_305_GT)) %>%
#   filter(LOP868_305_DP >= 10) %>%
#   filter(LOP868_305_ID >= 5) %>%
#   filter(LOP868_305_ID / LOP868_305_DP >= 0.1) %>%
#   group_by(LOP868_305_GT) %>%
#   tally
# print(summary305r)
#
# # Code to make supplementary figures from generated datasets
#
# Fig. S7
af.004 <- vcf_all_filt %>%
  filter(LOP868_004_DP >= 10) %>%
  filter(!is.na(LOP868_004_DP)) %>%
  filter(hi_intro_004) %>%
  mutate(hi_intro_call = ifelse(LOP868_004_GT == nonhi_gt, FALSE, TRUE)) %>%
  ggplot(., aes(x = LOP868_004_ID / LOP868_004_DP)) +
  geom_histogram(binwidth = 0.01) +
  labs(y = "Count", x = "") +
  ggtitle("LOP868.004") +
  theme_bw()

af.064 <- vcf_all_filt %>%
  filter(LOP868_064_DP >= 10) %>%
  filter(!is.na(LOP868_064_DP)) %>%
  filter(hi_intro_064) %>%
  mutate(hi_intro_call = ifelse(LOP868_064_GT == nonhi_gt, FALSE, TRUE)) %>%
  ggplot(., aes(x = LOP868_064_ID / LOP868_064_DP)) +
  geom_histogram(binwidth = 0.01) +
  labs(y = "Count", x = "") +
  ggtitle("LOP868.064") +
  scale_y_continuous(limits = c(0,200)) +
  theme_bw()

af.305 <- vcf_all_filt %>%
  filter(LOP868_064_DP >= 10) %>%
  filter(!is.na(LOP868_305_DP)) %>%
  filter(hi_intro_305) %>%
  mutate(hi_intro_call = ifelse(LOP868_305_GT == nonhi_gt, FALSE, TRUE)) %>%
  ggplot(., aes(x = LOP868_305_ID / LOP868_305_DP)) +
  geom_histogram(binwidth = 0.01) +
  labs(y = "Count", x = "HI Allele Depth Fraction") +
  ggtitle("LOP868.305") +
  scale_y_continuous(limits = c(0,200)) +
  theme_bw()

s7 <- ggarrange(af.004, af.064, af.305, nrow = 3)
ggsave("FigS7.png", plot = s7, width = 6.5, height = 4, units = "in", device = "png")
ggsave("FigS7.pdf", plot = s7, width = 6.5, height = 4, units = "in", device = "pdf")

# Fig. S8
p1 <- vcf_all_filt %>%
  filter(!is.na(LOP868_004_GT)) %>%
  mutate(`Locus Type` = ifelse(hi_intro_004, "Introgression", "Not Introgression")) %>%
  # filter(hi_intro_004) %>%
  ggplot(., aes(x = LOP868_DP, fill = `Locus Type`)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0,125)) +
  labs(x="Alca Tarma Depth",y="Density",fill="Locus Type") +
  guides(fill = F) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text.align = 0,
        legend.title.align = 0.5)

p2 <- vcf_all_filt %>%
  filter(!is.na(LOP868_064_GT)) %>%
  mutate(`Locus Type` = ifelse(hi_intro_064, "Introgression", "Not Introgression")) %>%
  # filter(hi_intro_004) %>%
  ggplot(., aes(x = LOP868_DP, fill = `Locus Type`)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0,125)) +
  labs(x="Alca Tarma Depth",y="Density",fill="Locus Type") +
  guides(fill = F) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text.align = 0,
        legend.title.align = 0.5)

p3 <- vcf_all_filt %>%
  filter(!is.na(LOP868_305_GT)) %>%
  mutate(`Locus Type` = ifelse(hi_intro_305, "Introgression", "Not Introgression")) %>%
  # filter(hi_intro_004) %>%
  ggplot(., aes(x = LOP868_DP, fill = `Locus Type`)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0,125)) +
  labs(x="Alca Tarma Depth",y="Density",fill="Locus Type") +
  guides(fill = F) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text.align = 0,
        legend.title.align = 0.5)

p4 <- vcf_all_filt %>%
  filter(!is.na(LOP868_004_GT)) %>%
  mutate(`Locus Type` = ifelse(hi_intro_004, "Introgression", "Not Introgression")) %>%
  # filter(hi_intro_004) %>%
  ggplot(., aes(x = LOP868_004_DP, fill = `Locus Type`)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0,125)) +
  labs(x="Dihaploid Depth",y="Density",fill="Locus Type") +
  guides(fill = F) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text.align = 0,
        legend.title.align = 0.5)

p5 <- vcf_all_filt %>%
  filter(!is.na(LOP868_064_GT)) %>%
  mutate(`Locus Type` = ifelse(hi_intro_064, "Introgression", "Not Introgression")) %>%
  # filter(hi_intro_004) %>%
  ggplot(., aes(x = LOP868_064_DP, fill = `Locus Type`)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0,125)) +
  labs(x="Dihaploid Depth",y="Density",fill="Locus Type") +
  guides(fill = F) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text.align = 0,
        legend.title.align = 0.5)

p6 <- vcf_all_filt %>%
  filter(!is.na(LOP868_305_GT)) %>%
  mutate(`Locus Type` = ifelse(hi_intro_305, "Introgression", "Not Introgression")) %>%
  # filter(hi_intro_004) %>%
  ggplot(., aes(x = LOP868_305_DP, fill = `Locus Type`)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0,125)) +
  labs(x="Dihaploid Depth",y="Density",fill="Locus Type") +
  guides(fill = F) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text.align = 0,
        legend.title.align = 0.5)

s8 <- grid.arrange(p1,p4,p2,p5,p3,p6, ncol = 2)
ggsave("FigS8.png", width = 6.5, height = 4, units = "in", plot = s8, device = "png")
ggsave("FigS8.pdf", width = 6.5, height = 4, units = "in", plot = s8, device = "pdf")

# Fig. S9--need to include code that makes this dataset first!

print("Calculating SNP run lengths to make Fig. S9")

calc_runlengths <- function(df) {
  chex <- which(df$hi_intro_call == T) # the whole thing. Runs pretty fast.
  chex.df <- data.frame("CHROM" = rep(NA, length(chex)),
                        "POS" = rep(NA, length(chex)),
                        "runlen_snp" = rep(NA, length(chex)),
                        "min_runlen_bp" = rep(NA, length(chex)),
                        "max_runlen_bp" = rep(NA, length(chex)),
                        "converted_markers" = rep(NA, length(chex)))
  ct = 0
  for (i in chex) {

    ct = ct + 1

    # recursive look back to determine if also an introgression call.
    # note, this does not account for the inducer also being heterozygous. would have to add this to the code.
    # even with a liberal interpretation that an incoming locus transmitted Alca Tarma alleles across the tract, what are the minimum and maximum lengths of the introgression events?
    decrement = 0
    converted_markers = 1

    while (1) {

      # if chromosomes are not equal, the run is ended. Break out of the while loop.
      if (df$CHROM[i] != df$CHROM[i-decrement-1]) {
        break
      }

      # if chromosomes are equal, start lookin' backwards
      # first, if call behind is also an introgression, then add 1 to converted SNP and keep going backwards
      else if (df$hi_intro_call[i-decrement-1]) {
        converted_markers <- converted_markers + 1
        decrement <- decrement + 1
        next
      }

      # if call behind is heterozygous in either haploid inducer, don't tally a converted marker, but keep going backwards
      else if (df$hi_info[i-decrement-1] == "ambiguous") {
        decrement <- decrement + 1
        next
      }

      # call behind is not an introgression and haploid inducer information is unambiguous. Break out of the while loop.
      else {break}
    }

    # next, look ahead. this is inefficient, as it will run the seed/extend  for each seed, regardless of whether multiple
    # seeds are in the same introgression tract.
    increment = 0
    while (1) {

      # if chromosomes are not equal, the run is ended. break out of the while loop.
      if (df$CHROM[i] != df$CHROM[i+increment+1]) {break}

      else if (df$hi_intro_call[i+increment+1]) {
        converted_markers <- converted_markers + 1
        increment <- increment + 1
        next
      }

      else if (df$hi_info[i+increment+1] == "ambiguous") {
        increment <- increment + 1
        next
      }

      else {break}
    } # end look ahead while loop

    # introgression run length in number of adjacent SNP loci
    runlength_SNP <- increment + decrement + 1

    # tract length bounded by possible introgression markers
    runlength_min_bp <- ifelse(df$POS[i+increment] - df$POS[i-decrement] == 0, 1, df$POS[i+increment] - df$POS[i-decrement])

    # tract length up to next unambiguous introgression-negative marker
    runlength_max_bp <- ifelse(df$POS[i+increment+1] - df$POS[i-decrement-1] == 0, 1, df$POS[i+increment+1] - df$POS[i-decrement-1])

    # record SNP run length and distance in physical bp of putative introgression tracts
    chex.df$CHROM[ct] <- df$CHROM[i]
    chex.df$POS[ct] <- df$POS[i]
    chex.df$runlen_snp[ct] <- runlength_SNP
    chex.df$min_runlen_bp[ct] <- runlength_min_bp
    chex.df$max_runlen_bp[ct] <- runlength_max_bp
    chex.df$converted_markers[ct] <- converted_markers

  } # closes for loop

  return(chex.df)

} # closes calc_runlengths function call

# compute run lengths of introgression loci
RL_004 <- vcf_all_filt %>%
  filter(LOP868_004_DP >= 10) %>%
  filter(!is.na(LOP868_004_GT)) %>%
  mutate(hi_intro_call = ifelse(LOP868_004_GT == nonhi_gt, FALSE, TRUE)) %>%
  mutate(ind_call = ifelse(IVP101_GT == "0/1" | PL4_GT == "0/1", "het", "hom")) %>%
  select(CHROM, POS, hi_intro_call, hi_info, ind_call, starts_with("LOP868_004"), starts_with("LOP868"), starts_with("IVP101"), starts_with("PL4")) %>%
  calc_runlengths(.) %>%
  mutate(dh = "LOP868_004")

RL_064 <- vcf_all_filt %>%
  filter(LOP868_064_DP >= 10) %>%
  filter(!(is.na(LOP868_064_GT))) %>%
  mutate(hi_intro_call = ifelse(LOP868_064_GT == nonhi_gt, FALSE, TRUE)) %>%
  select(CHROM, POS, hi_intro_call, hi_info, starts_with("LOP868"), starts_with("IVP101"), starts_with("PL4")) %>%
  calc_runlengths(.) %>%
  mutate(dh = "LOP868_064")

RL_305 <- vcf_all_filt %>%
  filter(LOP868_305_DP >= 10) %>%
  filter(!(is.na(LOP868_305_GT))) %>%
  mutate(hi_intro_call = ifelse(LOP868_305_GT == nonhi_gt, FALSE, TRUE)) %>%
  select(CHROM, POS, hi_info, hi_intro_call, starts_with("LOP868"), starts_with("IVP101"), starts_with("PL4")) %>%
  calc_runlengths(.) %>%
  mutate(dh = "LOP868_305")

RL_all <- rbind(RL_004, RL_064, RL_305) %>%
  mutate(run_start = NA) %>%
  mutate(run_end = NA)

# collapse runs of adjacent plausible SNP to genomic interval
RL_collapsed <- data.frame("CHROM" = NA,
                           "POS" = NA,
                           "runlen_snp" = NA,
                           "min_runlen_bp" = NA,
                           "max_runlen_bp" = NA,
                           "converted_markers" = NA,
                           "dh" = NA,
                           "run_start" = NA,
                           "run_end" = NA)

rl_all_rowct <- 1 # counts rows traversed in RL_all
rl_collapsed_rowct <- 2 # rows traversed in RL_collapsed

# while(rl_all_rowct < 1000) {
while (rl_all_rowct < nrow(RL_all)) {

  runlen <- RL_all$converted_markers[rl_all_rowct]
  # print(paste("runlen", runlen))

  if (runlen == 1) {
    RL_collapsed <- bind_rows(RL_collapsed, RL_all[rl_all_rowct,])
    # print("rows bound")
    # print(nrow(RL_collapsed))
    # print(RL_all[rowct,])
    # print(runlen)
    RL_collapsed$run_start[rl_collapsed_rowct] <- RL_all$POS[rl_all_rowct] # incorrect indexing here, after adding from longer run lengths, is now off
    RL_collapsed$run_end[rl_collapsed_rowct] <- RL_all$POS[rl_all_rowct] # incorrect indexing here
    rl_all_rowct <- rl_all_rowct + 1
    rl_collapsed_rowct <- rl_collapsed_rowct + 1
  }

  else {
    # print(RL_all[rowct,])
    # print(runlen)
    RL_collapsed <- bind_rows(RL_collapsed, RL_all[rl_all_rowct,])
    # print(RL_collapsed)
    # print(nrow(RL_collapsed))
    # print(paste("first rowct:",rowct))
    # print(paste("runlen: ", runlen))
    RL_collapsed$run_start[rl_collapsed_rowct] <- RL_all$POS[rl_all_rowct]
    # print("multi marker start added")
    # print(RL_all$POS[rowct])
    RL_collapsed$run_end[rl_collapsed_rowct] <- RL_all$POS[rl_all_rowct+runlen-1] # this line seems to be the culprit
    # print("multi marker end added")
    # print(paste("multi marker run end: ", RL_all$POS[rowct+runlen-1]))
    rl_all_rowct <- rl_all_rowct + runlen
    rl_collapsed_rowct <- rl_collapsed_rowct + 1
    # print(RL_all$POS[rowct+runlen])
  }
}

s9a <- ggplot(RL_collapsed, aes(y = converted_markers, x = runlen_snp)) +
  geom_point(alpha = 0.3) +
  labs(y = "Parental Markers\nwith HI alleles", x = "Parental Markers Traversed") +
  

s9b <- ggplot(RL_collapsed, aes(y = converted_markers, x = runlen_snp)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  labs(y = "Parental Markers\nwith HI alleles", x = "log10 Parental Markers Traversed") +
  # scale_y_continuous(limits = c(0,2000))
  scale_x_log10()

s9 <- grid.arrange(s9a,s9b, nrow = 2)
ggsave("FigS9.png", width = 6.5, height = 6.5, units = "in", device = "png", plot = s9)
ggsave("FigS9.pdf", width = 6.5, height = 6.5, units = "in", device = "pdf", plot = s9)

# Figures for Fig. 4A and S3 chromosome panels
# Note, since I can't make this entire figure until I have the population coverage data, I may want to move appropriate tables
# to a results folder, and then call an Rscript from there to make figures. 

# columns of summarized depth files
dp_cols <- c("chrom", "start", "end", "mean_dp", "median_dp")

dp_files <- dir(pattern = "summarized-depth-4x_LOP868", path = "data/merged", full.names = T)
bins_10k <- dir(pattern = "10k_windows", path = "../ref/", full.names = T) %>% 
  read_tsv(., col_names = c("chrom", "start", "end")) %>% 
  merge(read_tsv(dp_files, col_names = dp_cols), .) %>% 
  mutate(ncon = replace(median_dp, N_content >= 0.1, NA)) %>% 
  mutate(gc_norm_rd = median_dp / GC_content)

# mean shift cluster random subset of median 10k window read depth values to find latent CN states
set.seed(5)
mtx <- t(matrix(data = c(sample(na.omit(bins_10k$gc_norm), 1000), rep(0,1000)), ncol=2))
clusters <- msClustering(mtx, kernel = "epanechnikovKernel", h= quantile(dist(mtx[1,]), 0.25))
sc <- clusters$components[1,1]

# Fig. 4A, chr06
fig4a <- bins_10k %>% 
  filter(chrom == "chr06" & N_content <= 0.3) %>% 
  ggplot(., aes(x = start, y = gc_norm_rd)) + 
  geom_point(alpha = 0.3, size = 0.8) +
  labs(x = "", y = "Alca Tarma Copy Number") +
  ggtitle("A") +
  scale_y_continuous(limits = c(0, 220), breaks = c(0, scL*.25, scL*.5, scL*0.75, scL, scL*1.25, scL*1.5, scL*1.75, scL*2), labels = c(0,1,2,3,4,5,6,7,8)) +
  scale_x_continuous(limits = c(0,6e7), breaks = seq(0, 6e7, by = 2e7), labels = seq(0, 60, by = 20)) +
  facet_wrap(~chrom, nrow = 12, strip.position = "right") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray30", size = 0.2),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    # plot.title = element_text(size = 24, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")
  )

# Fig. 4B, zoomed in view of chr06
fig4b <- bins_10k %>% 
  filter(chrom == "chr06" & N_content <= 0.3 & start >= 3.2e7 & start <= 3.4e7) %>% 
  ggplot(., aes(x = start, y = gc_norm_rd)) + 
  geom_point() +
  geom_line() +
  labs(x = "", y = "Alca Tarma Copy Number") +
  ggtitle("B") +
  scale_x_continuous(breaks = seq(3.2e7, 3.4e7, by = 2.5e5), labels = seq(32,34,by=0.25)) +
  scale_y_continuous(limits = c(0, 220), breaks = c(0, scL*.25, scL*.5, scL*0.75, scL, scL*1.25, scL*1.5, scL*1.75, scL*2), labels = c(0,1,2,3,4,5,6,7,8)) +
  facet_wrap(~chrom, nrow = 12, strip.position = "right") +
  theme(
    # panel.grid.major.x = element_line(color = "gray30", size = 0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray30", size = 0.2),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    # plot.title = element_text(size = 24, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")
  )
