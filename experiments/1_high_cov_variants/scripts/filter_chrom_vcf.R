#!/share/comailab/kramundson/miniconda3/envs/ximena/bin/Rscript
# Kirk Amundson
# filter_chrom_vcf.R

# USAGE: filter_chrom_vcf.R data/calls/{chrom}-calls.vcf
# Is designed to be run in Snakemake workflow

library(tidyverse)
library(egg)
library(stringr)

args <- commandArgs(trailingOnly = T)

filter_vcf <- function(file) {
  print(paste("Now filtering", file))
  
  # parse header from VCF
  vcf_header <- system(paste("grep '#C'", file), intern = T) %>% 
    str_replace("^#C", "C") %>% 
    str_replace_all("[0-9]x_", "") %>% 
    str_split(pattern = "\t")
  
  # read VCF into R
  vcf <- read_tsv(file, col_names = vcf_header[[1]],
                  comment = "#", na = c("NA", ".", "./.", "./././."))
  
  print("VCF was successfully read into R, proceeding with parsing")
  
  # figure out how many different types of INFO fields there are. The parsing approach used here will only work if they're the same
  n_info <- vcf %>% 
    mutate(INFO = str_replace_all(INFO, "=\\.+", "")) %>% 
    select(INFO) %>% 
    unique() %>% 
    length()
  
  print(paste("R determined that the number of unique INFO columns was", n_info))
  
  # here, should add a check to make sure that splitting INFO using only the first row is safe
  # print(table(vcf_inf$INFO))
  stopifnot(n_info == 1)
  
  # parse INFO column, contains information about sites, aggregated across all samples
  info <- str_split(vcf$INFO[1], ";")[[1]] %>% 
    str_replace("=.+", "")
  
  vcf_inf <- vcf %>% 
    mutate(INFO = str_replace_all(INFO, "[A-Za-z]*=", "")) %>% 
    separate(INFO, into = info, sep = ";", convert = T)
  
  
  print("Implementing Site Filters")
  # implement site filters
  vcf_filter <- vcf_inf %>% 
    filter(NUMALT == 1) %>% 
    filter(CIGAR == "1X") %>% 
    mutate(MQM = as.numeric(MQM)) %>% 
    filter(MQM >= 50) %>% 
    filter(MQMR >= 50) %>% 
    mutate(MQ.diff = abs(MQM-MQMR)) %>% 
    filter(MQ.diff < 10) %>% 
    mutate(SAF = as.numeric(SAF)) %>% 
    mutate(AO = as.numeric(AO)) %>% 
    mutate(RPPR = as.numeric(RPPR)) %>% 
    filter(RPPR <= 20) %>% 
    mutate(RPP = as.numeric(RPP)) %>% 
    filter(RPP <= 20) %>% 
    mutate(EPP = as.numeric(EPP)) %>% 
    filter(EPP <= 20) %>% 
    mutate(EPPR = as.numeric(EPPR)) %>% 
    filter(EPPR <= 20) %>% 
    mutate(SAP = as.numeric(SAP)) %>% 
    filter(SAP <= 20) %>% 
    mutate(SRP = as.numeric(SRP)) %>% 
    filter(SRP <= 20)
  
  print(paste("Number of loci retained after site filters:", nrow(vcf_filter)))
  
  # parse FORMAT column; this field contains sample-specific attributes
  print("parsing FORMAT")
  attributes <- str_split(names(table(vcf_inf$FORMAT)), ":", simplify = F)
  print(attributes)
  
  vcf_format_split <- split(vcf_filter, f=vcf_filter$FORMAT)
  
  # Note: in past runs, only the 2nd class of FORMAT is kept after implementing sample-specific
  # filters, as LOP868-specific values are all converted to NA. For function, I'll split all
  # entries in the FORMAT-split list and print out table(LOP868_GT) to validate that all entries
  # in the table are missing (NA). If they are not, you need to go back and filter the [[1]] calls.
  vcf_rejoin <- vcf_format_split[[2]] %>% 
    separate(LOP868, sep = ":", convert = T,
             into = paste("LOP868", attributes[[2]], sep = "_")) %>% 
    separate(PL4, sep = ":", convert = T,
             into = paste("PL4", attributes[[2]], sep = "_")) %>% 
    separate(IVP101, sep = ":", convert = T,
             into = paste("IVP101", attributes[[2]], sep = "_")) %>% 
    separate(LOP868_004, sep = ":", convert = T,
             into = paste("LOP868_004", attributes[[2]], sep = "_")) %>% 
    separate(LOP868_064, sep = ":", convert = T,
             into = paste("LOP868_064", attributes[[2]], sep = "_")) %>% 
    separate(LOP868_305, sep = ":", convert = T,
             into = paste("LOP868_305", attributes[[2]], sep = "_")) %>%
    # filter(LOP868_DP >= 20 & LOP868_DP <= 70) %>% # probably didn't have the coverage to have SNP in 
    # filter(IVP101_DP >= 20 & IVP101_DP <= 60) %>%
    # filter(PL4_DP >= 20 & PL4_DP <= 80) %>%
    # filter(LOP868_GQ >= 35) %>%
    # filter(IVP101_GQ >= 35) %>%
    # filter(PL4_GQ >= 35) %>%
    filter(LOP868_GT %in% c("0/0/0/0", "1/1/1/1")) %>% 
    filter(!(LOP868_GT == "0/0/0/0" & IVP101_GT == "0/0")) %>%
    filter(!(LOP868_GT == "1/1/1/1" & IVP101_GT == "1/1")) %>%
    filter(!(LOP868_GT == "0/0/0/0" & PL4_GT == "0/0")) %>%
    filter(!(LOP868_GT == "1/1/1/1" & PL4_GT == "1/1")) %>% 
    mutate(LOP868_ID = ifelse(LOP868_GT == "0/0/0/0", LOP868_AO, LOP868_RO)) %>% 
    mutate(LOP868_004_ID = ifelse(LOP868_GT == "0/0/0/0", LOP868_004_AO, LOP868_004_RO)) %>% 
    mutate(LOP868_064_ID = ifelse(LOP868_GT == "0/0/0/0", LOP868_064_AO, LOP868_064_RO)) %>% 
    mutate(LOP868_305_ID = ifelse(LOP868_GT == "0/0/0/0", LOP868_305_AO, LOP868_305_RO)) %>% 
    filter(LOP868_ID == 0) %>% 
    filter(IVP101_DP >= 10) %>%
    filter(PL4_DP >= 10) %>%
    filter(LOP868_DP >= 10) %>%
    filter(DP <= 1.6 * 215) %>% # new, makes current SNP list obsolete
    mutate(hi_info = ifelse(IVP101_GT == "0/1" | PL4_GT == "0/1", "ambiguous", "unambiguous")) %>%
    mutate(perHiReads = LOP868_004_ID / LOP868_004_DP) %>%
    mutate(nonhi_gt = str_replace(LOP868_GT, "[0-9]/[0-9]/", "")) %>%
    mutate(hi_intro_004 = ifelse(LOP868_004_GT == nonhi_gt, FALSE, TRUE)) %>%
    mutate(hi_intro_064 = ifelse(LOP868_064_GT == nonhi_gt, FALSE, TRUE)) %>% 
    mutate(hi_intro_305 = ifelse(LOP868_305_GT == nonhi_gt, FALSE, TRUE)) %>% 
    mutate(IVP101_AB = IVP101_AO / IVP101_DP) %>%
    mutate(PL4_AB = PL4_AO / PL4_DP) %>%
    filter(between(IVP101_AB, 0.4, 0.6) | between(IVP101_AB, 0.95, 1) | between(IVP101_AB, 0, 0.05)) %>%
    filter(between(PL4_AB, 0.4, 0.6) | between(PL4_AB, 0.95, 1) | between(PL4_AB, 0, 0.05))
  
  # write filtered calls to disk as a tab-delimited gzip file
  writeout <- paste0(gsub("calls.vcf", "filtered-calls.vcf", file), '.gz')
  print(paste("Writing out filtered calls to:", writeout))
  write.table(arrange(vcf_rejoin, POS), gzfile(writeout),
              quote=F, sep = '\t', eol = '\n',
              na = "NA", row.names = F, col.names = T)
}

filter_vcf(args[1])
