library(here)
library(tidyverse)
library(fuzzyjoin)
library(padr)

## The following code processes results from the GWAS analysis and merges these with the NMT and
## control genes.

##    FILES REQUIRED:
##          Results from the GWAS analysis (25 RData files)
##          Coordinates for the NMT and control genes of interest (NMT_coords.txt & controls_coords.txt)
##          N.B.: Relevant functions are in the GWAS.R script

##    STRUCTURE OF CODE:
##              (1) Process NMT and control data
##              (2) Process GWAS data
##              (3) Merge GWAS results with NMT and control data, exporting sig SNPs
##              (4) Get summary statistics

source(here("analysis", "GWAS.R"))
setwd("~/Box Sync/Lampropeltis project/PANTHEROPHIS_introgression")


# (1) Process NMT and control data ----------------------------------------

# NMT and control files only have ranges so we need to pad data such that
# each position is assigned a row in the dataframe
nmts <- read_tsv("Gene_trees/Nmt_coords.txt", col_names = FALSE) %>% 
  rename(seqname = X1, gene = X2) %>% 
  separate(col = seqname, into = c("chrom", "number"), sep = "\\:") %>% # may want to check that col has unique entries for each row `nrow(nmts) == length(unique(nmts$number))` before moving to next line
  separate(col = number, into = c("start", "end"), sep = "-")
nmts$start <- as.numeric(nmts$start)
nmts$end <- as.numeric(nmts$end)

nmt_gathered <-
  nmts %>% 
  # add unique values for each row
  mutate(gene_group = row_number()) %>% 
  # gather start and end cols
  gather(key = "pos_name", value = "pos",
         start:end) %>% 
  # arrange ascending by row val
  arrange(gene_group, pos) %>% 
  # pad range values
  pad_int(by = "pos", group = "gene_group", step = 1) %>% 
  group_by(gene_group) %>% 
  # fill NA cells with duplicate values; default direction is down
  # need to specify which cols to fill; everything is all cols
  fill(everything()) %>% 
  dplyr::select(-pos_name) %>% 
  mutate(category = "NMT")

cont <- read_tsv("Gene_trees/controls_coords.txt", col_names = FALSE) %>% 
  rename(seqname = X1, gene = X2) %>% 
  separate(col = seqname, into = c("chrom", "number"), sep = "\\:") %>% # may want to check that col has unique entries for each row `nrow(nmts) == length(unique(nmts$number))` before moving to next line
  separate(col = number, into = c("start", "end"), sep = "-")
cont$start <- as.numeric(cont$start)
cont$end <- as.numeric(cont$end)

cont_gathered <-
  cont %>% 
  # add unique values for each row
  mutate(gene_group = row_number()) %>% 
  # gather start and end cols
  gather(key = "pos_name", value = "pos",
         start:end) %>% 
  # arrange ascending by row val
  arrange(gene_group, pos) %>% 
  # pad range values
  pad_int(by = "pos", group = "gene_group", step = 1) %>% 
  group_by(gene_group) %>% 
  # fill NA cells with duplicate values; default direction is down
  # need to specify which cols to fill; everything is all cols
  fill(everything()) %>% 
  dplyr::select(-pos_name) %>% 
  mutate(category = "cont")


# (2) Process GWAS data ---------------------------------------------------

# Iterate through all the RData files, read them in, calculate which SNPs are 
# significant, and export these data as a single file
files <- list.files("GWAS/", pattern = "RData")
file_no <- 1:length(files)

dat <-
  file_no %>% 
  lapply(function(x) {
    load(here("GWAS", files[x]))
    gwas <- gwas_stats(gwas)
    return(gwas)
  }) %>% 
  dplyr::bind_rows()

# Save GWAS results for all SNPs
write_tsv(dat, here("GWAS", "GWAS_results.txt"), col_names = TRUE)

# Retrieve only significant outliers based on two alpha thresholds
# sig_snps_0.01 <- dat %>%
#   dplyr::filter(signed.logp > -log10(0.01)) # 144,823 at alpha=0.01
# sig_snps_0.01$chrom <- as.character(sig_snps_0.01$chrom)

sig_snps_0.05 <- dat %>% 
  dplyr::filter(signed.logp > -log10(0.05)) # 954,013 at alpha=0.05
sig_snps_0.05$chrom <- as.character(sig_snps_0.05$chrom)

# Get summary statistics (number of outliers per chrom)
# chrom_sites_0.01 <- sig_snps_0.01 %>% group_by(chrom) %>% count() %>% ungroup()
# chrom_sites_0.01 %>% summarize(min_outliers = min(n),
#                                              max_outliers = max(n),
#                                              mean_outliers = mean(n))

chrom_sites_0.05 <- sig_snps_0.05 %>% group_by(chrom) %>% count() %>% ungroup()
chrom_sites_0.05 %>% summarize(min_outliers = min(n),
                                             max_outliers = max(n),
                                             mean_outliers = mean(n))


# (3) Join GWAS results with NMTs and cont genes --------------------------

# Do below for both alpha thresholds
nmtdat <- left_join(sig_snps_0.05, nmt_gathered) %>% 
  dplyr::select(chrom, pos, zscore, pos.Mb, signed.logp, r2, category) %>%
  distinct() %>%
  replace_na(list(category = "non-NMT"))
# write_tsv(nmtdat, here("GWAS", "GWAS_NMT_sigsnps_0.05.txt")) # change file name accordingly

nmtdat %>% 
  filter(category == "NMT") %>% 
  summarize(n()) # 128 SNPs at 0.01 / 1169 SNPs at 0.05

### Same as above but for control genes
contdat <- left_join(sig_snps_0.05, cont_gathered) %>% 
  dplyr::select(chrom, pos, zscore, pos.Mb, signed.logp, r2, category) %>%
  distinct() %>%
  replace_na(list(category = "non-control"))
# write_tsv(contdat, here("GWAS", "GWAS_control_sigsnps_0.05.txt"))

contdat %>% 
  filter(category == "cont") %>% 
  summarize(n()) # 529 SNPs at 0.01 / 5592 SNPs at 0.05


# (4) Get summary statistics ----------------------------------------------

# Which chroms have both Nmt and control genes?
chroms <- inner_join(cont %>% dplyr::select(chrom) %>% distinct(), 
                     nmts %>% dplyr::select(chrom) %>% distinct())
