library(tidyverse)
library(here)
library(fuzzyjoin)

## The following code generates summary statistics from the Dinvestigate (fstat) analysis.

##    FILES REQUIRED:
##          Results from Dinvestigate analysis (emoryi_slowinskii_guttatus_localFstats__50_25.txt)
##          Genomic locations for NMT and control genes (Nmt_coords.txt and control_coords2.txt)

##    STRUCTURE OF CODE:
##              (1) Import data
##              (2) Get summary statistics
##              (3) Bootstrap test


# (1) Import data ---------------------------------------------------------

# fstat data
dat <- read_tsv(here("data", "emoryi_slowinskii_guttatus_localFstats__50_25.txt"), col_names = TRUE) # 120,023

# Import Nmt information
nmts <- read_tsv(here("data", "Nmt_coords.txt"), col_names = FALSE) %>% 
  separate(col = X1, sep = ":", into = c("chr", "sites")) %>% 
  separate(col = sites, sep = "-", into = c("start", "end")) %>% 
  rename(gene = X2)
nmts$start <- as.numeric(nmts$start)
nmts$end <- as.numeric(nmts$end) # 167 genes

# Retrieve gene names from old control coordinates file
old_cont <- read_tsv(here("data", "controls_coords.txt"), col_names = FALSE) %>% 
  separate(col = X1, sep = ":", into = c("chr", "sites")) %>% 
  separate(col = sites, sep = "-", into = c("start", "end")) %>% 
  rename(gene = X2)
old_cont$start <- as.numeric(old_cont$start)
old_cont$end <- as.numeric(old_cont$end)

# Import control gene information
cont <- read_tsv(here("data", "control_coords2.txt"), col_names = FALSE) %>% 
  rename("chr" = X1, "start" = X2, "end" = X3)
cont$start <- as.numeric(cont$start)
cont$end <- as.numeric(cont$end) # 139 genes

cont <- left_join(cont, old_cont)


# (2) Get summary statistics ----------------------------------------------

# Get total numbers of SNPs for cont and nmt datasets
cont %>% 
  mutate(length = end-start) %>% 
  summarize(sum(length)) # 9819152

nmts %>% 
  mutate(length = end-start) %>% 
  summarize(sum(length)) # 2348450

# Now, let's get a genome-wide mean estimate for fdM for nmts and cont genes
nmtchr <- unique(nmts$chr) # 87 unique chroms
nmtdat <- dat %>% filter(chr %in% nmtchr) # reduce dataset size prior to joining
nmtjoin <-
  fuzzyjoin::fuzzy_inner_join(nmts, nmtdat, 
                              by = c("chr" = "chr",
                                     "start" = "windowStart",
                                     "end" = "windowEnd"),
                              match_fun = list(`==`, `>=`, `<=`))
write_tsv(nmtjoin, here("data", "Joined_nmts.txt"), col_names = TRUE) # 191 obs
nmtjoin %>% summarize(meanfdM = mean(f_dM),
                      minfdM = min(f_dM),
                      maxfdM = max(f_dM)) # 0.0546 mean

contchr <- unique(cont$chr) # 73 unique chroms
contdat <- dat %>% filter(chr %in% contchr) # reduce dataset size
contjoin <-
  fuzzyjoin::fuzzy_inner_join(cont, contdat, 
                              by = c("chr" = "chr",
                                     "start" = "windowStart",
                                     "end" = "windowEnd"),
                              match_fun = list(`==`, `>=`, `<=`))
write_tsv(contjoin, here("data", "Joined_cont.txt"), col_names = TRUE)
contjoin %>% summarize(meanfdM = mean(f_dM),
                       minfdM = min(f_dM),
                       maxfdM = max(f_dM)) # 0.0607 mean


# (3) Bootstrap for confidence intervals ----------------------------------

# Calculate mean and then subtract means to get bootstrap sample
bootstrap_rep <- function(boot_dat, prop_samp) {
  # Take random samples ~10% of dataframe sizes
  nmt <-
    boot_dat %>% 
    filter(gene_set == "nmt") %>% 
    sample_n(size = round(prop_samp*n()), replace = TRUE) %>% 
    summarize(mean_fdm = mean(f_dM)) %>% 
    pull(mean_fdm)

  cont <-
    boot_dat %>% 
    filter(gene_set == "control") %>% 
    sample_n(size = round(prop_samp*n()), replace = TRUE) %>% 
    summarize(mean_fdm = mean(f_dM)) %>% 
    pull(mean_fdm)
  
  # Calculate difference between means
  bs_samp <- nmt-cont
  
  return(bs_samp)
}

nmtjoin <- read_tsv(here("data", "Joined_nmts.txt"), col_names = TRUE) %>% 
  mutate(gene_set = "nmt") # 191
contjoin <- read_tsv(here("data", "Joined_cont.txt"), col_names = TRUE) %>% 
  mutate(gene_set = "control") # 67
boot_dat <- bind_rows(nmtjoin, contjoin)

n_reps = 10000

vec <-
  replicate(n_reps, bootstrap_rep(boot_dat, prop_samp = 0.3))
quantile(vec, probs = c(0.025, 0.975)) # -0.05308854, 0.04219925
mean(vec) # -0.006219514
