library(tidyverse)
library(here)
library(fuzzyjoin)

## The following code generates summary statistics from the ABBA-BABA analysis.

##    FILES REQUIRED:
##          Results from ABBA-BABA analysis (emoryi_slowinskii_guttatus_localFstats__50_25.txt)
##          Genomic locations for NMT and control genes (Nmt_coords.txt and controls_coords.txt)

##    STRUCTURE OF CODE:
##              (1) Import data
##              (2) Get summary statistics


# (1) Import data ---------------------------------------------------------

# ABBA-BABA data
dat <- read_tsv(here("data", "emoryi_slowinskii_guttatus_localFstats__50_25.txt"), col_names = TRUE) # 120,023

# Import Nmt information
nmts <- read_tsv(here("data", "Nmt_coords.txt"), col_names = FALSE) %>% 
  separate(col = X1, sep = ":", into = c("chr", "sites")) %>% 
  separate(col = sites, sep = "-", into = c("start", "end")) %>% 
  rename(gene = X2)
nmts$start <- as.numeric(nmts$start)
nmts$end <- as.numeric(nmts$end) # 167 genes

# Import control gene information
cont <- read_tsv(here("data", "controls_coords.txt"), col_names = FALSE) %>% 
  separate(col = X1, sep = ":", into = c("chr", "sites")) %>% 
  separate(col = sites, sep = "-", into = c("start", "end")) %>% 
  rename(gene = X2)
cont$start <- as.numeric(cont$start)
cont$end <- as.numeric(cont$end) # 142 genes


# (2) Get summary statistics ----------------------------------------------

# Get total numbers of SNPs for cont and nmt datasets
cont %>% 
  mutate(length = end-start) %>% 
  summarize(sum(length)) # 10049438

nmts %>% 
  mutate(length = end-start) %>% 
  summarize(sum(length)) # 2348450

# Now, let's get a genome-wide mean estimate for fdM for nmts and cont genes
nmtchr <- unique(nmts$chr) # 74 unique chroms
nmtdat <- dat %>% filter(chr %in% nmtchr) # reduce dataset size
nmtjoin <-
  fuzzyjoin::fuzzy_inner_join(nmts, nmtdat, 
                              by = c("chr" = "chr",
                                     "start" = "windowStart",
                                     "end" = "windowEnd"),
                              match_fun = list(`==`, `>=`, `<=`))
write_tsv(nmtjoin, here("ABBABABA", "Joined_nmts.txt"), col_names = TRUE)
nmtjoin %>% summarize(meanfdM = mean(f_dM),
                      minfdM = min(f_dM),
                      maxfdM = max(f_dM)) # 0.0546 mean

contchr <- unique(cont$chr) # 87 unique chroms
contdat <- dat %>% filter(chr %in% contchr) # reduce dataset size
contjoin <-
  fuzzyjoin::fuzzy_inner_join(cont, contdat, 
                              by = c("chr" = "chr",
                                     "start" = "windowStart",
                                     "end" = "windowEnd"),
                              match_fun = list(`==`, `>=`, `<=`))
write_tsv(contjoin, here("ABBABABA", "Joined_cont.txt"), col_names = TRUE)
contjoin %>% summarize(meanfdM = mean(f_dM),
                       minfdM = min(f_dM),
                       maxfdM = max(f_dM)) # 0.0607 mean

# TODO run statistical test!

