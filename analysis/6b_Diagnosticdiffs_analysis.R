library(tidyverse)
library(seqinr)
library(vcfR)
library(here)
library(adegenet)
library(algatr)

## The following code calculates the admixture index (fixed differences) between individuals
## at the emoryi-slowinskii contact zone based on their mitochondrial haplotype.
## This code is modified from Chambers et al. (2023). doi: https://doi.org/10.1093/sysbio/syac056

##    FILES REQUIRED:
##          vcfs for NMTs and control loci (cznmtsnps.vcf & czcontsnps2.vcf)
##          fastas for NMTs and control loci (cznmtsnps.min4.fasta & czcontsnps2.min4.fasta) <- vcfs converted to fasta using vcf2phylip.py
##          mitotype assignments (cz_mitotypes.txt)
##          vcfs for LD-pruned NMTs and control loci (cznmtsnps_ldp.vcf & czcontsnps2_ldp.vcf) <- generated using `ld_pruning.sh` script
##          fastas for LD-pruned NMTs and control loci (cznmtsnps_ldp.fasta & czcontsnps2_ldp.fasta) <- vcfs converted to fasta using vcf2phylip.py

##    STRUCTURE OF CODE:
##              (1) Read in input files
##              (2) Calculate diagnostic diffs between reference groups
##              (3) Calculate per locus allele frequencies
##              (4) Individual-based statistical test
##              (5) Run a DAPC on un-admixed vs admixed individuals

# Load relevant functions
source(here("analysis", "Diagnosticdiffs.R"))


# (1) Read in input files -------------------------------------------------

### Import mitotype assignments for 30 individuals
mitotypes <- read_tsv(here("data", "cz_mitotypes.txt"), col_names = c("INDV", "mitotype"))

### Run above function for NMTs and control loci
nmts <- process_seq(vcf_file = here("data", "cznmtsnps.vcf"),
                    fasta_file = here("data", "cznmtsnps.min4.fasta"),
                    mitotypes = mitotypes) # 30 inds of 38,553 vars (38,551 SNPs)
cont <- process_seq(vcf_file = here("data", "czcontsnps2.vcf.gz"),
                    fasta_file = here("data", "czcontsnps2.min4.fasta"),
                    mitotypes = mitotypes) # 30 inds of 163,960 vars (163,958 SNPs)


# (2) Calculate diagnostic diffs between reference inds -------------------

# Determine fixed differences between reference groups
dataset = cont # cont or nmts
dataset_name = "cont" # "cont" or "nmts"
# If `save_file` set to TRUE, will save as e.g. "pure_allele_dict_cont.rda"
pure_allele_dict_noambig_cont <- allele_dict(dataset, dataset_name, save_file = TRUE, output_path = here("data", "/")) # only need to run once; 26,928 for control

dataset = nmts # cont or nmts
dataset_name = "nmts" # "cont" or "nmts"
# If `save_file` set to TRUE, will save as e.g. "pure_allele_dict_nmts.rda"
pure_allele_dict_noambig_nmts <- allele_dict(dataset, dataset_name, save_file = FALSE) # only need to run once; 5,954 for NMTs


# (3) Calculate per-locus allele frequencies ------------------------------

# TODO this is no longer necessary because of LD-pruning?
# If above has already been run, can start below; will load object as `pure_allele_dict_noambig`:
dataset = nmts # cont or nmts
dataset_name = "nmts" # "cont" or "nmts"
load(paste0(here("data"), "/pure_allele_dict_", dataset_name, ".rda"))
pure_allele_dict_noambig_nmts <- pure_allele_dict_noambig

dataset = cont # cont or nmts
dataset_name = "cont" # "cont" or "nmts"
load(paste0(here("data"), "/pure_allele_dict_", dataset_name, ".rda"))
pure_allele_dict_noambig_cont <- pure_allele_dict_noambig

# Set a threshold value for diagnostic diffs for each sample (relative to ref inds)
threshold = 0.5 # change to 0.5, 0.75, 0.95

# If `save_file` set to TRUE, will save file as e.g. "summary_nmts_0.5.txt" in `output_path`
# Run below for both nmts and cont gene sets
results <- freq_data(dataset = cont_ldp, 
                     dataset_name = "cont", 
                     pure_allele_dict = pure_allele_dict_noambig_cont, 
                     threshold, 
                     save_file = FALSE, 
                     output_path = paste0(here("data"), "/"))

# Look at results
results$summary


# (4) Logistic regression -------------------------------------------------

# Import NMTs and control SNPs that have been LD-pruned; see `ld_pruning.sh` for code to LD prune data
nmts_ldp <- process_seq(vcf_file = here("data", "cznmtsnps_ldp.vcf"),
                        fasta_file = here("data", "cznmtsnps_ldp.fasta"),
                        mitotypes = mitotypes) # 30 inds, 14,330 SNPs (+INDV +mitotype cols)
cont_ldp <- process_seq(vcf_file = here("data", "czcontsnps2_ldp.vcf"),
                        fasta_file = here("data", "czcontsnps2_ldp.fasta"),
                        mitotypes = mitotypes) # 30 inds, 59,208 SNPs (+INDV +mitotype cols)

freqs_cont <-
  cont_ldp %>% 
  tidyr::pivot_longer(names_to = "locus", values_to = "value", cols = -c(INDV, mitotype)) %>%
  dplyr::filter(value != 'N') %>% 
  dplyr::inner_join(pure_allele_dict_noambig_cont, by = "locus") %>% 
  dplyr::mutate(gene_set = "cont",
                match_type = case_when(value == emoryi ~ "emoryi",
                                       value == slowinskii ~ "slowinskii",
                                       value == het ~ "het")) %>% 
  dplyr::filter(match_type != "het") %>% 
  dplyr::mutate(classification = case_when(match_type == mitotype ~ 1,
                                           match_type != mitotype ~ 0)) # 6,993 loci considered

freqs_nmts <-
  nmts_ldp %>% 
  tidyr::pivot_longer(names_to = "locus", values_to = "value", cols = -c(INDV, mitotype)) %>%
  dplyr::filter(value != 'N') %>% 
  dplyr::inner_join(pure_allele_dict_noambig_nmts, by = "locus") %>% 
  dplyr::mutate(gene_set = "nmts",
                match_type = case_when(value == emoryi ~ "emoryi",
                                       value == slowinskii ~ "slowinskii",
                                       value == het ~ "het")) %>% 
  dplyr::filter(match_type != "het") %>% 
  dplyr::mutate(classification = case_when(match_type == mitotype ~ 1,
                                           match_type != mitotype ~ 0)) # 1,414 loci considered

freqs <- bind_rows(freqs_cont, freqs_nmts)

mod <- glm(classification ~ INDV + gene_set, 
           family = binomial(link = "logit"),
           data = freqs)

summary(mod) # gene_setnmts are statistically significant (p-value = 7.30e-09)

# Odds ratio
exp(coef(mod)) # 1.085002 for gene_setnmts; i.e., 8.5% increase of matches in N-mts compared to control genes

# Save freq data for data visualization
results <- freq_data(dataset = cont_ldp, 
                     dataset_name = "cont_ldp", 
                     pure_allele_dict = pure_allele_dict_noambig_cont, 
                     threshold, 
                     save_file = TRUE, 
                     output_path = paste0(here("data"), "/"))
freqs_to_plot_cont <- results$freqs
save(freqs_to_plot_cont, file = paste0(here("data"), "/freq_data_cont_ldp_", threshold, ".rda"))

results <- freq_data(dataset = nmts_ldp, 
                     dataset_name = "nmts_ldp", 
                     pure_allele_dict = pure_allele_dict_noambig_nmts, 
                     threshold, 
                     save_file = TRUE, 
                     output_path = paste0(here("data"), "/"))
freqs_to_plot_nmts <- results$freqs
save(freqs_to_plot_nmts, file = paste0(here("data"), "/freq_data_nmts_ldp_", threshold, ".rda"))


# (5) Prepare data for DAPC -----------------------------------------------

### Separate datasets into least admixed and admixed individuals
# Import pop structure results from Marshall et al. (2021)
qmat <- read_tsv(here("data", "sample_qvalues.txt"))
# There are 13 slow inds and 17 emoryi inds
pure_slow <-
  qmat %>% 
  filter(mtclade == "SLOW") %>% 
  slice_max(slowinskii, n = 4)
admix_slow <-
  qmat %>% 
  filter(mtclade == "SLOW") %>% 
  slice_min(slowinskii, n = 9)
pure_emoryi <-
  qmat %>% 
  filter(mtclade == "EMNO") %>% 
  slice_max(emoryi, n = 4)
admix_emoryi <-
  qmat %>% 
  filter(mtclade == "EMNO") %>% 
  slice_min(emoryi, n = 13)

# Un-admixed samples
unadmixed <- bind_rows(pure_slow, pure_emoryi)
admixed <- bind_rows(admix_slow, admix_emoryi)

### ------------ N-mt genes ------------ 
# Import vcfs for N-mt and control gene SNPs that have been LD-pruned
nmts_vcf <- read.vcfR(here("data", "cznmtsnps_ldp.vcf"))
names_nmts <- names_helper(nmts_vcf)

# Convert to dosage matrix
nmts_dos <- vcf_to_dosage(nmts_vcf)
rownames(nmts_dos) <- names_nmts$INDV
# Remove any columns that have NAs
nmts_dos_nonas <- nmts_dos[, !apply(is.na(nmts_dos), 2, any)] # 3041 remain

# Subset into unadmixed and admixed inds
nmts_pure <-
  as.data.frame(nmts_dos_nonas) %>% 
  rownames_to_column(var = "INDV") %>% 
  filter(INDV %in% unadmixed$sampleID) %>% 
  column_to_rownames(var = "INDV") %>% 
  # mutate(mitotype = unadmixed_mitotypes$mitotype) %>%
  as.matrix()
# dplyr::select(where(~n_distinct(.) > 1)) # remove monomorphic loci

nmts_admix_slow <-
  as.data.frame(nmts_dos_nonas) %>% 
  rownames_to_column(var = "INDV") %>% 
  filter(INDV %in% admix_slow$sampleID) %>% 
  column_to_rownames(var = "INDV") %>% 
  as.matrix()

nmts_admix_emoryi <-
  as.data.frame(nmts_dos_nonas) %>% 
  rownames_to_column(var = "INDV") %>% 
  filter(INDV %in% admix_emoryi$sampleID) %>% 
  column_to_rownames(var = "INDV") %>% 
  as.matrix()

# Get pop assignment for pure/unadmixed individuals
# We're using the mitotypes object because the ordering is consistent
# with the vcf and dosage matrix
unadmixed_mitotypes <-
  mitotypes %>% 
  filter(INDV %in% unadmixed$sampleID)
# Check ordering
all(rownames(nmts_pure) == unadmixed_mitotypes$INDV)

### ------------ Control genes ------------ 
# Import vcfs for N-mt and control gene SNPs that have been LD-pruned
cont_vcf <- read.vcfR(here("data", "czcontsnps2_ldp.vcf"))
names_cont <- names_helper(cont_vcf)

# Convert to dosage matrix
cont_dos <- vcf_to_dosage(cont_vcf)
rownames(cont_dos) <- names_cont$INDV
# Remove any columns that have NAs
cont_dos_nonas <- cont_dos[, !apply(is.na(cont_dos), 2, any)] # 10368 remain

# Subset into unadmixed and admixed inds
cont_pure <-
  as.data.frame(cont_dos_nonas) %>% 
  rownames_to_column(var = "INDV") %>% 
  filter(INDV %in% unadmixed$sampleID) %>% 
  column_to_rownames(var = "INDV") %>% 
  as.matrix()
# dplyr::select(where(~n_distinct(.) > 1)) # remove monomorphic loci

cont_admix_slow <-
  as.data.frame(cont_dos_nonas) %>% 
  rownames_to_column(var = "INDV") %>% 
  filter(INDV %in% admix_slow$sampleID) %>% 
  column_to_rownames(var = "INDV") %>% 
  as.matrix()

cont_admix_emoryi <-
  as.data.frame(cont_dos_nonas) %>% 
  rownames_to_column(var = "INDV") %>% 
  filter(INDV %in% admix_emoryi$sampleID) %>% 
  column_to_rownames(var = "INDV") %>% 
  as.matrix()

# Check ordering of mitotypes
all(rownames(cont_pure) == unadmixed_mitotypes$INDV)

### ------------ Lysosome genes ------------ 
# Import vcf for lysosome SNPs
lyso_vcf <- read.vcfR(here("data", "cz_lyso_snps.vcf"))
names_lyso <- names_helper(lyso_vcf)

# Convert to dosage matrix
lyso_dos <- vcf_to_dosage(lyso_vcf) # 222,100
rownames(lyso_dos) <- names_lyso$INDV
# Remove any columns that have NAs
lyso_dos_nonas <- lyso_dos[, !apply(is.na(lyso_dos), 2, any)] # 35,862 remain

# Subset into unadmixed and admixed inds
lyso_pure <-
  as.data.frame(lyso_dos_nonas) %>% 
  rownames_to_column(var = "INDV") %>% 
  filter(INDV %in% unadmixed$sampleID) %>% 
  column_to_rownames(var = "INDV") %>% 
  as.matrix()
# dplyr::select(where(~n_distinct(.) > 1)) # remove monomorphic loci

lyso_admix_slow <-
  as.data.frame(lyso_dos_nonas) %>% 
  rownames_to_column(var = "INDV") %>% 
  filter(INDV %in% admix_slow$sampleID) %>% 
  column_to_rownames(var = "INDV") %>% 
  as.matrix()

lyso_admix_emoryi <-
  as.data.frame(lyso_dos_nonas) %>% 
  rownames_to_column(var = "INDV") %>% 
  filter(INDV %in% admix_emoryi$sampleID) %>% 
  column_to_rownames(var = "INDV") %>% 
  as.matrix()

# Check ordering of mitotypes; should be TRUE
all(rownames(lyso_pure) == unadmixed_mitotypes$INDV)


# (6) Run DAPC ------------------------------------------------------------

### ------------ Nmts ------------ 
dp = dapc(nmts_pure,
          grp = unadmixed_mitotypes$mitotype,
          n.da = 1,
          n.pca = 4)
pred = predict.dapc(dp,
                    newdata = nmts_pure)
pred_dat = as.data.frame(pred$ind.scores)
# Add mitotype on to predictions
pred_dat <- left_join(pred_dat %>% 
                        rownames_to_column(var = "INDV"),
                      unadmixed_mitotypes)

### How well does the 'pure' DAPC predict mitotypes of admixed individuals?
pred_slow <- predict.dapc(dp,
                          newdata = nmts_admix_slow)
pred_slow = as.data.frame(pred_slow$ind.scores)
pred_emoryi <- predict.dapc(dp,
                            newdata = nmts_admix_emoryi)
pred_emoryi = as.data.frame(pred_emoryi$ind.scores)

### ------------ Control genes ------------ 
dp_cont = dapc(cont_pure,
               grp = unadmixed_mitotypes$mitotype,
               n.da = 1,
               n.pca = 4)
pred_cont = predict.dapc(dp_cont,
                    newdata = cont_pure)
pred_dat_cont = as.data.frame(pred_cont$ind.scores)
# Add mitotype on to predictions
pred_dat_cont <- left_join(pred_dat_cont %>% 
                             rownames_to_column(var = "INDV"),
                           unadmixed_mitotypes)

### How well does the 'pure' DAPC predict mitotypes of admixed individuals?
pred_slow_cont <- predict.dapc(dp_cont,
                          newdata = cont_admix_slow)
pred_slow_cont = as.data.frame(pred_slow_cont$ind.scores)
pred_emoryi_cont <- predict.dapc(dp_cont,
                            newdata = cont_admix_emoryi)
pred_emoryi_cont = as.data.frame(pred_emoryi_cont$ind.scores)

## ------------ Lysosome SNPs ------------ 
dp_lyso = dapc(lyso_pure,
               grp = unadmixed_mitotypes$mitotype,
               n.da = 1,
               n.pca = 4)
pred_lyso = predict.dapc(dp_lyso,
                         newdata = lyso_pure)
pred_dat_lyso = as.data.frame(pred_lyso$ind.scores)
# Add mitotype on to predictions
pred_dat_lyso <- left_join(pred_dat_lyso %>% 
                             rownames_to_column(var = "INDV"),
                           unadmixed_mitotypes)

### How well does the 'pure' DAPC predict mitotypes of admixed individuals?
pred_slow_lyso <- predict.dapc(dp_lyso,
                               newdata = lyso_admix_slow)
pred_slow_lyso = as.data.frame(pred_slow_lyso$ind.scores)
pred_emoryi_lyso <- predict.dapc(dp_lyso,
                                 newdata = lyso_admix_emoryi)
pred_emoryi_lyso = as.data.frame(pred_emoryi_lyso$ind.scores)


# (7) Run MWU test --------------------------------------------------------

# First, combine DAPC results from above
## ------------ Nmts ------------ 
pred_slow <- left_join(pred_slow %>% rownames_to_column(var = "INDV"), mitotypes)
pred_slow <- pred_slow %>% 
  mutate(category = "admix_slow",
         dataset = "nmt")

pred_emoryi <- left_join(pred_emoryi %>% rownames_to_column(var = "INDV"), mitotypes)
pred_emoryi <- pred_emoryi %>% 
  mutate(category = "admix_emoryi",
         dataset = "nmt")

dapc_nmts <- bind_rows(pred_slow, pred_emoryi)

## ------------ Control genes ------------ 
pred_slow_cont <- left_join(pred_slow_cont %>% rownames_to_column(var = "INDV"), mitotypes)
pred_slow_cont <- pred_slow_cont %>% 
  mutate(category = "admix_slow",
         dataset = "cont")

pred_emoryi_cont <- left_join(pred_emoryi_cont %>% rownames_to_column(var = "INDV"), mitotypes)
pred_emoryi_cont <- pred_emoryi_cont %>% 
  mutate(category = "admix_emoryi",
         dataset = "cont")
dapc_cont <- bind_rows(pred_slow_cont, pred_emoryi_cont)

## ------------ Lysosomes ------------ 
pred_slow_lyso <- left_join(pred_slow_lyso %>% rownames_to_column(var = "INDV"), mitotypes)
pred_slow_lyso <- pred_slow_lyso %>% 
  mutate(category = "admix_slow",
         dataset = "cont")

pred_emoryi_lyso <- left_join(pred_emoryi_lyso %>% rownames_to_column(var = "INDV"), mitotypes)
pred_emoryi_lyso <- pred_emoryi_lyso %>% 
  mutate(category = "admix_emoryi",
         dataset = "cont")
dapc_lyso <- bind_rows(pred_slow_lyso, pred_emoryi_lyso)

### ------------ Run MWU test on Nmts vs control genes ------------  
dapc_admixed <- bind_rows(dapc_nmts, dapc_cont)

results <- wilcox.test(LD1 ~ dataset, dapc_admixed, exact = TRUE)
results

### ------------ Run MWU test on Nmts vs lysosomes ------------  
dapc_admixed_lyso <- bind_rows(dapc_nmts, dapc_lyso)

results_lyso <- wilcox.test(LD1 ~ dataset, dapc_admixed_lyso, exact = TRUE)
results_lyso
