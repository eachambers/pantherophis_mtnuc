library(tidyverse)
library(seqinr)
library(vcfR)
library(here)

## The following code calculates the admixture index (fixed differences) between individuals
## at the emoryi-slowinskii contact zone based on their mitochondrial haplotype.
## This code is modified from Chambers et al. (2023). doi: https://doi.org/10.1093/sysbio/syac056

##    FILES REQUIRED:
##          vcfs for NMTs and control loci (cznmtsnps.vcf & czcontsnps.vcf)
##          fastas for NMTs and control loci (cznmtsnps.min4.fasta & czcontsnps.min4.fasta) <- vcfs converted to fasta using vcf2phylip.py
##          mitotype assignments (cz_mitotypes.txt)
##          vcfs for LD-pruned NMTs and control loci (cznmtsnps_ldp.vcf & czcontsnps_ldp.vcf) <- generated using `ld_pruning.sh` script
##          fastas for LD-pruned NMTs and control loci (cznmtsnps_ldp.min1.fasta & czcontsnps_ldp.min1.fasta) <- vcfs converted to fasta using vcf2phylip.py

##    STRUCTURE OF CODE:
##              (1) Read in input files
##              (2) Calculate diagnostic diffs between reference groups
##              (3) Calculate per locus allele frequencies
##              (4) Individual-based statistical test

# Load relevant functions
source(here("analysis", "Diagnosticdiffs.R"))


# (1) Read in input files -------------------------------------------------

### Import mitotype assignments for 30 individuals
mitotypes <- read_tsv(here("data", "cz_mitotypes.txt"), col_names = FALSE) %>% 
  rename(INDV = X1, mitotype = X2)

### Run above function for NMTs and control loci
nmts <- process_seq(vcf_file = here("data", "cznmtsnps.vcf"),
                    fasta_file = here("data", "cznmtsnps.min4.fasta"),
                    mitotypes = mitotypes) # 30 inds of 38,553 vars (38,551 SNPs)
cont <- process_seq(vcf_file = here("data", "czcontsnps.vcf"),
                    fasta_file = here("data", "czcontsnps.min4.fasta"),
                    mitotypes = mitotypes) # 30 inds of 163,960 vars (163,958 SNPs)


# (2) Calculate diagnostic diffs between reference inds -------------------

# Determine fixed differences between reference groups
dataset = cont # cont or nmts
dataset_name = "cont" # "cont" or "nmts"
# If `save_file` set to TRUE, will save as e.g. "pure_allele_dict_nmts.rda"
pure_allele_dict_noambig_cont <- allele_dict(dataset, dataset_name, save_file = FALSE) # only need to run once; 5,954 for NMTs and 27,790 for control

dataset = nmts # cont or nmts
dataset_name = "nmts" # "cont" or "nmts"
# If `save_file` set to TRUE, will save as e.g. "pure_allele_dict_nmts.rda"
pure_allele_dict_noambig_nmts <- allele_dict(dataset, dataset_name, save_file = FALSE) # only need to run once; 5,954 for NMTs and 27,790 for control


# (3) Calculate per-locus allele frequencies ------------------------------

# TODO none of this is no longer necessary because of LD-pruning
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
results <- freq_data(dataset = cont, 
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
                        fasta_file = here("data", "cznmtsnps_ldp.min1.fasta"),
                        mitotypes = mitotypes) # 30 inds, 14,330 SNPs (+INDV +mitotype cols)
cont_ldp <- process_seq(vcf_file = here("data", "czcontsnps_ldp.vcf"),
                        fasta_file = here("data", "czcontsnps_ldp.min1.fasta"),
                        mitotypes = mitotypes) # 30 inds, 61,697 SNPs (+INDV +mitotype cols)

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

summary(mod) # gene_setnmts are statistically significant (p-value = 3.49e-07)

# Odds ratio
exp(coef(mod)) # 1.074336 for gene_setnmts; i.e., 7% increase of matches in N-mts compared to control genes

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


# ===========================================================================
# Take a look at the per-individual differences in proportions
ind_stats <- ind_stats(nmts_ldp, 
                       cont_ldp, 
                       pure_allele_dict_nmts = pure_allele_dict_noambig_nmts, 
                       pure_allele_dict_cont = pure_allele_dict_noambig_cont)

ind_stats %>% 
  ggplot(aes(x = prop_diff_ind, fill = mitotype)) +
  geom_bar(stat = "count") +
  annotate("text", x = 0.06, y = 1.5, size = 8, label = "More matches in nmts") +
  annotate("text", x = -0.06, y = 1.5, size = 8, label = "More matches in controls") +
  scale_fill_manual(values = c("#eda157", "#1743d7"))

hist(ind_stats$prop_diff_ind)
t.test(ind_stats$prop_diff_ind)
# nmts = "#b35455", cont = "#6c8a9d"
# emoryi = "#eda157", slow = "#1743d7"

# Visualize another way
ind_stats %>% 
  ggplot(aes(x = reorder(INDV, prop_diff_ind, sum), y = prop_diff_ind, fill = mitotype)) +
  geom_bar(stat = "identity") +
  xlab("Individual") +
  ylab("Difference between proportion of matches \nbetween Nmts and control genes") +
  scale_fill_manual(values = c("#eda157", "#1743d7")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 17/30 individuals had props > 0; 5 of these inds had slow mitotypes and 12 had emoryi mitotypes
ind_stats %>% 
  filter(prop_diff_ind > 0) %>% 
  nrow()

# ind_stats %>% 
#   filter(mitotype == "emoryi") %>% 
#   filter(prop_diff_ind > 0) %>% 
#   nrow()
