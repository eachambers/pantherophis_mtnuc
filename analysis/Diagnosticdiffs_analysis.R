library(tidyverse)
library(seqinr)
library(vcfR)
library(here)

## The following code calculates the admixture index (fixed differences) between individuals
## at the emoryi-slowinskii contact zone based on their mitochondrial haplotype.
## This code is modified from Chambers et al. (2023). doi: https://doi.org/10.1093/sysbio/syac056

##    FILES REQUIRED:
##          vcfs for NMTs and control loci (cznmtsnps.vcf.gz & czcontsnps.vcf.gz)
##          mitotype assignments (cz_mitotypes.txt)
##          sequence data for NMTs and control loci <- converted to fasta using vcf2phylip.py script

##    STRUCTURE OF CODE:
##              (1) Read in input files
##              (2) Calculate diagnostic diffs between reference groups
##              (3) Calculate per locus allele frequencies

# Load relevant functions
source(here("Fixed_diffs", "Fixeddiff.R"))


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
dataset = nmts # cont or nmts
dataset_name = "nmts" # "cont" or "nmts"

# If `save_file` set to TRUE, will save as e.g. "pure_allele_dict_nmts.rda"
allele_dict <- allele_dict(dataset, dataset_name, save_file = TRUE) # only need to run once; 5,954 for NMTs and 27,790 for control


# (3) Calculate per-locus allele frequencies ------------------------------

# If above has already been run, can start below; will load object as `pure_allele_dict_noambig`:
load(paste0(here("data"), "/pure_allele_dict_", dataset_name, ".rda"))

# Set a threshold value for fixed diffs between samples (relative to ref inds)
threshold = 0.5 # change to 0.5, 0.75, 0.95

# If `save_file` set to TRUE, will save file as e.g. "summary_nmts_0.5.txt" in `output_path`
results <- freq_data(dataset, 
                     dataset_name, 
                     pure_allele_dict = pure_allele_dict_noambig, 
                     threshold, 
                     save_file = FALSE, 
                     output_path = paste0(here("data"), "/"))

# Look at results
results$summary
