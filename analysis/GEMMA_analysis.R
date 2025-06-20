library(here)
library(tidyverse)
library(vcfR)
library(algatr)

## The following code creates input files for GEMMA GWAS analysis and examines results.

##    FILES REQUIRED:
##          `cz_snps.fam` or `cz_snps_ldp.fam`: fam file, generated using Plink
##          `cz_snps.vcf.gz` or `cz_snps_ldp.vcf.gz`: variants file
##          `cz_snps.eigenvec` or `cz_snps_ldp.eigenvec`: results from Plink PCA
##          `cz_mitotypes.txt` mitotype assignments
##          `cz_snps_gemma.assoc.txt` results from GEMMA analysis

##    STRUCTURE OF CODE:
##              (1) Import data files
##              (2) Export phenotype and covariate files for input into GEMMA
##              (3) Examine GEMMA results


# (1) Import data ---------------------------------------------------------

ids <- read.table(here("data", "cz_snps_ldp.fam"), stringsAsFactors = FALSE)[[2]]

# Read in the vcf
vcf <- read.vcfR(here("data", "cz_snps_ldp.vcf.gz"))
# Convert to a dosage matrix
dos <- algatr::vcf_to_dosage(vcf)

# Read in PCA for covariates
pca <- read.table(here("data", "cz_snps_ldp.eigenvec")) %>% 
  rename("SampleID" = V2) %>% 
  mutate("short_sample" = SampleID)
pca$short_sample <- gsub("_sorted_dedup_RG.bam", "", pca$short_sample)

# Read in phenotypes (mitotypes)
pheno <- read.table(here("data", "cz_mitotypes.txt")) %>% 
  rename("short_sample" = V1, "mitotype" = V2)

# Join with mitotypes so ordering is consistent
dat <- left_join(pca, pheno)

# Make a Plink-formatted phenotype file
plink_pheno <- data.frame(family = 0,
                          SampleID = dat$SampleID,
                          mitotype = dat$mitotype) %>% 
  mutate(pheno = case_when(mitotype == "emoryi" ~ 0,
                           mitotype == "slowinskii" ~ 1)) %>% 
  dplyr::select(-mitotype)


# (2) Export data files ---------------------------------------------------

write.table(dat %>% dplyr::select(V3:V5), 
            here("data", "cz_snps_covar.txt"), 
            col.names = FALSE, row.names = FALSE)
write.table(dat %>% dplyr::select(mitotype), 
            here("data", "cz_snps_pheno.txt"), 
            col.names = FALSE, row.names = FALSE)
write.table(plink_pheno, here("data", "plink_pheno.txt"),
            col.names = FALSE, row.names = FALSE)


# (3) Examine GEMMA results -----------------------------------------------

# Read in GEMMA results
gwscan <- read.table(here("data", "output", "cz_snps_gemma.assoc.txt"), as.is = "rs", header = TRUE)

# Plot inflation
plot.inflation(gwscan$p_lrt)

# Look at raw (unadjusted) p-values
min(gwscan$p_wald)
length(which(gwscan$p_wald < 0.05))

# Set a p-value threshold based on number of observations (identical to Bonferroni correction I think)
thresh <- 0.05/nrow(gwscan)
# Now see how many outliers there are below that threshold
res <- gwscan %>% filter(p_wald < thresh)
length(res)

# Can also determine outliers using the p.adjust function with various methods
# Bonferroni correction (should be the same as above)
pval_bonf <- p.adjust(gwscan$p_lrt, method = "bonferroni")
length(which(pval_bonf < 0.05))

# False discovery rate correction
pval_fdr <- p.adjust(gwscan$p_score, method = "fdr")
length(which(pval_fdr < 0.05))

