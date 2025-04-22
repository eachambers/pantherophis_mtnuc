library(here)
library(tidyverse)
library(vcfR)

## The following code XXX

##    FILES REQUIRED:
##          XXX

##    STRUCTURE OF CODE:
##              (1) XXX


# Input data formatting ---------------------------------------------------

ids <- read.table(here("data", "cz_snps.fam"), stringsAsFactors = FALSE)[[2]]

# Read in the vcf
vcf <- read.vcfR(here("data", "cz_snps.vcf.gz"))
# Convert to a dosage matrix
dos <- algatr::vcf_to_dosage(vcf)

# Read in PCA for covariates
pca <- read.table(here("data", "cz_snps.eigenvec")) %>% rename("SampleID" = V2)
pca$V2 <- gsub("_sorted_dedup_RG.bam", "", pca$SampleID)

# Read in phenotypes (mitotypes)
pheno <- read.table(here("data", "cz_mitotypes.txt")) %>% 
  rename("SampleID" = V1, "mitotype" = V2)

# Join with mitotypes so ordering is consistent
dat <- left_join(pca, pheno)


# Export data -------------------------------------------------------------

write.table(dat %>% dplyr::select(V3:V5), 
            here("data", "cz_snps_covar.txt"), 
            col.names = FALSE, row.names = FALSE)
write.table(dat %>% dplyr::select(mitotype), 
            here("data", "cz_snps_pheno.txt"), 
            col.names = FALSE, row.names = FALSE)

write.table(pruned_pos[c("ID","POS","CHR")],map.file,sep = " ",
            quote = FALSE,row.names = FALSE,col.names = FALSE)

