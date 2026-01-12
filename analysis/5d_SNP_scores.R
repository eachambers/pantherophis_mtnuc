library(vegan)
library(tidyverse)
library(cowplot)
library(here)
library(algatr)
library(vcfR)
theme_set(theme_cowplot())

## The following code generates SNP scores along PC1 and conducts a MWU test
## to assess whether there is an absolute loading difference between N-mt SNPs 
## and those in control genes.

##    FILES REQUIRED:
##          N-mt and control gene SNPs for samples at contact zone ("cznmtsnps.vcf" and "czcontsnps2.vcf.gz")
##          Ordering of samples in vcf ("bams")
##          Assignment of samples to mitotype ("cz_mitotypes.txt")

##    STRUCTURE OF CODE:
##              (1) Import data
##              (2) Run PCA and extract loadings
##              (3) Visualize results (Fig. 3A)


# (1) Import data ---------------------------------------------------------

nmtvcf <- vcfR::read.vcfR(here("data", "cznmtsnps.vcf"))
contvcf <- vcfR::read.vcfR(here("data", "czcontsnps2.vcf.gz"))

nmtdos <- algatr::vcf_to_dosage(nmtvcf)
contdos <- algatr::vcf_to_dosage(contvcf)
# Combine dosage matrices
dos <- bind_cols(as.data.frame(nmtdos), as.data.frame(contdos))
# Impute to the median
gen <- algatr::simple_impute(dos)

bams = read.table(here("data/bams"), col.names = "filename")
bams <- 
  bams %>% 
  mutate(bamname = filename) %>% 
  separate(filename, into = c("sampleID", "tmp1", "tmp2", "tmp3"), sep = "_") %>% 
  dplyr::select(bamname, sampleID)

mitotypes <- read_tsv(here("data", "cz_mitotypes.txt"), col_names = c("sampleID", "mitotype"))

# Merge bams with mitotypes
metadata <- left_join(bams, mitotypes) %>% dplyr::select(bamname, mitotype)


# (2) Run PCA and extract loadings ----------------------------------------

pca <- prcomp(gen, center = TRUE, scale. = TRUE)
head(pca$rotation)
# Extract PC1
pc1_ind <- pca$x[, 1]  # first principal component for individuals
snp_pc1 <- as.data.frame(pca$rotation[, 1]) %>% 
    rename(SNP_scores = "pca$rotation[, 1]") %>% 
    rownames_to_column("SNP")  # loadings of SNPs on PC1
abs_snp_pc1 <- abs(snp_pc1)   # often used for magnitude comparison

# Filter Nmt and control SNPs
abs_snp_pc1 <- snp_pc1 %>% 
    mutate(abs_snpscore = abs(snp_pc1$SNP_scores))
# Filter Nmt and control SNPs
abs_nmt_snpscores <- abs_snp_pc1 %>% filter(SNP %in% colnames(nmtdos))
abs_cont_snpscores <- abs_snp_pc1 %>% filter(SNP %in% colnames(contdos))
wilcox.test(abs_nmt_snpscores$SNP_scores, abs_cont_snpscores$SNP_scores, alternative = "greater")


# (3) Visualize results ---------------------------------------------------

# Two categories on x-axis, SNP scores on y
p_scores <- 
abs_snp_pc1 %>% 
    mutate(snp_set = case_when(SNP %in% colnames(nmtdos) ~ "Nmt",
                                SNP %in% colnames(contdos) ~ "Control")) %>% 
    ggplot(aes(x = snp_set, y = abs_snpscore)) + # SNP_scores
    geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_minimal() +
  ylab("PC1 SNP loading") +
  xlab("SNP set") +
  ggtitle("Comparison of SNP loadings (Mann–Whitney U test)")

p_scores
ggsave(here("outputs/SNP_scores_PC1_abs.pdf"), width = 9, height = 6.3)
