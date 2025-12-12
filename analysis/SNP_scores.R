library(vegan)
library(tidyverse)
library(cowplot)
library(here)
library(algatr)
library(vcfR)
theme_set(theme_cowplot())


# (1) Import data ---------------------------------------------------------

# vcf <- vcfR::read.vcfR(here("data", "cz_snps.vcf.gz"))
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


# (2) Calculate distances ----------------------------------------------

pca <- prcomp(gen, center = TRUE, scale. = TRUE)
head(pca$rotation)
# Extract PC1
pc1_ind <- pca$x[, 1]  # first principal component for individuals
snp_pc1 <- as.data.frame(pca$rotation[, 1]) %>% 
    rename(SNP_scores = "pca$rotation[, 1]") %>% 
    rownames_to_column("SNP")  # loadings of SNPs on PC1
abs_snp_pc1 <- abs(snp_pc1)   # often used for magnitude comparison

# Filter Nmt and control SNPs
nmt_snpscores <- snp_pc1 %>% filter(SNP %in% colnames(nmtdos))
cont_snpscores <- snp_pc1 %>% filter(SNP %in% colnames(contdos))

wilcox.test(nmt_snpscores$SNP_scores, cont_snpscores$SNP_scores, alternative = "greater")

# Try with absolute values
abs_snp_pc1 <- snp_pc1 %>% 
    mutate(abs_snpscore = abs(snp_pc1$SNP_scores))
# Filter Nmt and control SNPs
abs_nmt_snpscores <- abs_snp_pc1 %>% filter(SNP %in% colnames(nmtdos))
abs_cont_snpscores <- abs_snp_pc1 %>% filter(SNP %in% colnames(contdos))
wilcox.test(abs_nmt_snpscores$SNP_scores, abs_cont_snpscores$SNP_scores, alternative = "greater")

# Compute genotypic correlation distance among SNPs
# R <- cor(dos, use = "pairwise.complete.obs")
# gendist <- as.dist(1 - abs(cor(dos, use = "pairwise.complete.obs")))
# rownames(gendist) == rownames(metadata)

# plink --vcf cznmtcontsnps.vcf.gz --distance square ibs --out cznmtcontsnps_dist --const-fid --allow-extra-chr --autosome-num 95
# gendist <- algatr::gen_dist(plink_file = here("data/cznmtcontsnps_dist.mibs"), plink_id_file = here("data/cznmtcontsnps_dist.mibs.id"), dist_type = "plink")


# (6) Visualize results ---------------------------------------------------

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
