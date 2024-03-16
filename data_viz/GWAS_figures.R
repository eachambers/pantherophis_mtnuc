library(here)
library(tidyverse)
library(cowplot)
library(gt)
library(gtExtras)
theme_set(theme_cowplot())

## The following code generates Fig. XX, the Manhattan plot of GWAS results, as well as
## mean p-values per gene per chromosome for outliers.

##    FILES REQUIRED:
##          Significant SNPs for NMT and control genes ("GWAS_control_sigsnps_0.05.txt" & "GWAS_NMT_sigsnps_0.05.txt"), generated using `GWAS_analysis.R`
##          Non-significant SNPs from GWAS analysis ("GWAS_results.txt"), generated using `GWAS_analysis.R`
##          `cont_gathered` and `nmt_gathered` from GWAS_analysis.R script

##    STRUCTURE OF CODE:
##              (1) Read in NMT and control gene data and merge them for plotting
##              (2) Read in non-significant GWAS SNPs
##              (3) Build Manhattan plot (Fig. XX)
##              (4) Build nice table of results
##              (5) Build mean p-values for genes (Fig. XX)

setwd("~/Box Sync/Lampropeltis project/PANTHEROPHIS_introgression")


# (1) Read in NMT and control gene data -----------------------------------

# If continuing from GWAS_analysis.R script, no need to import the following files
contdat <- read_tsv("GWAS/GWAS_control_sigsnps_0.05.txt", col_names = TRUE)
nmtdat <- read_tsv("GWAS/GWAS_NMT_sigsnps_0.05.txt", col_names = TRUE)

# Combine cont and nmts together so they aren't double-plotted
contdat <- contdat %>% 
  dplyr::rename(cont_cat = category)
nmtdat <- nmtdat %>% 
  dplyr::rename(nmt_cat = category)
sig_snps <- full_join(contdat, nmtdat) # 954,013 if alpha=0.05


# (2) Read in non-sig GWAS SNPs -------------------------------------------

# If continuing from GWAS_analysis.R script, no need to import the following file
dat <- read_tsv("GWAS/GWAS_results.txt", col_names = TRUE)
alpha = 0.05
nonsig_snps <-
  dat %>% 
  dplyr::filter(signed.logp < -log10(alpha)) # 25,263,575 at alpha=0.05

# Now, categorize nonsig snps that are in NMTs and control genes
# `nmt_gathered` and `cont_gathered` are from GWAS_analysis.R script
nonsig_snps <- left_join(nonsig_snps, nmt_gathered) %>% 
  dplyr::select(chrom, pos, zscore, pos.Mb, signed.logp, r2, category) %>%
  distinct() %>%
  replace_na(list(category = "non-NMT")) %>% 
  dplyr::rename(nmt_cat = category)

nonsig_snps <- left_join(nonsig_snps, cont_gathered) %>% 
  dplyr::select(chrom, pos, zscore, pos.Mb, signed.logp, r2, nmt_cat, category) %>%
  distinct() %>%
  replace_na(list(category = "non-control")) %>% 
  dplyr::rename(cont_cat = category)

# Get stats
nonsig_snps %>% dplyr::filter(nmt_cat == "NMT") %>% count() # 36,045 at alpha=0.05
nonsig_snps %>% dplyr::filter(cont_cat == "cont") %>% count() # 152,351 at alpha=0.05

# Verify that there are no SNPs that are categorized as both NMT and controls:
nonsig_snps %>% dplyr::filter(nmt_cat == "NMT" & cont_cat == "cont") # should have 0 rows


# (3) Build Manhattan plot ------------------------------------------------

# Determine max y value
max_y = max(sig_snps$signed.logp)
thresh_y = max_y+max_y*0.05
min_y = min(nonsig_snps$signed.logp)

# For simplicity, let's only plot a single chromosome that contains both N-mts and control outliers:
sig_stats <-
  sig_snps %>% 
  dplyr::filter(cont_cat == "cont" | nmt_cat == "NMT") %>% 
  group_by(chrom, nmt_cat, cont_cat) %>% 
  count()

# Combine cat cols into a single col
sig_snps <-
  sig_snps %>% 
  dplyr::mutate(status = case_when(nmt_cat == "NMT" ~ "sig_NMT",
                                   cont_cat == "cont" ~ "sig_cont",
                                   nmt_cat == "non-NMT" & cont_cat == "non-control" ~ "sig_neither"))
nonsig_snps <-
  nonsig_snps %>% 
  dplyr::mutate(status = case_when(nmt_cat == "NMT" ~ "nonsig_NMT",
                                   cont_cat == "cont" ~ "nonsig_cont",
                                   nmt_cat == "non-NMT" & cont_cat == "non-control" ~ "nonsig_neither"))

# Interesting chroms for us to look at: NW_023010694.1, NW_023010706.1, NW_023010796.1, NW_023010908.1
chrom = "NW_023010694.1"

p <-
  sig_snps %>% 
  dplyr::filter(chrom == chrom) %>% 
  ggplot(aes(x = pos.Mb, y = signed.logp, color = status)) +
  geom_point(size = 0.25, alpha = 0.25) +
  ylab("-log10(p-value)") +
  xlab("position (Mb)") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_y_continuous(expand = c(0,0), limits = c(0, thresh_y)) +
  scale_x_continuous(expand = c(0,0))

p +
  geom_point(data = nonsig_snps %>% 
               dplyr::filter(chrom == chrom), size = 0.25, alpha = 0.25) +
  ylab("-log10(p-value)") +
  xlab("position (Mb)") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black", linewidth = 0.6) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("sig_neither" = "darkgrey", "sig_NMT" = "brown", "sig_cont" = "skyblue4",
                                "nonsig_neither" = "lightgrey", "nonsig_NMT" = "brown2", "nonsig_cont" = "skyblue1")) # export 6x3


# (4) Build nice table of sig SNPs ----------------------------------------

d <- max(abs(min(sig_snps$signed.logp)), abs(max(sig_snps$signed.logp)))

sig_snps_subset <-
  sig_snps %>%
  dplyr::filter(chrom %in% cont_sig$chrom | chrom %in% nmt_sig$chrom) %>% 
  arrange(desc(signed.logp)) %>% 
  top_n(n = 50, wt = signed.logp)

min <- abs(min(sig_snps_subset$signed.logp))
max <- abs(max(sig_snps_subset$signed.logp))
sig_snps_subset %>% 
  gt::gt() %>% 
  gtExtras::gt_hulk_col_numeric(signed.logp, trim = TRUE, domain = c(-max, max))


# (5) Build figure with mean significance values per gene -----------------

# Retrieve p-values for Nmts and max per-gene mean p-value
nmt_vals <-
  sig_snps %>% 
  dplyr::filter(status == "sig_NMT") %>% 
  left_join(nmt_gathered) # 771 SNPs; 54 chroms

nmt_sig <-
  nmt_vals %>% 
  group_by(chrom, gene_group) %>% 
  summarize(mean_pval = mean(signed.logp)) %>% 
  mutate(dataset = "nmt")

cont_vals <-
  sig_snps %>% 
  dplyr::filter(status == "sig_cont") %>% 
  left_join(cont_gathered) # 4838 SNPs; 64 chroms

cont_sig <-
  cont_vals %>% 
  group_by(chrom, gene_group) %>% 
  summarize(mean_pval = mean(signed.logp)) %>% 
  mutate(dataset = "cont")

# Get mean (chromosome-wide) p-values from chroms that have both Nmts and cont genes
chrom_means <- sig_snps %>% 
  dplyr::filter(status == "sig_neither") %>% 
  group_by(chrom) %>% 
  summarize(mean_pval = mean(signed.logp)) %>% 
  # Only consider chroms that contain either Nmt or control genes (or both)
  dplyr::filter(chrom %in% cont_sig$chrom | chrom %in% nmt_sig$chrom) %>%
  mutate(dataset = "other")


# Plot results ------------------------------------------------------------

dat <- bind_rows(cont_sig, nmt_sig, chrom_means)

# Which chroms have highest (>2) -logpvals? Save result to file
dat %>% 
  filter(mean_pval > 2) %>% 
  arrange(desc(mean_pval)) %>% 
  write_tsv(file = "GWAS/top_sig_chroms.txt", col_names = TRUE)

top_nmts <- dat %>% 
  filter(mean_pval > 2) %>% 
  arrange(desc(mean_pval)) %>% 
  dplyr::filter(dataset == "nmt")

top_cont <- dat %>% 
  filter(mean_pval > 2) %>% 
  arrange(desc(mean_pval)) %>% 
  dplyr::filter(dataset == "cont")

# What genes have the strongest associations?
nmt_gathered %>% 
  dplyr::filter(gene_group %in% top_nmts$gene_group) %>% 
  dplyr::select(-pos) %>% 
  distinct() %>% 
  write_tsv(file = "GWAS/top_nmts.txt", col_names = TRUE)

cont_gathered %>% 
  dplyr::filter(gene_group %in% top_cont$gene_group) %>% 
  dplyr::select(-pos) %>% 
  distinct() %>% 
  write_tsv(file = "GWAS/top_cont.txt", col_names = TRUE)

dat %>% 
  ggplot(aes(x = chrom, y = mean_pval, color = dataset)) +
  theme(panel.grid.major = element_line(color = "lightgrey",
                                        linewidth = 0.25),
        # axis.text.x = element_text(angle = 90, size = 6),
        axis.text.x = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.position = "none") +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c("other" = "darkgrey", "nmt" = "brown", "cont" = "skyblue4")) +
  ylab("Mean -log10(p-value)") +
  xlab("Scaffold") # export 10x5
