library(here)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

## The following code generates Fig. 3A, the Manhattan plot of GWAS results (Fig. S4), 
## calculates mean p-values per gene per chromosome for outliers, and generates 
## quantile-quantile plots (Fig. S5).

##    FILES REQUIRED:
##          Significant SNPs for NMT and control genes ("GWAS_control_sigsnps_*.txt" & "GWAS_NMT_sigsnps_*.txt"), generated using `GWAS_analysis.R`
##          Non-significant SNPs from GWAS analysis ("GWAS_results.txt"), generated using `GWAS_analysis.R`
##          `cont_gathered` and `nmt_gathered` objects from GWAS_analysis.R script

##    STRUCTURE OF CODE:
##              (1) Read in NMT and control gene data and merge them for plotting
##              (2) Build mean p-values for genes (Fig. 3A)
##              (3) Read in non-significant GWAS SNPs
##              (4) Build Manhattan plot (Fig. S2)
##              (5) Build Q-Q plot (Fig. S5)


# (1) Read in NMT and control gene data -----------------------------------

# If continuing from GWAS_analysis.R script, no need to import the following files
contdat <- read_tsv(here("data", "GWAS_control_sigsnps_0.05.txt"), col_names = TRUE)
nmtdat <- read_tsv(here("data", "GWAS_NMT_sigsnps_0.05.txt"), col_names = TRUE)

# Combine cont and nmts together so they aren't double-plotted
contdat <- contdat %>% 
  dplyr::rename(cont_cat = category)
nmtdat <- nmtdat %>% 
  dplyr::rename(nmt_cat = category)
sig_snps <- full_join(contdat, nmtdat) # 954,013 if alpha=0.05


# (2) Build figure with mean significance values per gene -----------------

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


# Plot results
dat <- bind_rows(cont_sig, nmt_sig, chrom_means)

# Which chroms have highest (>2) -logpvals? Save result to file
dat %>% 
  filter(mean_pval > 2) %>% 
  arrange(desc(mean_pval)) %>% 
  write_tsv(file = here("data", "top_sig_chroms.txt"), col_names = TRUE)

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
  write_tsv(file = here("data", "top_nmts.txt"), col_names = TRUE)

cont_gathered %>% 
  dplyr::filter(gene_group %in% top_cont$gene_group) %>% 
  dplyr::select(-pos) %>% 
  distinct() %>% 
  write_tsv(file = here("data", "top_cont.txt"), col_names = TRUE)

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


# (3) Read in data for Manhat and qq plots --------------------------------

# If continuing from GWAS_analysis.R script, no need to import the following files
gwas <- read_tsv(here("data", "GWAS_results.txt"), col_names = TRUE)
contdat <- read_tsv(paste0(here("data"), "/GWAS_control_sigsnps_", alpha, ".txt"), col_names = TRUE)
nmtdat <- read_tsv(paste0(here("data"), "/GWAS_NMT_sigsnps_", alpha, ".txt"), col_names = TRUE)

alpha = 0.01

# Combine cont and nmts together so they aren't double-plotted
contdat <- contdat %>% 
  dplyr::rename(cont_cat = category)
nmtdat <- nmtdat %>% 
  dplyr::rename(nmt_cat = category)
sig_snps <- full_join(contdat, nmtdat) %>% # 954,013 if alpha=0.05; 144,823 if alpha=0.01
  dplyr::mutate(category = case_when(nmt_cat == "NMT" ~ "NMT",
                                     cont_cat == "cont" ~ "control",
                                     nmt_cat == "non-NMT" & cont_cat == "non-control" ~ "other")) %>% 
  dplyr::select(-c(nmt_cat, cont_cat))

# Let's only look at most significantly associated mean genes (from above analysis)
mostsig <- read_tsv(here("data", "top_sig_chroms.txt"), col_names = TRUE)
mostsigchroms <- unique(mostsig$chrom)

# Extract above scaffolds from GWAS results
subset <-
  gwas %>% 
  dplyr::filter(chrom %in% mostsigchroms) # 1,846,955 obs

# Assign categories to the subset data according to outliers and categories
# We want SNPs to be categorized only if they're below set p-value and
# we want to categorize them based on whether they belong in N-mt genes, 
# control genes, or elsewhere
nonsigsubset <-
  subset %>% 
  dplyr::filter(signed.logp < -log10(alpha)) # 1,827,678 obs

sigsubset <-
  sig_snps %>% 
  dplyr::filter(chrom %in% mostsigchroms) # 19,277 obs


# (4) Build Manhattan plot ------------------------------------------------

axis_set <- nonsigsubset %>% 
  group_by(chrom) %>% 
  summarize(center = mean(pos))

# Set levels
sigsubset$chrom <- factor(sigsubset$chrom, levels = (unique(sigsubset$chrom)))

p_sig <-
  ggplot() +
  geom_point(data = sigsubset %>% filter(category == "other"), aes(x = pos, y = signed.logp), col = "thistle", size = 1.4, alpha = 0.75) +
  xlab(NULL) +
  ylab("-log(p-value)") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, (max(sigsubset$signed.logp) + 0.25))) +
  geom_point(data = sigsubset %>% filter(category != "other"), aes(x = pos, y = signed.logp, col = category), size = 1.4, alpha = 0.75) +
  scale_color_manual(values = c("NMT" = "brown", "control" = "skyblue4")) +
  geom_point(data = nonsigsubset, aes(x = pos, y = signed.logp, col = chrom), size = 1.4, alpha = 0.75) +
  scale_color_manual(values = rep(c("grey54","lightgrey"), ceiling(length(unique(sigsubset$chrom))/2))[1:length(unique(sigsubset$chrom))]) +
  theme(axis.text.x = element_text(angle = 60, size = 4, vjust = 0.5))

p <-
  ggplot() +
  geom_point(data = nonsigsubset, aes(x = pos, y = signed.logp), col = "grey", size = 1.4, alpha = 0.75) +
  ylab("-log(p-value)") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, (max(sigsubset$signed.logp) + 0.25))) +
  # scale_color_manual(values = rep(c("grey54","lightgrey"), ceiling(length(unique(nonsigsubset$chrom))/2))[1:length(unique(nonsigsubset$chrom))]) +
  # xlab(NULL) +
  # scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center) +
  # theme(axis.text.x = element_text(angle = 60, size = 4, vjust = 0.5)) +
  facet_grid(~chrom, scales = "free_x")

p +
  geom_point(data = sigsubset %>% filter(category == "other"), aes(x = pos, y = signed.logp), col = "tan1", size = 1.4, alpha = 0.75) +
  geom_point(data = sigsubset %>% filter(category != "other"), aes(x = pos, y = signed.logp, col = category), size = 1.4, alpha = 0.75) +
  scale_color_manual(values = c("NMT" = "brown", "control" = "skyblue4")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


# (5) Build qqplot --------------------------------------------------------

# Subset data to retrieve p-values for all SNPs within N-mt genes and control genes
gwas_cont <- left_join(cont_gathered, gwas) %>% 
  na.omit() # 157,943 SNPs with GWAS results
gwas_nmts <- left_join(nmt_gathered, gwas) %>% 
  na.omit() # 37,214 SNPs with GWAS results

# Build Q-Q plot
qqplot(gwas_cont$signed.logp, gwas_nmts$signed.logp)
abline(0,1, col = "red") # export 8x6

## Let's take gene-wide average p-values and re-build the qqplot
gwas_cont_avg <-
  gwas_cont %>% 
  group_by(gene_group) %>% 
  summarize(mean_signedpval = mean(signed.logp))
gwas_nmts_avg <-
  gwas_nmts %>% 
  group_by(gene_group) %>% 
  summarize(mean_signedpval = mean(signed.logp))

qqplot(gwas_cont_avg$mean_signedpval, gwas_nmts_avg$mean_signedpval)
abline(0,1, col = "red") # save 8x6
