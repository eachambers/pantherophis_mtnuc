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
sig_snps <- full_join(contdat, nmtdat) # 954,013 if alpha=0.05; 144,823 if alpha=0.01


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

# We're only going to build a Manhattan plot for the top ten scaffolds that contain
# either N-mt or control genes
mostsig <- read_tsv(here("data", "top_sig_chroms.txt"), col_names = TRUE)
mostsigchroms <- unique(mostsig$chrom)

# Extract most sig scaffolds from GWAS results
subset <-
  gwas %>% 
  dplyr::filter(chrom %in% mostsigchroms) # 1,846,955 obs

alpha = "0.01"
contdat <- read_tsv(paste0(here("data"), "/GWAS_control_sigsnps_", alpha, ".txt"), col_names = TRUE)
nmtdat <- read_tsv(paste0(here("data"), "/GWAS_NMT_sigsnps_", alpha, ".txt"), col_names = TRUE)

# Combine cont and nmts together so they aren't double-plotted
contdat <- contdat %>% 
  dplyr::rename(cont_cat = category)
nmtdat <- nmtdat %>% 
  dplyr::rename(nmt_cat = category)
sig_snps <- full_join(contdat, nmtdat) %>% # 954,013 if alpha=0.05; 144,823 if alpha=0.01
  dplyr::mutate(Category = case_when(nmt_cat == "NMT" ~ "N-mt",
                                     cont_cat == "cont" ~ "Control",
                                     nmt_cat == "non-NMT" & cont_cat == "non-control" ~ "Outlier")) %>% 
  dplyr::select(-c(nmt_cat, cont_cat))

# Assign categories to the subset data according to outliers and categories
# We want SNPs to be categorized only if they're below set p-value and
# we want to categorize them based on whether they belong in N-mt genes, 
# control genes, or elsewhere
alpha = 0.01 # needs to be numeric now
nonsigsubset <-
  subset %>% 
  dplyr::filter(signed.logp < -log10(alpha)) %>% 
  dplyr::select(chrom, pos, signed.logp) %>% 
  mutate(Category = "Non-outlier") # 1,827,678 obs

sigsubset <-
  sig_snps %>% 
  dplyr::filter(chrom %in% mostsigchroms) %>% 
  dplyr::select(chrom, pos, signed.logp, Category) # 19,277 obs

# We need to ensure ordering is correct so we're going to combine
# datasets together for ease
plot_dat <- bind_rows(nonsigsubset, sigsubset)
plot_dat <-
  plot_dat %>% 
  arrange(chrom, pos) %>% 
  mutate(full_pos = 1:nrow(plot_dat))


# (4) Build Manhattan plot ------------------------------------------------

# Get center for each scaffold for adding breaks to x axis
axis_set <- plot_dat %>% 
  group_by(chrom) %>% 
  summarize(center = mean(full_pos))

# Set levels
plot_dat$chrom <- factor(plot_dat$chrom, levels = (unique(plot_dat$chrom)))

# Sloppy way of doing this but it works
odd_chroms <- c("NW_023010710.1", "NW_023010746.1", "NW_023010870.1", "NW_023012548.1", "NW_023036219.1")
even_chroms <- c("NW_023010730.1", "NW_023010793.1", "NW_023010881.1", "NW_023031296.1", "NW_023041095.1")

# Build Manhattan plot
ggplot() +
  geom_scattermore(data = plot_dat %>% 
                     filter(Category == "Outlier"), 
                   aes(x = full_pos, y = signed.logp), 
                   col = "goldenrod1", size = 1.4, alpha = 0.75) +
  xlab(NULL) +
  ylab("-log(p-value)") +
  scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center, expand = c(0, 0), guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, (max(plot_dat$signed.logp) + 0.25))) +
  geom_point(data = plot_dat %>% 
               filter(Category == "N-mt" | Category == "Control"), 
             aes(x = full_pos, y = signed.logp, col = Category), size = 1.4, alpha = 0.75) +
  scale_color_manual(values = c("N-mt" = "brown", "Control" = "skyblue4")) +
  geom_scattermore(data = plot_dat %>% 
                     filter(Category == "Non-outlier") %>% 
                     filter(chrom %in% odd_chroms),
                   aes(x = full_pos, y = signed.logp), col = "lightgrey", size = 1.4, alpha = 0.75) +
  geom_scattermore(data = plot_dat %>% 
                     filter(Category == "Non-outlier") %>% 
                     filter(chrom %in% even_chroms),
                   aes(x = full_pos, y = signed.logp), col = "grey54", size = 1.4, alpha = 0.75) +
  # scale_color_manual(values = rep(c("grey54","lightgrey"), ceiling(length(unique(nonsigsubset$chrom))/2)[1:length(unique(nonsigsubset$chrom))]) +
  theme(axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5)) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black", linewidth = 0.6) # export 10x6


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
