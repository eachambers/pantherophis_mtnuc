library(tidyverse)
library(cowplot)
library(here)
library(grid)
library(gridExtra)
theme_set(theme_cowplot())

## The following code creates Figs. 3C and 3D.

##    FILES REQUIRED:
##          RData objects with allele freqs for:
##              - Control loci (freq_data_cont.rda)
##              - Nmts (freq_data_nmts.rda)
##          Diagnostic difference summary files, one for each threshold (e.g., "summary_cont_0.5.txt)
##          Top significantly associated chromosomes ("top_sig_chroms.txt") from `GWAS_figures.R` script

##    STRUCTURE OF CODE:
##              (1) Load input data
##              (2) Tidy data
##              (3) Build bar plot for each SNP (Fig. XX)
##              (4) Build summary bar plots for Nmts and control genes (Fig. XX)
##              (5) Bar plots for specific chromosomes (Fig. XX)


# (1) Load input data -----------------------------------------------------

threshold = 0.5
load(file = paste0(here("data"), "/freq_data_cont_", threshold, ".rda")) # will load as freqs
freqs_cont <- freqs
load(file = paste0(here("data"), "/freq_data_nmts_", threshold, ".rda")) # will load as freqs
freqs_nmts <- freqs


# (2) Tidy data -----------------------------------------------------------

tidy_freqs <- function(freqs) {
  tidy_dat <-
    freqs %>%
    pivot_longer(names_to = "Assignment", values_to = "fraction", cols = -c(mitotype, locus)) %>% 
    mutate(status = case_when(mitotype == "emoryi" & Assignment == "frac_emoryi" ~ "match",
                              mitotype == "emoryi" & Assignment == "frac_slow" ~ "mismatch",
                              mitotype == "emoryi" & Assignment == "frac_het" ~ "het",
                              mitotype == "slowinskii" & Assignment == "frac_slow" ~ "match",
                              mitotype == "slowinskii" & Assignment == "frac_emoryi" ~ "mismatch",
                              mitotype == "slowinskii" & Assignment == "frac_het" ~ "het")) %>% 
    group_by(locus, status) %>% 
    summarize(summed_prop = sum(fraction)) %>% 
    distinct()
  tidy_dat$status <- factor(tidy_dat$status, levels = c("match", "het", "mismatch"))
  
  return(tidy_dat)
}

tidy_cont <- tidy_freqs(freqs_cont)
tidy_nmts <- tidy_freqs(freqs_nmts)


# TODO which of the following figs did we actually keep?

# (3) Build barplot for each SNP ------------------------------------------

p_cont <-
  tidy_cont %>% 
  ggplot(aes(x = locus, y = summed_prop, fill = status)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c(het = "lightgrey", match = "#318b84", mismatch = "#b57430"),
                    labels = c(match = "Match", het = "Heterozygous", mismatch = "Mismatch")) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Control loci")

p_nmts <-
  tidy_nmts %>% 
  ggplot(aes(x = locus, y = summed_prop, fill = status)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c(het = "lightgrey", match = "#318b84", mismatch = "#b57430"),
                    labels = c(match = "Match", het = "Heterozygous", mismatch = "Mismatch")) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("NMTs")

plots <- cowplot::plot_grid(p_nmts, p_cont, ncol = 1)

# Common x and y labs
ylab <- textGrob("Proportion of individuals", 
                   gp = gpar(col = "black", fontsize = 15), rot = 90)
xlab <- textGrob("Position (SNP)", 
                   gp = gpar(col = "black", fontsize = 15))

# Add common labels to plot
grid.arrange(arrangeGrob(plots, left = ylab, bottom = xlab)) # export 10x8


# (4) Build Fig. 3C: summarized barplots -------------------------------------------

## Import the summary data
files <- intersect(list.files(here("data"), pattern = c("summary"), full.names = TRUE),
                   list.files(here("data"), pattern = c("0.5"), full.names = TRUE))

file_nos <- 1:length(files)
dat <-
  file_nos %>% 
  lapply(function(x) {
    dat <- read_tsv(files[x], col_names = TRUE)
    slow <-
      dat %>% 
      as.data.frame() %>% 
      dplyr::filter(mitotype == "slowinskii") %>% 
      dplyr::rename("match" = "slow_match", "mismatch" = "slow_mismatch", "het" = "slow_het") %>% 
      dplyr::select(mitotype, match, mismatch, het, dataset, threshold)
    emoryi <-
      dat %>% 
      as.data.frame() %>% 
      dplyr::filter(mitotype == "emoryi") %>% 
      dplyr::rename("match" = "emoryi_match", "mismatch" = "emoryi_mismatch", "het" = "emory_het") %>% 
      dplyr::select(mitotype, match, mismatch, het, dataset, threshold)
    results <- dplyr::bind_rows(slow, emoryi)
    return(results)
  }) %>% 
  dplyr::bind_rows()

dat <-
  dat %>% 
  pivot_longer(names_to = "status", values_to = "SNPs", cols = 2:4)
dat$status <- factor(dat$status, levels = c("match", "het", "mismatch"))

dat %>% 
  filter(dataset == "cont") %>% 
  ggplot(aes(x = mitotype, y = SNPs, fill = status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(het = "lightgrey", match = "#318b84", mismatch = "#b57430"),
                    labels = c(match = "Match", het = "Heterozygous", mismatch = "Mismatch")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") +
  xlab("Mitotype of individual") +
  ylab("Number of SNPs") # export 3x3

dat %>% 
  filter(dataset == "nmts") %>% 
  ggplot(aes(x = mitotype, y = SNPs, fill = status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(het = "lightgrey", match = "#318b84", mismatch = "#b57430"),
                    labels = c(match = "Match", het = "Heterozygous", mismatch = "Mismatch")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") +
  xlab("Mitotype of individual") +
  ylab("Number of SNPs") # export 3x3


# Build proportional bar plots --------------------------------------------

dat <-
  dat %>% 
  group_by(dataset, mitotype) %>% 
  mutate(total_snps = sum(SNPs),
         prop_snps = SNPs/total_snps)

# Now plot results as grouped barplot
dat %>% 
  ggplot(aes(x = dataset, y = prop_snps, fill = status, group = dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(het = "lightgrey", match = "#318b84", mismatch = "#b57430"),
                    labels = c(match = "Match", het = "Heterozygous", mismatch = "Mismatch")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") +
  xlab("Mitotype of individual") +
  ylab("Proportion of SNPs") +
  facet_grid(~mitotype) # export 3x3

# Not broken up by mitotype, still proportional:
dat %>% 
  ungroup() %>% 
  group_by(dataset, status) %>% 
  mutate(total_snps_all = sum(SNPs)) %>% 
  ungroup() %>% 
  dplyr::select(dataset, status, total_snps_all) %>% 
  distinct() %>% 
  group_by(dataset) %>% 
  mutate(sum_snps_all = sum(total_snps_all),
         prop_snps_all = total_snps_all/sum_snps_all) %>% 
  ggplot(aes(x = dataset, y = prop_snps_all, fill = status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(het = "lightgrey", match = "#318b84", mismatch = "#b57430"),
                    labels = c(match = "Match", het = "Heterozygous", mismatch = "Mismatch")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") +
  xlab("Mitotype of individual") +
  ylab("Proportion of SNPs") # export 3x3


# (5) Fig. 3D: Barplots for specific chromosomes -----------------------------------

# Read in top chromosome results from GWAS results
top_chroms <- read_tsv(here("data", "top_sig_chroms.txt"), col_names = TRUE)
chroms <- top_chroms %>% dplyr::select(chrom) %>% distinct()

cont <-
  tidy_cont %>% 
  separate(locus, into = c("chromprefix", "chromnum", "site"), sep = "_") %>% 
  unite(chrom, chromprefix:chromnum, sep = "_") %>% 
  dplyr::filter(chrom %in% chroms$chrom) %>% 
  pivot_wider(names_from = status, values_from = summed_prop)

cont <-
  cont %>% 
  arrange(chrom, site) %>% 
  dplyr::mutate(order = 1:nrow(cont)) %>% 
  pivot_longer(names_to = "status", values_to = "summed_prop", cols = -c(chrom, site, order))
cont$order <- as.factor(cont$order)

nmt <-
  tidy_nmts %>% 
  separate(locus, into = c("chromprefix", "chromnum", "site"), sep = "_") %>% 
  unite(chrom, chromprefix:chromnum, sep = "_") %>% 
  dplyr::filter(chrom %in% chroms$chrom) %>% 
  pivot_wider(names_from = status, values_from = summed_prop)

nmt <-
  nmt %>% 
  arrange(chrom, site) %>% 
  dplyr::mutate(order = 1:nrow(nmt)) %>% 
  pivot_longer(names_to = "status", values_to = "summed_prop", cols = -c(chrom, site, order))
nmt$order <- as.factor(nmt$order)

# Build the plots
p_nmt <-
  nmt %>% 
  ggplot(aes(x = order, y = summed_prop, fill = status)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 8)) +
  scale_fill_manual(values = c(het = "lightgrey", match = "#318b84", mismatch = "#b57430"),
                    labels = c(match = "Match", het = "Heterozygous", mismatch = "Mismatch")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  ggtitle("NMTs") +
  facet_wrap(~chrom, scales = "free_x", strip.position = "bottom", nrow = 1)

p_cont <-
  cont %>% 
  ggplot(aes(x = order, y = summed_prop, fill = status)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 8)) +
  scale_fill_manual(values = c(het = "lightgrey", match = "#318b84", mismatch = "#b57430"),
                    labels = c(match = "Match", het = "Heterozygous", mismatch = "Mismatch")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  ggtitle("Control loci") +
  facet_wrap(~chrom, scales = "free_x", strip.position = "bottom", nrow = 1)

plots <- cowplot::plot_grid(p_nmt, p_cont, ncol = 1)

### Combine plots together
# Common x and y labs
ylab <- textGrob("Proportion of individuals", 
                 gp = gpar(col = "black", fontsize = 15, face = "bold"), rot = 90)
xlab <- textGrob("Position (SNP)", 
                 gp = gpar(col = "black", fontsize = 15, face = "bold"))

# Add common labels to plot
grid.arrange(arrangeGrob(plots, left = ylab, bottom = xlab)) # export 7x4.5

