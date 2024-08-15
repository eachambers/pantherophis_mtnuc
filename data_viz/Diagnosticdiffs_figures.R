library(tidyverse)
library(cowplot)
library(here)
library(grid)
library(gridExtra)
theme_set(theme_cowplot())

## The following code creates Figs. 3C and 3D.

##    FILES REQUIRED:
##          RData objects with allele freqs for:
##              - Control loci (freq_data_cont_ldp.rda)
##              - Nmts (freq_data_nmts_ldp.rda)
##          Diagnostic difference summary files, one for each threshold (e.g., "summary_cont_ldp_0.5.txt)
##          Top significantly associated chromosomes ("top_sig_chroms.txt") from `GWAS_figures.R` script

##    STRUCTURE OF CODE:
##              (1) Load input data
##              (2) Tidy data
##              (3) Build summary bar plots for Nmts and control genes (Fig. 3C)
##              (4) Bar plots for specific chromosomes (Fig. 3D)


# (1) Load input data -----------------------------------------------------

threshold = 0.5
load(file = paste0(here("data"), "/freq_data_cont_ldp_", threshold, ".rda"))
load(file = paste0(here("data"), "/freq_data_nmts_ldp_", threshold, ".rda"))


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

tidy_cont <- tidy_freqs(freqs_to_plot_cont)
tidy_nmts <- tidy_freqs(freqs_to_plot_nmts)


# (3) Fig. 3C: summarized barplots ----------------------------------------

## Import the summary data
files <- intersect(list.files(here("data"), pattern = c("summary"), full.names = TRUE),
                   list.files(here("data"), pattern = c("0.5"), full.names = TRUE))

dat <-
  1:length(files) %>% 
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

dattoplot <-
  dat %>% 
  pivot_longer(cols = match:het, names_to = "category", values_to = "number_sites") %>% 
  group_by(dataset, category) %>% 
  mutate(total_number_sites = sum(number_sites)) %>% 
  dplyr::select(-mitotype) %>% 
  ungroup() %>% 
  group_by(dataset) %>% 
  mutate(prop_sites = total_number_sites/sum(number_sites)) %>% 
  dplyr::select(-number_sites) %>% 
  distinct()

dattoplot$category <- factor(dattoplot$category, levels = c("mismatch", "het", "match"))

dattoplot %>% 
  ggplot(aes(x = dataset, y = prop_sites, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(het = "lightgrey", match = "#318b84", mismatch = "#b57430"),
                    labels = c(match = "Match", het = "Heterozygous", mismatch = "Mismatch")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none",
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  xlab("Gene category") +
  ylab("Proportion of SNPs") # export 4x4


# (4) Fig. 3D: Barplots for specific chromosomes -----------------------------------

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

