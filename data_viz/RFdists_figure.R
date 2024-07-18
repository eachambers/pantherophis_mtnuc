library(tidyverse)
library(here)
library(cowplot)
library(ggtree)

## The following code takes in results from RF distances analysis and plots stacked barplots 
## based on proportions of N-mt and control gene trees that match to the species, 
## mitochondrial, or minor trees (Fig. 2C).

##    FILES REQUIRED:
##          `cont_results` and `nmt_results` objects from RFdists_analysis.R script

##    STRUCTURE OF CODE:
##              (1) Process data
##              (2) Build plot


# (1) Process data --------------------------------------------------------

nmt_plot <-
  nmt_results %>% 
  rename("minor" = "RF_nmt_minor", # TODO CHANGE INPUT FILE NAMES
         "mito" = "RF_nmt_mito",
         "spp" = "RF_nmt_spp") %>% 
  pivot_longer(names_to = "comp_tree", values_to = "RF_dist", cols = -c(order, avg_bs)) %>% 
  group_by(comp_tree) %>% 
  dplyr::filter(RF_dist == 1) %>% 
  summarize(nmt_matches = n())

cont_plot <-
  cont_results %>% 
  rename("minor" = "RF_cont_minor", # TODO CHANGE INPUT FILE NAMES
         "mito" = "RF_cont_mito",
         "spp" = "RF_cont_spp") %>% 
  pivot_longer(names_to = "comp_tree", values_to = "RF_dist", cols = -c(order, avg_bs)) %>% 
  group_by(comp_tree) %>% 
  dplyr::filter(RF_dist == 1) %>% 
  summarize(cont_matches = n())

plot_dat <- left_join(nmt_plot, cont_plot) %>% # joins by comp_tree
  pivot_longer(names_to = "dataset",
               values_to = "no_matches",
               cols = -c(comp_tree))


# (2) Build plot ----------------------------------------------------------

plot_dat %>% 
  ggplot(aes(x = dataset, y = no_matches, fill = comp_tree)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.75) +
  ylab("Number of matching topologies") +
  scale_fill_manual(values = c("mito" = "#6b8577", "minor" = "#2b5ba4", "spp" = "#95627a")) +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("Control genes", "N-mt genes")) # export 6x6
