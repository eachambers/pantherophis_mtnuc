library(tidyverse)
library(here)
library(phangorn)
library(ape)
library(cowplot)

## The following code takes in gene trees from two sources: (1) N-mt genes (n=167) and 
## (2) control genes (n=142) and compares them (by calculating Robinson-Foulds [RF] distances 
## and average bootstrap support values) to three topologies: (1) the species tree, 
## (2) the mitochondrial tree, and (3) the ILS tree. It also calculates weighted RF distances 
## (i.e., taking into account branch lengths).

## Some of the following code (and functions, in `RFdists.R`), particularly code to order 
## loci according to RFs, comes directly from the genesortR package available here: 
## https://github.com/mongiardino/genesortR.

##    FILES REQUIRED:
##          Gene trees from N-mt genes ("nmt_3taxa_genetrees.txt") and control genes 
##          ("cont_3taxa_genetrees.txt")
##          Topologies for species, mitochondrial, and minor trees ("panther_speciestree.txt", 
##          "panther_mitotree.txt", and "panther_minortree.txt")

##    STRUCTURE OF CODE:
##              (1) Read in input files
##              (2) Calculate average bootstrap values
##              (3) Calculate RF distances
##              (4) Join together results and get summary stats
##              (5) Logistic regression

# Load relevant functions
source(here("analysis", "RFdists.R"))


# (1) Read in input files -------------------------------------------------

# The following will produce an object that contains (multi)phylo objects for each of the N-mt gene trees (`nmt_gt_file` 
# arg for file name), control gene trees (`cont_gt_file` arg for file name; default is NULL), the species tree,  
# the mitochondrial tree, and the minor tree. The function assumes that the species, mito, and minor trees have the same prefix 
# (specified using `prefix`) and are suffixed with "_speciestree.txt", "_mitotree.txt", and "_minortree.txt".
dat <- read_files(data_path = paste0(here("data"), "/"), 
                  nmt_gt_file = "nmt_3taxa_genetrees.txt", 
                  cont_gt_file = "cont_3taxa_genetrees.txt", 
                  prefix = "panther")


# (2) Calculate average bootstrap values ----------------------------------

# Calculate average bootstrap values for N-mt gene trees
bs_nmt <- calc_avg_bs(dat[["nmt_gt"]])
bs_nmt %>% dplyr::filter(avg_bs >= 80) %>% tally() # 133 of 167 trees have avg bs >= 80%

# Calculate average bootstrap values for N-mt gene trees
bs_cont <- calc_avg_bs(dat[["cont_gt"]])
bs_cont %>% dplyr::filter(avg_bs >= 80) %>% tally() # 125 of 142 trees have avg bs >= 80%


# (3) Calculate RF distances and average bootstrap values -----------------

# The `calc_rfs()` function calculates RF distances given two sets of trees, specified using the `comp_gt` and `comp_to` args. 
# If `weighted = TRUE` (defaults to FALSE), function will calculate weighted RF distances (i.e., taking into account branch
# lengths). Be sure your trees are rooted; there is no check for this within the function. Can check by doing so:
all(ape::is.rooted(dat[["nmt_gt"]]))
all(ape::is.rooted(dat[["cont_gt"]]))

RF_nmt_spp <- calc_rfs(comp_gt = dat[["nmt_gt"]], comp_to = dat[["spptree"]], weighted = FALSE)
RF_nmt_mito <- calc_rfs(comp_gt = dat[["nmt_gt"]], comp_to = dat[["mitotree"]], weighted = FALSE)
RF_nmt_minor <- calc_rfs(comp_gt = dat[["nmt_gt"]], comp_to = dat[["minortree"]], weighted = FALSE)

RF_cont_spp <- calc_rfs(comp_gt = dat[["cont_gt"]], comp_to = dat[["spptree"]], weighted = FALSE)
RF_cont_mito <- calc_rfs(comp_gt = dat[["cont_gt"]], comp_to = dat[["mitotree"]], weighted = FALSE)
RF_cont_minor <- calc_rfs(comp_gt = dat[["cont_gt"]], comp_to = dat[["minortree"]], weighted = FALSE)


# (4) Join results and get summary stats ----------------------------------

# Combine results together
nmt_spp <- turn_tibble_helper(RF_nmt_spp) %>% 
  dplyr::rename(RF_nmt_spp = value)
nmt_mito <- turn_tibble_helper(RF_nmt_mito) %>% 
  dplyr::rename(RF_nmt_mito = value)
nmt_minor <- turn_tibble_helper(RF_nmt_minor) %>% 
  dplyr::rename(RF_nmt_minor = value)

nmt_results <- full_join(nmt_spp, nmt_mito, by = "order") # join based on the order column
nmt_results <- full_join(nmt_results, nmt_minor, by = "order") # join based on the order column
nmt_results <- full_join(nmt_results, bs_nmt %>% rename("order" = tree_no), by = "order") # join based on the order column

# Remove trees with avg bootstrap < 80%
nmt_results <- nmt_results %>% 
  dplyr::filter(avg_bs >= 80) # 133 obs

### Do the same for control genes
cont_spp <- turn_tibble_helper(RF_cont_spp) %>% 
  dplyr::rename(RF_cont_spp = value)
cont_mito <- turn_tibble_helper(RF_cont_mito) %>% 
  dplyr::rename(RF_cont_mito = value)
cont_minor <- turn_tibble_helper(RF_cont_minor) %>% 
  dplyr::rename(RF_cont_minor = value)

cont_results <- full_join(cont_spp, cont_mito, by = "order") # join based on the order column
cont_results <- full_join(cont_results, cont_minor, by = "order") # join based on the order column
cont_results <- full_join(cont_results, bs_cont %>% rename("order" = tree_no), by = "order") # join based on the order column

# Remove trees with avg bootstrap < 80%
cont_results <- cont_results %>% 
  dplyr::filter(avg_bs >= 80) # 125 obs

#### Get summary statistics
# How many matches are there in N-mts?
nmt_results %>% dplyr::filter(RF_nmt_spp == 1) %>% summarize(matches = n(),
                                                             prop_gts = matches/nrow(nmt_results)*100) # 66.2%
nmt_results %>% dplyr::filter(RF_nmt_mito == 1) %>% summarize(matches = n(),
                                                             prop_gts = matches/nrow(nmt_results)*100) # 26.3%
nmt_results %>% dplyr::filter(RF_nmt_minor == 1) %>% summarize(matches = n(),
                                                             prop_gts = matches/nrow(nmt_results)*100) # 7.52%

# How many matches are there in control genes?
cont_results %>% dplyr::filter(RF_cont_spp == 1) %>% summarize(matches = n(),
                                                             prop_gts = matches/nrow(cont_results)*100) # 69.6%
cont_results %>% dplyr::filter(RF_cont_mito == 1) %>% summarize(matches = n(),
                                                              prop_gts = matches/nrow(cont_results)*100) # 28.0%
cont_results %>% dplyr::filter(RF_cont_minor == 1) %>% summarize(matches = n(),
                                                             prop_gts = matches/nrow(cont_results)*100) # 2.40%

# (5) Logistic regression -------------------------------------------------

cont_results_mito <- cont_results %>% 
  mutate(gene_set = "control") %>% 
  rename("RF_mito" = RF_cont_mito,
         "RF_minor" = RF_cont_minor,
         "RF_spp" = RF_cont_spp)
nmt_results_mito <- nmt_results %>% 
  mutate(gene_set = "nmt") %>% 
  rename("RF_mito" = RF_nmt_mito,
         "RF_minor" = RF_nmt_minor,
         "RF_spp" = RF_nmt_spp)

logres_dat <- bind_rows(cont_results_mito, nmt_results_mito)

# Classification is mitotype, gene set is Nmt or control
mod <- glm(RF_mito ~ gene_set, 
           family = binomial(link = "logit"),
           data = logres_dat)

# mod <- glm(RF_spp ~ gene_set,
#            family = binomial(link = "logit"),
#            data = logres_dat)
# 
# mod <- glm(RF_minor ~ gene_set,
#            family = binomial(link = "logit"),
#            data = logres_dat)

summary(mod) # gene_setnmt is not statistically significant

# Odds ratio
exp(coef(mod)) # 0.9642857 for gene_setnmt; i.e., 4% decrease in topology matches in N-mts compared to control genes

