library(here)
library(tidyverse)
library(fuzzyjoin)
library(padr)
# library(GWASTools) # BiocManager::install("GWASTools")
# library(scattermore) # efficient plotting of points for Manhattan plot

## The following code processes results from the GWAS analysis and merges these with the NMT and
## control genes.

##    FILES REQUIRED:
##          Results from the association analysis (25 RData files)
##          Coordinates for the NMT and control genes of interest (NMT_coords.txt & controls_coords.txt)

##    STRUCTURE OF CODE:
##              (1) Process NMT and control data
##              (2) Process association analysis data
##              (3) Merge association analysis results with NMT and control data, exporting sig SNPs
##              (4) Get summary statistics

# Load relevant functions
source(here("analysis", "GWAS_functions.R"))


# (1) Process NMT and control data ----------------------------------------

# NMT and control files only have ranges so we need to pad data such that
# each position is assigned a row in the dataframe
nmts <- read_tsv(here("data", "Nmt_coords.txt"), col_names = c("seqname", "gene")) %>% 
  separate(col = seqname, into = c("chrom", "number"), sep = "\\:") %>% # may want to check that col has unique entries for each row `nrow(nmts) == length(unique(nmts$number))` before moving to next line
  separate(col = number, into = c("start", "end"), sep = "-")
nmts$start <- as.numeric(nmts$start)
nmts$end <- as.numeric(nmts$end)

nmt_gathered <-
  nmts %>% 
  # add unique values for each row
  mutate(gene_group = row_number()) %>% 
  # gather start and end cols
  gather(key = "pos_name", value = "pos",
         start:end) %>% 
  # arrange ascending by row val
  arrange(gene_group, pos) %>% 
  # pad range values
  pad_int(by = "pos", group = "gene_group", step = 1) %>% 
  group_by(gene_group) %>% 
  # fill NA cells with duplicate values; default direction is down
  # need to specify which cols to fill; everything is all cols
  fill(everything()) %>% 
  dplyr::select(-pos_name) %>% 
  mutate(category = "NMT")

cont <- read_tsv(here("data", "control_coords2.txt"), col_names = TRUE)
cont$start <- as.numeric(cont$start)
cont$end <- as.numeric(cont$end)

cont_gathered <-
  cont %>% 
  # add unique values for each row
  mutate(gene_group = row_number()) %>% 
  # gather start and end cols
  gather(key = "pos_name", value = "pos",
         start:end) %>% 
  # arrange ascending by row val
  arrange(gene_group, pos) %>% 
  # pad range values
  pad_int(by = "pos", group = "gene_group", step = 1) %>% 
  group_by(gene_group) %>% 
  # fill NA cells with duplicate values; default direction is down
  # need to specify which cols to fill; everything is all cols
  fill(everything()) %>% 
  dplyr::select(-pos_name) %>% 
  mutate(category = "cont")

gathered <- bind_rows(nmt_gathered, cont_gathered)
nmtcont <- bind_rows(nmts %>% mutate(category = "nmt"), 
                     cont %>% mutate(category = "control"))


# (2) Process genome annotations ------------------------------------------

gff <- read.delim(here("data/GCF_001185365.1_UNIGE_PanGut_3.0_genomic.gff"), 
                  comment.char = '#', 
                  sep = "\t",
                  header = FALSE)
colnames(gff) <- c("chrom", "source", "type", "start", "end", 
                   "score", "strand", "phase", "attributes")

# Retain only gene and extract gene names from attributes col
# There are 25300 genes in the genome
gff_genes <- gff %>% 
  filter(type == "gene") %>% 
  dplyr::select(chrom, start, end, attributes) %>% 
  separate(attributes, into = c("ID", "Dbxref", "Name", "gbkey", "gene", "gene_biotype"), sep = ";") %>% 
  mutate(gene = str_remove_all(gene, "gene=")) %>% 
  dplyr::select(chrom, start, end, gene)

# The COX genes need to be retrieved using LOC numbers
cox <- read_tsv(here("data", "COXgene_LOC.txt"), col_names = c("gene_name", "gene"))
# Replace these 14 occurrences with LOC numbers instead of gene names
cox_regions <- left_join(cox, nmtcont %>% rename(gene_name = gene)) %>% 
  dplyr::select(-gene_name)
# Remove COX genes from `nmtcont` and bind with `cox_regions`
nmtcont <-
  nmtcont %>% 
  filter(!gene %in% cox$gene_name) %>% 
  bind_rows(cox_regions)

gff_genes <- left_join(gff_genes, nmtcont) %>% 
  mutate(snp_in_nmt = case_when(category == "nmt" ~ 1,
                                .default = 0))

# Verify that all were found
nrow(gff_genes %>% filter(category == "control")) # 139 genes; CORRECT
nrow(gff_genes %>% filter(category == "nmt")) # 167 genes; CORRECT


# (3) Process GWAS data ---------------------------------------------------

# Iterate through all the RData files (exclude the traits object), read them in, 
# calculate absolute beta and determine which SNPs are within Nmts or control 
# genes, and export these data as a single file
files <- list.files(here("data/GWAS"), pattern = "RData") %>% 
  stringr::str_subset(., "traits_etc_0.RData", negate = TRUE)

dat <-
  1:length(files) %>% 
  lapply(function(x) {
    load(here("data", "GWAS", files[x]))
    gwas <- gwas %>% 
      mutate(abs_beta = abs(beta)) %>% 
      left_join(gathered) %>% 
      mutate(snp_in_nmt = case_when(category == "NMT" ~ 1,
                                    .default = 0))
    return(gwas)
  }) %>% 
  dplyr::bind_rows()


# (4) Combine GWAS results with genes -------------------------------------

# Make sure column types are consistent for fuzzy joining
dat <- dat %>% 
  mutate(chrom = as.character(chrom),
         pos = as.numeric(pos),
         abs_beta = as.numeric(abs_beta))

gff_genes <- gff_genes %>% 
  mutate(chrom = as.character(chrom),
         start = as.numeric(start),
         end = as.numeric(end))

# Calculate max and mean absolute per-gene beta values
result <- fuzzy_inner_join(dat, gff_genes,
                           by = c("chrom" = "chrom",
                                  "pos" = "start",
                                  "pos" = "end"),
                           match_fun = list(`==`, `>=`, `<=`)) %>%
  group_by(GENE) %>%
  summarize(max_absbeta = max(abs_beta, na.rm = TRUE),
            mean_absbeta = mean(abs_beta, na.rm = TRUE)) %>%
  ungroup()

# Calculate max absolute per-gene beta values with +/-2kb flanking regions
gff_genes_flank <- gff_genes %>% 
  mutate(start_flank = start - 2000,
         end_flank = end + 2000)

result_flank <- fuzzy_inner_join(dat, gff_genes_flank,
                                 by = c("chrom" = "chrom",
                                        "pos" = "start_flank",
                                        "pos" = "end_flank"),
                                 match_fun = list(`==`, `>=`, `<=`)) %>%
  group_by(GENE) %>%
  summarize(max_absbeta = max(abs_beta, na.rm = TRUE),
            mean_absbeta = mean(abs_beta, na.rm = TRUE)) %>%
  ungroup()

# Export results



# (5) Mann-Whitney U test ---------------------------------------------------

# If starting here, retrieve input data:
