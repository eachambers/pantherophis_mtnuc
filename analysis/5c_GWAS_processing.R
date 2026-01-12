library(here)
library(tidyverse)
library(fuzzyjoin)
library(padr)
library(data.table)
library(GenomicRanges) # BiocManager::install("GenomicRanges")
# library(GWASTools) # BiocManager::install("GWASTools")
# library(scattermore) # efficient plotting of points for Manhattan plot

## The following code processes results from the GWAS analysis and merges these with the NMT and
## control genes.

##    FILES REQUIRED:
##          Results from the association analysis (25 RData files)
##          Coordinates for the NMT and control genes of interest (NMT_coords.txt & controls_coords.txt)

##    STRUCTURE OF CODE:
##              (1) Process NMT and control data
##              (2) Process genome annotations
##              (3) Process association analysis results
##              (4) Combine association analysis results with genome annotations
##              (5) Run Mann-Whitney U test on absolute betas
##              (6) Run Mann-Whitney U test on N-mt vs control genes
##              (7) Adjusted p-values from association analysis


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
# Replace these 14 occurrences with locus numbers instead of gene names (i.e., retain gff formatting)
cox_regions <- left_join(cox, nmtcont %>% dplyr::rename(gene_name = gene)) %>% na.omit()

# Remove COX genes from `nmtcont` and bind with `cox_regions`
nmtcont <-
  nmtcont %>% 
  filter(!gene %in% cox$gene_name) %>% 
  bind_rows(cox_regions %>% dplyr::select(-gene_name))

# Do a check to ensure that all relevant genes were found
nrow(gff_genes %>% left_join(nmtcont) %>% filter(category == "control")) # 139 genes; CORRECT
nrow(gff_genes %>% left_join(nmtcont) %>% filter(category == "nmt")) # 167 genes; CORRECT

# Export gff_genes for sanity
write_tsv(gff_genes, here("data", "gff_genes.txt"), col_names = TRUE)


# (3) Process GWAS data ---------------------------------------------------

# Iterate through all the RData files (exclude the traits object), read them in, 
# calculate absolute beta and determine which SNPs are within Nmts or control 
# genes, and export these data as a single file
files <- list.files(here("data/GWAS"), pattern = "RData") %>% 
  stringr::str_subset(., "traits_etc_0.RData", negate = TRUE) # should be 25

# files <- list.files(here("data/old_GWAS/12_2023_emslow"), pattern = "RData") %>% 
#   stringr::str_subset(., "traits_etc_0.RData", negate = TRUE) # should be 25

dat <-
  1:length(files) %>% 
  lapply(function(x) {
    load(here("data", "GWAS", files[x]))
    gwas <- gwas %>% 
      mutate(abs_beta = abs(beta))
    return(gwas)
  }) %>% 
  dplyr::bind_rows() # 26,988,595


# (4) Combine GWAS results with genes -------------------------------------

# Convert to GRanges
snps_gr <- GRanges(seqnames = dat$chrom,
                   ranges = IRanges(start = dat$pos, end = dat$pos),
                   abs_beta = dat$abs_beta)
genes_gr <- GRanges(seqnames = gff_genes$chrom,
                    ranges = IRanges(start = gff_genes$start, end = gff_genes$end),
                    gene_id = gff_genes$gene)

# Find overlaps between SNPs and genes; ignore the warning message
hits <- findOverlaps(snps_gr, genes_gr) # 13,976,120

# Convert to df to extract SNP-gene pairs
df <- data.frame(gene_id = mcols(genes_gr)$gene_id[subjectHits(hits)],
                 abs_beta = mcols(snps_gr)$abs_beta[queryHits(hits)]) # 13,976,120

nmt_correct_names <- nmtcont %>% filter(category == "nmt")

# Add flanking regions of +/-2Kb to genes
genes_gr_flank <- resize(genes_gr, 
                         width = width(genes_gr) + 4000, # add 2kb to each side
                         fix = "center")                 # keep gene centered

# Find overlaps between SNPs and genes
hits_flank <- findOverlaps(snps_gr, genes_gr_flank)

# Extract SNP-gene pairs
df_flank <- data.frame(gene_id = mcols(genes_gr_flank)$gene_id[subjectHits(hits_flank)],
                 abs_beta = mcols(snps_gr)$abs_beta[queryHits(hits_flank)])

# Of the 13M SNPs found in genes, how many are within N-mt genes? And within control genes?
df_flank %>% mutate(row_id = row_number()) %>% filter(gene_id %in% nmt_correct_names$gene) %>% nrow() # 51,256
df_flank %>% mutate(row_id = row_number()) %>% filter(gene_id %in% nmt_correct_names$gene) %>% dplyr::select(gene_id) %>% distinct() %>% nrow() # 167
df_flank %>% mutate(row_id = row_number()) %>% filter(gene_id %in% cont$gene) %>% nrow() # 167,036
df_flank %>% mutate(row_id = row_number()) %>% filter(gene_id %in% cont$gene) %>% dplyr::select(gene_id) %>% distinct() %>% nrow() # 139

max_beta_flank <- df_flank %>%
  group_by(gene_id) %>% 
  summarize(max_abs_beta = max(abs_beta)) %>% 
  # Add col with whether it's nmt or not
  mutate(is_gene_nmt = case_when(gene_id %in% nmt_correct_names$gene ~ 1,
                                 .default = 0))

# Export 23638 rows
write_tsv(max_beta_flank, here("data", "abs_max_beta_2kbflank.txt"))


# (5) Mann-Whitney U test: N-mt genes vs all other SNPs ---------------------

# If starting here, retrieve input data:
# max_beta_flank <- read_tsv(here("data", "abs_max_beta_2kbflank.txt"))

# Run Mann-Whitney test; 0 is automatically assigned the first group and 1 is the second
# which means that alternative = "less" is testing whether non-NMTs have significantly lower
# abs_beta values than NMT genes
mwu <- wilcox.test(data = max_beta_flank, max_abs_beta ~ is_gene_nmt, alternative = "less")

# Could manually set the levels on the is_gene_nmt column such that a 1 is first
max_beta_flank$is_gene_nmt <- factor(max_beta_flank$is_gene_nmt, levels = c("1", "0"))
# Check
levels(max_beta_flank$is_gene_nmt)
# Run test again, now specifying 'greater' because NMT will appear first and you want to test
# whether NMT genes have higher abs beta values:
mwu <- wilcox.test(data = max_beta_flank, max_abs_beta ~ is_gene_nmt, alternative = "greater")


# (6) Mann-Whitney U test: N-mt vs. control genes ---------------------------

# Extract N-mt and control genes from the max beta df
relevant_genes <- bind_rows(nmt_correct_names %>% dplyr::select(-category), cont)
max_beta_flank_subset <- max_beta_flank %>% filter(gene_id %in% relevant_genes$gene) # 306 genes total
# Check levels for sanity
levels(max_beta_flank_subset$is_gene_nmt)

# Run MWU test
mwu <- wilcox.test(data = max_beta_flank_subset, max_abs_beta ~ is_gene_nmt, alternative = "greater")


# (7) Adjusted p-values from association analysis ---------------------------

sig = 0.05
threshold = -log(sig)

outliers <- dat %>% filter(logp.adj >= threshold)

sign = as.numeric(dat$zscore > 0) # returns 1s and 0s for each SNP
sign[sign == 0] = -1 # switches occurrences of 0s to -1s
dat$signed.logp = dat$logp*sign # assign signs to log p-values based on Z-score results
dat$pos.Mb = dat$pos/1e+6 # convert bp to Mb

  sig_snps <- dat %>% 
    dplyr::filter(signed.logp > -log10(sig))

  dat$signed.logpadj = dat$logp.adj*sign # assign signs to log p-values based on Z-score results

  sig_snps_adj <- dat %>% 
    dplyr::filter(signed.logpadj > -log10(sig))

# Convert back to original p-values
dat$pvals = 10^(-dat$logp)
dat$pvals_adjusted_after = p.adjust(dat$pvals, method = "fdr")

dat %>% filter(pvals_adjusted_after <= 0.05)

# Overall outliers
outliers <- dat %>% filter(logp.adj >= threshold)
