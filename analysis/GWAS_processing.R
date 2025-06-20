library(here)
library(tidyverse)
library(fuzzyjoin)
library(padr)
library(GWASTools) # BiocManager::install("GWASTools")
library(scattermore) # efficient plotting of points for Manhattan plot

## The following code processes results from the GWAS analysis and merges these with the NMT and
## control genes.

##    FILES REQUIRED:
##          Results from the GWAS analysis (25 RData files)
##          Coordinates for the NMT and control genes of interest (NMT_coords.txt & controls_coords.txt)

##    STRUCTURE OF CODE:
##              (1) Process NMT and control data
##              (2) Process GWAS data
##              (3) Merge GWAS results with NMT and control data, exporting sig SNPs
##              (4) Get summary statistics

# Load relevant functions
source(here("analysis", "GWAS_functions.R"))


# (1) Process NMT and control data ----------------------------------------

# NMT and control files only have ranges so we need to pad data such that
# each position is assigned a row in the dataframe
nmts <- read_tsv(here("data", "Nmt_coords.txt"), col_names = FALSE) %>% 
  rename(seqname = X1, gene = X2) %>% 
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

cont <- read_tsv(here("data", "controls_coords.txt"), col_names = FALSE) %>% 
  rename(seqname = X1, gene = X2) %>% 
  separate(col = seqname, into = c("chrom", "number"), sep = "\\:") %>% # may want to check that col has unique entries for each row `nrow(nmts) == length(unique(nmts$number))` before moving to next line
  separate(col = number, into = c("start", "end"), sep = "-")
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


# (2) Process GWAS data ---------------------------------------------------

# Iterate through all the RData files, read them in, calculate which SNPs are 
# significant, and export these data as a single file
files <- list.files(here("data"), pattern = "RData")
file_no <- 1:length(files)

dat <-
  file_no %>% 
  lapply(function(x) {
    load(here("GWAS", files[x]))
    gwas <- gwas_stats(gwas)
    return(gwas)
  }) %>% 
  dplyr::bind_rows()

# Save GWAS results for all SNPs
write_tsv(dat, here("data", "GWAS_results.txt"), col_names = TRUE)

# Retrieve only significant outliers based on two alpha thresholds
# sig_snps_0.01 <- dat %>%
#   dplyr::filter(signed.logp > -log10(0.01)) # 144,823 at alpha=0.01
# sig_snps_0.01$chrom <- as.character(sig_snps_0.01$chrom)

sig_snps_0.05 <- dat %>% 
  dplyr::filter(signed.logp > -log10(0.05)) # 954,013 at alpha=0.05
sig_snps_0.05$chrom <- as.character(sig_snps_0.05$chrom)

# Get summary statistics (number of outliers per chrom)
# chrom_sites_0.01 <- sig_snps_0.01 %>% group_by(chrom) %>% count() %>% ungroup()
# chrom_sites_0.01 %>% summarize(min_outliers = min(n),
#                                              max_outliers = max(n),
#                                              mean_outliers = mean(n))

chrom_sites_0.05 <- sig_snps_0.05 %>% group_by(chrom) %>% count() %>% ungroup()
chrom_sites_0.05 %>% summarize(min_outliers = min(n),
                                             max_outliers = max(n),
                                             mean_outliers = mean(n))


# (3) Join GWAS results with NMTs and cont genes --------------------------

# Do below for both alpha thresholds
nmtdat <- left_join(sig_snps_0.05, nmt_gathered) %>% 
  dplyr::select(chrom, pos, zscore, pos.Mb, signed.logp, r2, category) %>%
  distinct() %>%
  replace_na(list(category = "non-NMT"))
# write_tsv(nmtdat, here("GWAS", "GWAS_NMT_sigsnps_0.05.txt")) # change file name accordingly

nmtdat %>% 
  filter(category == "NMT") %>% 
  summarize(n()) # 128 SNPs at 0.01 / 1169 SNPs at 0.05

### Same as above but for control genes
contdat <- left_join(sig_snps_0.05, cont_gathered) %>% 
  dplyr::select(chrom, pos, zscore, pos.Mb, signed.logp, r2, category) %>%
  distinct() %>%
  replace_na(list(category = "non-control"))
# write_tsv(contdat, here("GWAS", "GWAS_control_sigsnps_0.05.txt"))

contdat %>% 
  filter(category == "cont") %>% 
  summarize(n()) # 529 SNPs at 0.01 / 5592 SNPs at 0.05


# (4) Get summary statistics ----------------------------------------------

# Which chroms have both Nmt and control genes?
chroms <- inner_join(cont %>% dplyr::select(chrom) %>% distinct(), 
                     nmts %>% dplyr::select(chrom) %>% distinct())


# RUNNING GEMMA -----------------------------------------------------------

# Preparing input data ----------------------------------------------------

# The following is adapted from the script available here: https://github.com/rcc-uchicago/genetic-data-analysis-2/blob/master/format.genotypes.for.gemma.R
library(stringr)

# Load the list of 1,934 samples that were used in Nicod et al, 2016. 
ids <- read.table("listof1934miceusedforanalysis.txt",
                  stringsAsFactors = FALSE)[[1]]

# Repeat for each autosomal chromosome.
for (i in 1:19) {
  cat(sprintf("chromosome %d\n",i))
  
  # Load the genotype data from the .RData file.
  cat(" - Loading genotype data from .RData file.\n")
  input.file <- sprintf("chr%d.prunedgen.final.maf001.0.98.RData",i)
  load(input.file)
  
  # Write the marker positions in a space-delimited text file with one
  # row for each marker, and three columns: (1) marker id, (2)
  # base-pair position, and (3) chromosome. See p. 12 of the GEMMA
  # manual for more information about this file. Here, since the SNP
  # id is not provided, we enter a dummy id of the form snp-X-Y, where
  # X is the chromosome and Y is the base-pair position.
  map.file <- sprintf("chr%02d.map.txt",i)
  cat(" - Writing marker positions to",map.file,"\n")
  pruned_pos$CHR <- i
  pruned_pos     <- transform(pruned_pos,ID = paste("snp",CHR,POS,sep = "-"))
  write.table(pruned_pos[c("ID","POS","CHR")],map.file,sep = " ",
              quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  # Align the columns of the genotype matrix so that they match up with
  # the 1,934 selected ids.
  geno.ids       <- unlist(nameList)
  geno.ids       <- str_replace(geno.ids,fixed("___"),"/")
  geno.ids       <- str_replace(geno.ids,fixed("_recal.reheadered.bam"),"")
  cols           <- match(ids,geno.ids)
  pruned_dosages <- pruned_dosages[,cols]
  
  # Write the mean genotypes ("dosages") as a space-delimited text
  # file in the format used by GEMMA, in which we have one row per
  # marker, and one column per sample (mouse). The first three columns
  # give the marker id, and the two alleles (here we give the
  # alternative allele first, but this isn't required).
  #
  # For some reason the dosages were stored as numbers between 0 and
  # 1, but ordinarily, for diploid organisms, the dosages would be
  # stored as numbers between 0 and 2, representing an expected allele
  # count. It appears that these are the expected counts of the
  # ALTERNATIVE allele (to be certain, we would have to contact the
  # authors to make sure, or see if this information is given
  # somewhere in the paper).
  geno.file <- sprintf("chr%02d.geno.txt",i)
  cat(" - Writing genotype data to",geno.file,"\n")
  pruned_dosages <- 2 * pruned_dosages
  pruned_dosages <- as.data.frame(pruned_dosages,check.names = FALSE)
  pruned_dosages <- round(pruned_dosages,digits = 3)
  pruned_dosages <- cbind(pruned_pos[c("ID","ALT","REF")],pruned_dosages)
  write.table(pruned_dosages,geno.file,sep = " ",quote = FALSE,
              row.names = FALSE,col.names = FALSE)
}


