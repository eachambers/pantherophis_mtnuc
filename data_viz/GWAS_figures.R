library(here)
library(tidyverse)
library(cowplot)
library(padr)
theme_set(theme_cowplot())

## The following code generates Fig. S5 quantile-quantile plots.

##    FILES REQUIRED:
##          Non-significant SNPs from GWAS analysis ("GWAS_results.txt"), generated using `GWAS_analysis.R`
##          `cont_gathered` and `nmt_gathered` objects from GWAS_analysis.R script; below re-generates them too

##    STRUCTURE OF CODE:
##              (1) Read in NMT and control gene data
##              (2) Process GWAS results
##              (5) Build Q-Q plot (Fig. S5)


# (1) Read in NMT and control gene data -----------------------------------

# If continuing from GWAS_analysis.R script, no need to import the following files
gwas <- read_tsv(here("data", "GWAS_results.txt"), col_names = TRUE)
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


# (2) Process GWAS results ------------------------------------------------

sign = as.numeric(gwas$zscore > 0) # returns 1s and 0s for each SNP
sign[sign == 0] = -1 # switches occurrences of 0s to -1s
gwas$signed.logp = gwas$logp*sign # assign signs to log p-values based on Z-score results
# dat$pos.Mb = dat$pos/1e+6 # convert bp to Mb

gwas$signed.logpadj = gwas$logp.adj*sign # assign signs to log p-values based on Z-score results

# Subset data to retrieve p-values for all SNPs within N-mt genes and control genes
gwas_cont <- left_join(cont_gathered, gwas) %>% 
  na.omit() # 157,943 SNPs with GWAS results
gwas_nmts <- left_join(nmt_gathered, gwas) %>% 
  na.omit() # 37,214 SNPs with GWAS results


# (3) Build qqplot --------------------------------------------------------

# Build Q-Q plot
qqplot(gwas_cont$signed.logp, gwas_nmts$signed.logp)
abline(0,1, col = "red") # export 8x6
pdf(here("outputs/QQplot_all.pdf"), width = 8, height = 6)
dev.off()

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