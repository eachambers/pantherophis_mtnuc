library(here)
library(tidyverse)

# Concatenate GWAS results --------------------------------------------------

files <- list.files(here("data/GWAS"), pattern = "RData") %>% 
  stringr::str_subset(., "traits_etc_0.RData", negate = TRUE) # should be 25

dat <-
  1:length(files) %>% 
  lapply(function(x) {
    load(here("data", "GWAS", files[x]))
    gwas <- gwas %>% 
      mutate(abs_beta = abs(beta))
    return(gwas)
  }) %>% 
  dplyr::bind_rows() # 26,988,595

write_tsv(dat, here("outputs/GWAS_results.txt"))


# Detect outliers -----------------------------------------------------------

sig = 0.05
threshold = -log(sig)

# No outliers found when looking at logp.adj column
dat %>% filter(logp.adj >= threshold) %>% nrow()
max(dat$logp.adj) # 0.3748067

# Combine significance and direction of effect into one value
# Based on code starting here: https://github.com/z0on/Multivariate_GWAS/blob/master/RDA_GWAS_dec2022.R#L42
sign = as.numeric(dat$zscore > 0) # returns 1s and 0s for each SNP
sign[sign == 0] = -1 # switches occurrences of 0s to -1s
dat$signed.logp = dat$logp*sign # assign signs to log p-values based on Z-score results
dat$pos.Mb = dat$pos/1e+6 # convert bp to Mb
dat %>% dplyr::filter(signed.logp > -log10(sig)) %>% nrow() # 1,006,910 outliers

# Do p-value adjustment after having combined chroms together
# Convert back to original p-values
dat$pvals = 10^(-dat$logp)
dat$padj_after = p.adjust(dat$pvals, method = "BH")
dat %>% filter(padj_after <= 0.05) %>% nrow() # 0 outliers