library(tidyverse)
library(here)
library(cowplot)
library(fuzzyjoin)
theme_set(theme_cowplot())

## The following code creates Fig. 2b. It also generates summary statistics from the ABBA-BABA analysis.

##    FILES REQUIRED:
##          Results from ABBA-BABA analysis (emoryi_slowinskii_guttatus_localFstats__50_25.txt)
##          Genomic locations for NMT and control genes (Nmt_coords.txt and controls_coords.txt)


# Import ABBA-BABA data ---------------------------------------------------

dat <- read_tsv(here("ABBABABA", "emoryi_slowinskii_guttatus_localFstats__50_25.txt"), col_names = TRUE) # 120,023

# Import NMT and control coords -------------------------------------------

# Import Nmt information
nmts <- read_tsv(here("ABBABABA", "Nmt_coords.txt"), col_names = FALSE) %>% 
  separate(col = X1, sep = ":", into = c("chr", "sites")) %>% 
  separate(col = sites, sep = "-", into = c("start", "end")) %>% 
  rename(gene = X2)
nmts$start <- as.numeric(nmts$start)
nmts$end <- as.numeric(nmts$end) # 167 genes

# Import control gene information
cont <- read_tsv(here("ABBABABA", "controls_coords.txt"), col_names = FALSE) %>% 
  separate(col = X1, sep = ":", into = c("chr", "sites")) %>% 
  separate(col = sites, sep = "-", into = c("start", "end")) %>% 
  rename(gene = X2)
cont$start <- as.numeric(cont$start)
cont$end <- as.numeric(cont$end) # 142 genes

# Get summary statistics --------------------------------------------------

# Get total numbers of SNPs for cont and nmt datasets
cont %>% 
  mutate(length = end-start) %>% 
  summarize(sum(length)) # 10049438

nmts %>% 
  mutate(length = end-start) %>% 
  summarize(sum(length)) # 2348450

# Now, let's get a genome-wide mean estimate for fdM for nmts and cont genes
nmtchr <- unique(nmts$chr) # 74 unique chroms
nmtdat <- dat %>% filter(chr %in% nmtchr) # reduce dataset size
nmtjoin <-
  fuzzyjoin::fuzzy_inner_join(nmts, nmtdat, 
                              by = c("chr" = "chr",
                                     "start" = "windowStart",
                                     "end" = "windowEnd"),
                              match_fun = list(`==`, `>=`, `<=`))
write_tsv(nmtjoin, here("ABBABABA", "Joined_nmts.txt"), col_names = TRUE)
nmtjoin %>% summarize(meanfdM = mean(f_dM),
                      minfdM = min(f_dM),
                      maxfdM = max(f_dM)) # 0.0546 mean

contchr <- unique(cont$chr) # 87 unique chroms
contdat <- dat %>% filter(chr %in% contchr) # reduce dataset size
contjoin <-
  fuzzyjoin::fuzzy_inner_join(cont, contdat, 
                              by = c("chr" = "chr",
                                     "start" = "windowStart",
                                     "end" = "windowEnd"),
                              match_fun = list(`==`, `>=`, `<=`))
write_tsv(contjoin, here("ABBABABA", "Joined_cont.txt"), col_names = TRUE)
contjoin %>% summarize(meanfdM = mean(f_dM),
                      minfdM = min(f_dM),
                      maxfdM = max(f_dM)) # 0.0607 mean


# -------------------------------------------------------------------------
# Build plot --------------------------------------------------------------

# Take every other line so there's no overlap in regions for visualizations
datsamp <- dat[seq(1, nrow(dat), 2), ] # 60012 obs
# Select one representative scaffold for visualizing
chrom = "NW_023010793.1"
# chrom = "NW_023010694.1"

subset <- datsamp %>% 
  dplyr::filter(chr == chrom)

p <-
  subset %>% 
  dplyr::filter(chr == chrom) %>% 
  ggplot(aes(x = windowStart, y = f_dM)) +
  geom_line(color = "grey") +
  ggtitle(chrom)

p +
  geom_rect(data = nmts %>% 
              filter(chr == chrom) %>% 
              dplyr::rename(windowStart = start, windowEnd = end), 
            aes(NULL, NULL, xmin = windowStart, xmax = windowEnd), 
            ymin = min(subset$f_dM), ymax = max(subset$f_dM),
            fill = "#914a32")

p2 <-
  p +
  geom_rug(data = nmts %>% 
             filter(chr == chrom) %>% 
             dplyr::rename(windowStart = start, windowEnd = end),
           aes(x = windowStart),
           sides = "t",
           inherit.aes = FALSE,
           length = unit(0.05, "npc"),
           color = "#914a32")

p2 +
  geom_rug(data = cont %>% 
             filter(chr == chrom) %>% 
             dplyr::rename(windowStart = start, windowEnd = end),
           aes(x = windowStart),
           sides = "t",
           inherit.aes = FALSE,
           length = unit(0.05, "npc"),
           color = "#2b5ba4") +
  xlab("Genomic position") +
  ylab("Mean fdM") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.0546, color = "#914a32", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0.0607, color = "#2b5ba4", linetype = "dashed", alpha = 0.5) # export 10x4


# Build boxplot -----------------------------------------------------------

contjoin <- contjoin %>% mutate(dataset = "control")
nmtjoin <- nmtjoin %>% mutate(dataset = "nmt")
joined <- bind_rows(contjoin, nmtjoin)

joined %>% 
  ggplot(aes(x = dataset, y = f_dM)) +
  geom_boxplot() +
  theme(axis.title.x = element_blank()) +
  ylab("Mean fdM") # export 5x5
