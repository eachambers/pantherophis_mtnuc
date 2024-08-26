library(tidyverse)
library(here)
library(cowplot)
library(scales)
theme_set(theme_cowplot())

## The following code creates Figures 2B and S2 using the results from the ABBA-BABA analysis.

##    FILES REQUIRED:
##          Joined_cont.txt and Joined_nmts.txt, produced using `ABBABABA.R`

##    STRUCTURE OF CODE:
##              (1) Import data
##              (2) Build Fig. 2B
##              (3) Build Fig. S3


# (1) Import data ---------------------------------------------------------

contjoin <- read_tsv(here("data", "Joined_cont.txt"), col_names = TRUE)
nmtjoin <- read_tsv(here("data", "Joined_nmts.txt"), col_names = TRUE)

contjoin <- contjoin %>% mutate(dataset = "control")
nmtjoin <- nmtjoin %>% mutate(dataset = "nmt")
joined <- bind_rows(contjoin, nmtjoin)


# (2) Build Fig 2B --------------------------------------------------------

joined %>% 
  ggplot(aes(x = dataset, y = f_dM)) +
  geom_boxplot() +
  theme(axis.title.x = element_blank()) +
  ylab("Mean fdM") # export 5x5


# (3) Build Fig S3 --------------------------------------------------------

dat <- data.frame(topology = c("BBAA", "ABBA", "BABA"),
                  score = c(1.59484e+06, 1.06098e+06, 519655))

dat %>% 
  ggplot(aes(x = topology, y = score)) + 
  geom_bar(stat = "identity") +
  scale_y_continuous(label = comma) +
  ylab("Shared derived alleles") +
  xlab("Topology") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) # export 6.4x6.4
