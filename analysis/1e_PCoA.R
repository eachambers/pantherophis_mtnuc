library(vegan)
library(tidyverse)
library(cowplot)
library(here)
theme_set(theme_cowplot())

## The following code generates Figure S2, the PCoA figure, color-coded based on NGSadmix K=4 population assignment.

##    FILES REQUIRED:
##          bams which has the correct ordering of individuals
##          inds2pops1 2-column tab-delimited table of individual assignments to populations; 
##            must be in the same order as samples in the bam list or vcf file
##          panther.ibsMat


# Import data -------------------------------------------------------------

bams = read.table(here("data/bams"), col.names = "filename")
bams <- 
  bams %>% 
  separate(filename, into = c("sampleID", "tmp1", "tmp2", "tmp3"), sep = "_") %>% 
  dplyr::select(sampleID)

# loading individual to population correspondences
i2p = read_tsv(here("data/inds2pops.txt"))
# Link pop assignments to sample ordering from bams; this will ensure
# ordering is consistent with the ibs matrix (i.e., bams file)
samps <- left_join(bams, i2p)

ma = as.matrix(read.table(here("data/panther.ibsMat"))) # actual data


# Performing PCoA ---------------------------------------------------------

# Constrained PCoA using population structure (distance-based RDA)
conds = data.frame(pop = samps$pop)
pp = vegan::capscale(ma ~ conds$pop)


# Plot the results with ggplot ------------------------------------------------

smry <- summary(pp)
df1 <- data.frame(smry$sites[,1:2]) # PC1 and PC2
# df2 <- data.frame(smry$species[,1:2]) # loadings for PC1 and PC2
df <- bind_cols(df1, samps)

colors = c("emoryi" = "#be9739", 
           "guttatus" = "#c7544a",
           "slowinskii" = "#75a3dd",
           "meahllmorum" = "#5d8252")

# parentals <- c("DBS 789", "JJB 6167", "JJB 7567", "TJH 3395")

pcoa_p <-
  ggplot(df, aes(x = MDS1, y = MDS2)) +
  # ggplot(df, aes(x = CAP1, y = CAP2)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgrey") +
  geom_point(aes(color = pop), size = 4) +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  # stat_ellipse(aes(x = MDS1, y = MDS2, fill = pop),
  #          geom = "polygon", alpha = 0.4) +
  # stat_ellipse(aes(x = CAP1, y = CAP2, fill = pop),
  #          geom = "polygon", alpha = 0.4) +
  scale_color_manual(values = colors) +
  # scale_fill_manual(values = colors) +
  geom_label(aes(color = pop), label = df$field_no, nudge_y = 0.1,
           size = 3) +
  ggtitle(paste0(smry$call))


# Export PDF --------------------------------------------------------------

pcoa_p
ggsave(here("outputs/PCoA_plot.pdf"), width = 8, height = 6)


# RUNNING PCOA WITHOUT GUTTATUS -------------------------------------------

# Having a divergent population can skew PCoA results; the following runs
# the PCoA without guttatus individuals

# Import data -------------------------------------------------------------

# There are nine guttatus samples; the following doesn't have these samples included
bams_nogutt = read.table("bams_nogutt")[,1] # should be 124 levels
goods_nogutt = c(1:length(bams_nogutt)) # should be 1 through 124

# loading individual to population correspondences
i2p_nogutt = read.table(here("data/inds2pops1_nogutt"), sep = "\t") # again, ensure 124 rows
row.names(i2p_nogutt) = i2p_nogutt[,1]
i2p_nogutt = i2p_nogutt[goods_nogutt,] 
site_nogutt = i2p_nogutt[,2] # this is just the pop assignment

ma_nogutt = as.matrix(read.table(here("data/myresult.ibsMat_nogutt"))) # V124 final


# Performing PCoA ---------------------------------------------------------

# Number based on the population assignment (1-4)
conds_nogutt = data.frame(cbind(site_nogutt))

pp0_nogutt = capscale(ma_nogutt ~ 1) # capscale is to conduct dbRDA
pp_nogutt = capscale(ma_nogutt ~ site_nogutt, conds_nogutt)
cmd_nogutt = pp0_nogutt # why?


# Plot the results with ggplot ------------------------------------------------

smry_nogutt <- summary(cmd_nogutt)
df1_nogutt <- data.frame(smry_nogutt$sites[,1:2]) # PC1 and PC2
# df2 <- data.frame(smry$species[,1:2]) # loadings for PC1 and PC2

# Merge population assignment with samples
new_df_nogutt <-
  df1_nogutt %>% 
  mutate(site_nogutt)

colors_nogutt = c("SLOW" = "#c7544a",
           "EMSO" = "#75a3dd",
           "EMNO" = "#5d8252",
           "HYBR" = "#a6a4a1")

pcoa_p_nogutt <-
  ggplot(df1_nogutt, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = site_nogutt), size = 4) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme(legend.position = "none",
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 25)) +
  stat_ellipse(aes(x = MDS1, y = MDS2, fill = site_nogutt),
               geom="polygon", alpha=0.4) +
  # coord_fixed() +
  scale_color_manual(values = colors_nogutt) +
  scale_fill_manual(values = colors_nogutt)


# Export PDF --------------------------------------------------------------

pcoa_p_nogutt
ggsave(here("outputs/PCoA_plot_nogutt.pdf"), width = 9, height = 6.3)


# RUNNING PCOA WITHOUT GUTTATUS or MEAHLLMORUM -------------------------------------------

# Having a divergent population can skew PCoA results; the following runs
# the PCoA without guttatus or meahllmorum individuals

# Import data -------------------------------------------------------------

# There is one guttatus and one meahllmorum sample
gutt_meah <- i2p %>% filter(pop == "guttatus" | pop == "meahllmorum") %>% pull(sampleID)
samps_sub <- samps %>% filter(!sampleID %in% gutt_meah)

# Remove guttatus and meahllmorum from IBS matrix
rownames(ma) <- samps$sampleID
colnames(ma) <- samps$sampleID
ma_sub <- ma[,!colnames(ma) %in% gutt_meah]
ma_sub <- ma_sub[!rownames(ma_sub) %in% gutt_meah, ]


# Performing PCoA ---------------------------------------------------------

# Constrained PCoA using population structure (distance-based RDA)
conds = data.frame(pop = samps_sub$pop)
pp = vegan::capscale(ma_sub ~ conds$pop)


# Plot the results with ggplot ------------------------------------------------

smry <- summary(pp)
df1 <- data.frame(smry$sites[,2:3]) # MDS1 and MDS2
df <- bind_cols(df1, samps_sub)

colors_sub = c("emoryi" = "#be9739", "slowinskii" = "#75a3dd")

ggplot(df, aes(x = MDS1, y = MDS2)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgrey") +
  geom_point(aes(color = pop), size = 4) +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  scale_color_manual(values = colors_sub) +
  geom_label(aes(color = pop), label = df$field_no, nudge_y = 0.1,
             size = 3) +
  ggtitle(paste0(smry$call))


# Export PDF --------------------------------------------------------------

ggsave(here("outputs/PCoA_plot_noguttmeah.pdf"), width = 9, height = 6.3)
