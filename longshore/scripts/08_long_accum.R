library(vegan)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)

### ~~~~~~~~~~ SPECIES ACCUMULATION CURVES ~~~~~~~~~~ ###

### FORMAT DATA ###

counts <- read.csv("06b_long_counts_sampleMerge_phyloseq.csv")
taxa <- read.csv("03_long_subsetTaxa_phyloseq.csv")

tab <- right_join(counts, taxa, by ='otuID')

tab1 <- tab %>%
  tibble::column_to_rownames("otuID")

tab2 <- t(tab1)

tab2 <- tab2[-c(46:51), ]

write.csv(tab2, "07_long_counts_subsetTaxa.csv")

### LOAD FORMATTED DATA ###

dat <- read.csv("07_long_counts_subsetTaxa.csv")
dat1 <- read.csv("../counts/07b_long_counts_subsetTaxa_specAccum.csv")

### SPECIES ACCUMULATION ###
# If you suspect your sampling is incomplete and you want to estimate total 
  # species richness, methods like "Chao" can be particularly useful.
# For large datasets, "exact" may be impractical, whereas "random" or 
  # "rarefaction" are generally more computationally efficient.
# If the focus is on rare species, "Lomolino" might provide interesting insights.
# If the order of sampling is relevant to your study (e.g., temporal changes 
  # in species richness), "collector" could be the method of choice.

spec_accum_exa <- specaccum(dat1, method = "exact") # Methods look exactly the same!
plot(spec_accum_exa)

### GGPLOT2 ###

# Extracting data
exa_data <- data.frame(Sites = spec_accum_exa$sites, Species = spec_accum_exa$richness)

# Plotting with ggplot2
f <- ggplot(exa_data, aes(x = Sites, y = Species)) +
  geom_line(size=0.5) +
  geom_point() +
  xlab("Number of Samples") +
  ylab("Cumulative Number of OTUs") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 375, 25), limits = c(0, 375)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 45, 5), limits = c(0, 45)) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.tag = element_text(size = 16, face = "bold"))
f

### ~~~~~~~~~~ READ COVERAGE CURVES ~~~~~~~~~~ ###

### FORMAT DATA ###

## Read data (plot d)
long <- read_csv("07c_long_reads_subsetTaxa.csv")

# Use dplyr to aggregate and sort
otu_reads <- long %>%
  group_by(ID) %>%
  summarise(TotalReads = sum(`ReadNumber`)) %>% # Ensure this matches your column name
  arrange(TotalReads)

# Add a column for cumulative species count
otu_reads <- otu_reads %>%
  mutate(CumulativeOTUs = cumsum(TotalReads))

## Custom plot data (plot e)
reads <- read.csv("05_long_reads_sampleMerge_phyloseq.csv")
taxa <- read.csv("03_long_subsetTaxa_phyloseq.csv")
all <- right_join(reads, taxa, by ='otuID')
all <- all[ ,-(47:52)]
all <- all %>%
  tibble::column_to_rownames("otuID")
readSums <- colSums(all)

tab <- right_join(counts, taxa, by ='otuID')
tab <- tab[ ,-c(47:52)]
tab <- tab %>%
  tibble::column_to_rownames("otuID")
countSums <- colSums(tab)

df <- data.frame(TotalReads = unlist(readSums), TotalTaxa = unlist(countSums))

### GGPLOT2 ###

## Cumulative read plot
# Plot the curve with ggplot2
d <- ggplot(otu_reads, aes(x = CumulativeOTUs, y = seq_along(ID))) +
  geom_line() +
  geom_point() +
  labs(x = "Cumulative Number of Reads", y = "Cumulative Number of OTUs") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 375, 25), limits = c(0, 375)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 201000, 40000), limits = c(0,208000)) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.tag = element_text(size = 16, face = "bold"))
d

## Manual plot
# Plotting with ggplot2
e <- ggplot(df, aes(x = TotalReads, y = TotalTaxa)) +
  geom_smooth(method = "gam", 
              color="red", size = 0.5, fill = "lightgrey") + 
  geom_point(size=0.5) +
  xlab("Number of Reads (per sample)") +
  ylab("Number of OTUs") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 50, 5), limits = c(0,50)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 30000, 7500), limits = c(0,33000)) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.tag = element_text(size = 16, face = "bold"))
e

### ~~~~~~~~~~ PREP PLOTS FOR MANUSCRIPT ~~~~~~~~~~ ###

lat <- d + e + f + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = "") 

lat # view multi-panel figure

ggsave("long-curves.png", plot = lat, width = 12, height = 4, dpi = 300)


full <- a + b + c+ d + e + f + plot_layout(ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = "") 

full # view multi-panel figure

ggsave("Fig.png", plot = full, width = 12, height = 8, dpi = 300)
