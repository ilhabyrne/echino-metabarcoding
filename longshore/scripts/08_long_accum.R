library(vegan)
library(tidyverse)
library(ggplot2)

### FORMAT DATA ###

counts <- read.csv("06b_long_counts_sampleMerge_phyloseq.csv")
taxa <- read.csv("03_long_subsetTaxa_phyloseq.csv")

tab <- right_join(counts, taxa, by ='otuID')

tab1 <- tab %>%
  tibble::column_to_rownames("otuID")

tab2 <- t(tab1)

tab2 <- tab2[-c(46:51), ]

write.csv(tab2, "07_long_counts_subsetTaxa.csv")

### LOAD DATA ###

dat <- read.csv("07_long_counts_subsetTaxa.csv")
dat1 <- read.csv("07b_long_counts_subsetTaxa_specAccum.csv")

### RAREFACTION ###
# If you suspect your sampling is incomplete and you want to estimate total 
  # species richness, methods like "Chao" can be particularly useful.
# For large datasets, "exact" may be impractical, whereas "random" or 
  # "rarefaction" are generally more computationally efficient.
# If the focus is on rare species, "Lomolino" might provide interesting insights.
# If the order of sampling is relevant to your study (e.g., temporal changes 
  # in species richness), "collector" could be the method of choice.

spec_accum_rar <- specaccum(dat1, method = "rarefaction")
plot(spec_accum_rar)

spec_accum_exa <- specaccum(dat1, method = "exact")
plot(spec_accum_exa)

### GGPLOT2 ###

# Extracting data
accum_data <- data.frame(Sites = spec_accum$sites, Species = spec_accum$richness)

# Plotting with ggplot2
ggplot(accum_data, aes(x = Sites, y = Species)) +
  geom_line() +
  geom_point() +
  xlab("Number of Samples") +
  ylab("Accumulated Number of OTUs") +
  theme_classic()
