library(vegan)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("~/Documents/GitHub/echino-metabarcoding/moore/curve")

### ~~~~~~~~~~ SPECIES ACCUMULATION CURVES ~~~~~~~~~~ ###

### FORMAT DATA ###

count <- read.csv("06b_moore_counts_sampleMerge_phyloseq.csv")
tax <- read.csv("moore_subsetTaxa_phyloseq.csv")

merged <- right_join(count, tax, by ='otuID')

merged <- merged %>%
  tibble::column_to_rownames("otuID")

merged <- t(merged)

merged <- merged[-c(16:21), ]

write.csv(merged, "07_moore_counts_subsetTaxa.csv")

### LOAD FORMATTED DATA ###

df <- read.csv("07_moore_counts_subsetTaxa_specAccum.csv") # numeric data only

### SPECIES ACCUMULATION ###
# If you suspect your sampling is incomplete and you want to estimate total 
# species richness, methods like "Chao" can be particularly useful.
# For large datasets, "exact" may be impractical, whereas "random" or 
# "rarefaction" are generally more computationally efficient.
# If the focus is on rare species, "Lomolino" might provide interesting insights.
# If the order of sampling is relevant to your study (e.g., temporal changes 
# in species richness), "collector" could be the method of choice.

spec_accum <- specaccum(df, method = "exact") # Methods look exactly the same!
plot(spec_accum)

### GGPLOT2 ###

# Extracting data (plot c)
dat <- data.frame(Sites = spec_accum$sites, Species = spec_accum$richness)

# Plotting with ggplot2
c <- ggplot(dat, aes(x = Sites, y = Species)) +
  geom_line(size=0.5) +
  geom_point() +
  xlab("Number of Samples") +
  ylab("Cumulative Number of OTUs") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1600, 100), limits = c(0, 1600)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 15, 3), limits = c(0, 15)) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.tag = element_text(size = 16, face = "bold"))
c

### ~~~~~~~~~~ READ COVERAGE CURVES ~~~~~~~~~~ ###

### FORMAT DATA ###

## Read data (plot a)
d1 <- read.csv("05_moore_reads_sampleMerge_phyloseq.csv")
d2 <- read.csv("moore_subsetTaxa_phyloseq.csv")
d3 <- left_join(d2, d1, by = "otuID")
d3 <- d3 %>%
  tibble::column_to_rownames("otuID")
d3 <- d3[,-(1:3)]

d3 <- as.data.frame(d3) # Ensure it's a data frame
d3$OTUs <- rownames(d3) # Convert row names to a column

# Transform the data to long format
d4 <- d3 %>%
  pivot_longer(
    -OTUs, # Keep the OTUs column as is
    names_to = "sampleID", # This will contain your original column names (samples)
    values_to = "readCount" # This will contain the read counts
  )

write.csv(d4, "07c_moore_reads_subsetTaxa.csv") # make sure to remove any zeros

moore <- read_csv("07c_moore_reads_subsetTaxa.csv")

# Use dplyr to aggregate and sort
moore1 <- moore %>%
  group_by(OTUs) %>%
  summarise(TotalReads = sum(`readCount`)) %>% # Ensure this matches your column name
  arrange(TotalReads)

# Add a column for cumulative species count
moore2 <- moore1 %>%
  mutate(CumulativeOTUs = cumsum(TotalReads))

## Custom plot data (plot b)
d3 <- left_join(d2, d1, by = "otuID")
d3 <- d3[,-(2:4)]
reads <- d3 %>%
  tibble::column_to_rownames("otuID")
readSums <- colSums(reads)

merged <- right_join(count, tax, by ='otuID')
merged <- merged[ ,-(17:19)]
counts <- merged %>%
  tibble::column_to_rownames("otuID")
countSums <- colSums(counts)

df <- data.frame(TotalReads = unlist(readSums), TotalTaxa = unlist(countSums))

write.csv(df, "moore_totalReads_totalOTUs.csv")

### GGPLOT2 ###

## Cumulative read plot
# Plot the curve with ggplot2
a <- ggplot(moore2, aes(x = CumulativeOTUs, y = seq_along(OTUs))) +
  geom_line() +
  geom_point() +
  labs(x = "Cumulative Number of Reads", y = "Cumulative Number of OTUs") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1600, 100), limits = c(0, 1600)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 330000, 50000), limits = c(0, 330000)) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.tag = element_text(size = 16, face = "bold"))
a

## Manual plot
# Plotting with ggplot2
b <- ggplot(df, aes(x = TotalReads, y = TotalTaxa)) +
  geom_smooth(method = "gam", 
              color="red", size = 0.5, fill = "lightgrey") + 
  geom_point(size=0.5) +
  xlab("Number of Reads (per sample)") +
  ylab("Number of OTUs (per sample)") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 275, 25), limits = c(0,275)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 47000, 7500), limits = c(0,47000)) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.tag = element_text(size = 16, face = "bold"))
b

### ~~~~~~~~~~ PREP PLOTS FOR MANUSCRIPT ~~~~~~~~~~ ###

mor <- a + b + c + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = "") 

mor # view multi-panel figure

ggsave("moore-curves.png", plot = mor, width = 12, height = 4, dpi = 300)
