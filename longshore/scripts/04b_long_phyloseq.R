## Load phyloseq
install.packages("phyloseq")
library(phyloseq)

## Load directory
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/longshore/phyloseq")

## Load files
otumat <- read.csv("otumat/06b_long_counts_sampleMerge_phyloseq.csv")
taxmat <- read.csv("taxonomy/03_long_subsetTaxa_phyloseq.csv")
sampledata <- read.csv("02_long_sampledata_phyloseq.csv")

## Convert OTU column to row names
otumat1 <- otumat[,-1]
rownames(otumat1) <- otumat[,1]

taxmat1 <- taxmat[,-1]
rownames(taxmat1) <- taxmat[,1]
taxmat2 <- as.matrix(taxmat1)

sampledata1 <- sampledata[,-1]
rownames(sampledata1) <- sampledata[,1]
sampledata2 <- sample_data(sampledata1)

### Create physeq object
OTU <- otu_table(otumat1, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat2)

physeq = phyloseq(OTU, TAX, sampledata2)


## Plot counts for all taxa across unique samples

### Load ggplot2
library(ggplot2)

### Create custom colour palette
mycolors <- c("black", "dimgray", "lightgrey", "palegreen4", "yellowgreen", 
              "dodgerblue3","deepskyblue3","lightblue","mediumpurple3", "mediumpurple1", 
              "thistle","lightpink4","lightpink3","lightpink2", "lightpink","beige")

### Stacked barplot
fig1 <- plot_bar(physeq, x = "site", fill = "phylum", facet_grid = "year") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Site", y = "Counts", color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 50, 10), limits = c(0, 50)) +
  theme_classic() +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1$data$site <- factor(fig1$data$site, levels = c("Fore & Aft Reef", "Bramble Reef","Brittomart Reef",
                                                              "Otter Reef","Yamacutta Reef","Eddy Reef","Hall-Thompson Reef",
                                                              "Gibson Reef","Sudbury Reef","Green Island","Pixie Reef",
                                                              "Tongue Reef","Undine Reef","Osterlund Reef",
                                                              "Lizard Island"))

fig1$data$phylum <- factor(fig1$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Echinodermata", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Porifera", "Sipuncula"))
fig1


## Plot Echinoderm counts across unique samples

### Subset data to show Echinodermata only
(echino <- subset_taxa(physeq, phylum == "Echinodermata"))

### Plot 2019 only (also remove any sites with zero reads)
(echino19 <- subset_samples(echino, year == "2019"))

### Define a vector of site names to remove
sites_to_remove <- c("Fore & Aft Reef", "Tongue Reef", "Pixie Reef",
                     "Brittomart Reef", "Yamacutta Reef")

### Subset the samples in echino19 that are not from the specified sites
echino19b <- subset_samples(echino19, !site %in% sites_to_remove)

### Create blues palette
library(RColorBrewer)
my_blu = brewer.pal(9,"Blues")[3:8]

### Stacked barplot
fig2 <- plot_bar(echino, x = "site", fill = "class", facet_grid = "year") + 
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Site", y = "Counts", color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 5, 1), limits = c(0, 5)) +
  theme_classic() +
  scale_fill_manual(values=my_blu) +
  scale_color_manual(values=my_blu) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))


fig2$data$site <- factor(fig2$data$site, levels = c("Bramble Reef","Otter Reef","Eddy Reef","Hall-Thompson Reef",
                                                      "Gibson Reef","Sudbury Reef","Green Island",
                                                      "Undine Reef","Osterlund Reef",
                                                      "Lizard Island"))

fig2$data$site <- factor(fig2$data$site, levels = c("Fore & Aft Reef", "Bramble Reef","Brittomart Reef",
                                                    "Otter Reef","Yamacutta Reef","Eddy Reef","Hall-Thompson Reef",
                                                    "Gibson Reef","Sudbury Reef","Green Island","Pixie Reef",
                                                    "Tongue Reef","Undine Reef","Osterlund Reef",
                                                    "Lizard Island"))

fig2$data$class <- factor(fig2$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig2
