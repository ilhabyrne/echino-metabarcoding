
## Load phyloseq
install.packages("phyloseq")
library(phyloseq)

## Plot all samples post-filtering

### Load directory
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/longshore/phyloseq")

### Load files
otumat <- read.csv("otumat/04_long_reads_percentThresh_phyloseq.csv")
taxmat <- read.csv("taxonomy/03_long_subsetTaxa_phyloseq.csv")
sampledata <- read.csv("01_long_sampledata_phyloseq.csv")

### Convert OTU column to row names
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

### Stacked barplot
library(ggplot2)

mycolors <- c("black", "dimgray", "lightgrey", "palegreen4", "yellowgreen", 
              "dodgerblue3","deepskyblue3","lightblue","mediumpurple3", "mediumpurple1", 
              "thistle","lightpink4","lightpink3","lightpink2", "lightpink","beige")

subset = subset_samples(physeq, tow != c("EC", "MC", "OC"))

fig1 <- plot_bar(subset, x = "site.rep", fill = "phylum", facet_grid = "year") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30000, 5000), limits = c(0, 30000)) +
  theme_classic() +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1$data$site.rep <- factor(fig1$data$site.rep, levels = c("Fore_1","Fore_2","Bram_1","Bram_2",
                                                            "Brit_1","Brit_2","Otte_1","Otte_2",
                                                            "Yama_1","Yama_2","Eddy_1","Eddy_2",
                                                            "Hall_1","Hall_2","Gibs_1","Gibs_2",
                                                            "Sudb_1","Sudb_2","Gree_1","Gree_2",
                                                            "Pixi_1","Pixi_2","Tong_1","Tong_2",
                                                            "Undi_1","Undi_2","Oste_1","Oste_2",
                                                            "Liza_1","Liza_2"))

fig1$data$phylum <- factor(fig1$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Echinodermata", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Porifera", "Sipuncula"))
fig1

## Plot reads across unique samples

### Stacked barplot
fig2 <- plot_bar(subset, x = "site", fill = "phylum", facet_grid = "year") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Site", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30000, 5000), limits = c(0, 30000)) +
  theme_classic() +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig2$data$site.year <- factor(fig2$data$site.year, levels = c("Fore & Aft Reef", "Bramble Reef","Brittomart Reef",
                                                              "Otter Reef","Yamacutta Reef","Eddy Reef","Hall-Thompson Reef",
                                                              "Gibson Reef","Sudbury Reef","Green Island","Pixie Reef",
                                                              "Tongue Reef","Undine Reef","Osterlund Reef",
                                                              "Lizard Island"))

fig2$data$phylum <- factor(fig2$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Echinodermata", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Porifera", "Sipuncula"))
fig2


## Plot Echinoderm reads across all samples

(echino <- subset_taxa(subset, phylum == "Echinodermata"))

### Load R Color Brewer
install.packages("RColorBrewer")
library(RColorBrewer)

### Create blues palette
my_blu = brewer.pal(9,"Blues")[3:8]

### Stacked barplot
fig3 <- plot_bar(echino, x = "site.rep", fill = "class", facet_grid = "year") + 
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Sample", y = "Reads",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 550, 50), limits = c(0, 550)) +
  theme_classic() +
  scale_fill_manual(values=my_blu) +
  scale_color_manual(values=my_blu) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig3$data$site.rep <- factor(fig3$data$site.rep, levels = c("Fore_1","Fore_2","Bram_1","Bram_2",
                                                            "Brit_1","Brit_2","Otte_1","Otte_2",
                                                            "Yama_1","Yama_2","Eddy_1","Eddy_2",
                                                            "Hall_1","Hall_2","Gibs_1","Gibs_2",
                                                            "Sudb_1","Sudb_2","Gree_1","Gree_2",
                                                            "Pixi_1","Pixi_2","Tong_1","Tong_2",
                                                            "Undi_1","Undi_2","Oste_1","Oste_2",
                                                            "Liza_1","Liza_2"))

fig3$data$class <- factor(fig3$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig3

## Plot Echinoderm reads across all unique samples

(echino <- subset_taxa(subset, phylum == "Echinodermata"))

### Load R Color Brewer
install.packages("RColorBrewer")
library(RColorBrewer)

### Create blues palette
my_blu = brewer.pal(9,"Blues")[3:8]

### Stacked barplot
fig4 <- plot_bar(echino, x = "site", fill = "class", facet_grid = "year") + 
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Sample", y = "Reads",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 900, 150), limits = c(0, 900)) +
  theme_classic() +
  scale_fill_manual(values=my_blu) +
  scale_color_manual(values=my_blu) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig4$data$site.rep <- factor(fig4$data$site.rep, levels = c("Fore_1","Fore_2","Bram_1","Bram_2",
                                                            "Brit_1","Brit_2","Otte_1","Otte_2",
                                                            "Yama_1","Yama_2","Eddy_1","Eddy_2",
                                                            "Hall_1","Hall_2","Gibs_1","Gibs_2",
                                                            "Sudb_1","Sudb_2","Gree_1","Gree_2",
                                                            "Pixi_1","Pixi_2","Tong_1","Tong_2",
                                                            "Undi_1","Undi_2","Oste_1","Oste_2",
                                                            "Liza_1","Liza_2"))

fig4$data$class <- factor(fig4$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig4

### Plot 2019 only (also remove any sites with zero reads)
(echino19 <- subset_samples(echino, year == "2019"))

fig3b <- plot_bar(echino19, x = "site.rep", fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Sample", y = "Reads",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 500, 50), limits = c(0, 500)) +
  theme_classic() +
  scale_fill_manual(values=my_blu) +
  scale_color_manual(values=my_blu) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig3b$data$site.rep <- factor(fig3b$data$site.rep, levels = c("Fore_1","Fore_2","Bram_1","Bram_2",
                                                            "Brit_1","Brit_2","Otte_1","Otte_2",
                                                            "Yama_1","Yama_2","Eddy_1","Eddy_2",
                                                            "Hall_1","Hall_2","Gibs_1","Gibs_2",
                                                            "Sudb_1","Sudb_2","Gree_1","Gree_2",
                                                            "Pixi_1","Pixi_2","Tong_1","Tong_2",
                                                            "Undi_1","Undi_2","Oste_1","Oste_2",
                                                            "Liza_1","Liza_2"))

fig3b$data$class <- factor(fig3b$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig3b

### Plot 2019 only (also remove any sites with zero reads)
(echino19 <- subset_samples(echino, year == "2019"))

### Define a vector of site names to remove
sites_to_remove <- c("Fore & Aft Reef", "Tongue Reef", "Pixie Reef",
                     "Brittomart Reef", "Yamacutta Reef")

### Subset the samples in echino19 that are not from the specified sites
echino19b <- subset_samples(echino19, !site %in% sites_to_remove)

fig4b <- plot_bar(echino19b, x = "site", fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Site", y = "Reads",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 900, 100), limits = c(0, 900)) +
  theme_classic() +
  scale_fill_manual(values=my_blu) +
  scale_color_manual(values=my_blu) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))


fig4b$data$site <- factor(fig4b$data$site, levels = c("Bramble Reef","Otter Reef","Eddy Reef","Hall-Thompson Reef",
                                                      "Gibson Reef","Sudbury Reef","Green Island",
                                                      "Undine Reef","Osterlund Reef",
                                                      "Lizard Island"))

fig4b$data$class <- factor(fig4b$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig4b
