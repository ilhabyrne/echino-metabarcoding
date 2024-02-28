## Load phyloseq
install.packages("phyloseq")
library(phyloseq)

## Load directory
setwd("~/Desktop/Undergraduate/hons_manuscript_2023/longshore/phyloseq")

## Load files
otumat <- read.csv("otumat/06b_long_counts_sampleMerge_phyloseq.csv")
taxmat <- read.csv("taxonomy/05_long_subsetTaxa_echinoBOLD_phyloseq.csv")
sampledata <- read.csv("metadata/02_long_sampledata_phyloseq.csv")

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
cols1 <- c("#9E0142", "#C2294A","#F46D43", "#FA9856", "#FDBE6E", "#FEE08B",
           "#F6FBB2", "#E6F598", "#BEE5A0", "#94D4A4", "#66C2A5", "#439BB5",
           "#4075B4", "#5E4FA2")

### Stacked barplot
fig1 <- plot_bar(physeq, x = "site.short", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Site", y = "OTU counts", color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 45, 5), limits = c(0, 45)) +
  theme_classic() +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 14, face="bold"))

fig1$data$site.short <- factor(fig1$data$site.short, levels = c("Fore-Aft", "Bramble","Brittomart",
                                                     "Otter","Yamacutta","Eddy","Hall-Thompson",
                                                     "Gibson","Sudbury","Green","Pixie",
                                                     "Tongue","Undine","Osterlund",
                                                     "Lizard"))

fig1$data$phylum <- factor(fig1$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Porifera", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Echinodermata","Sipuncula"))


fig1 <- fig1 + facet_wrap(~year, nrow=3, strip.position = "right") +
  theme(panel.spacing = unit(4, "mm"),
        strip.background = element_rect(fill="lightgrey"),
        strip.placement = "outside",
        strip.text = element_text(size=12, face="bold"))
fig1


## Plot Echinoderm counts across unique samples

### Subset data to show Echinodermata only
(echino <- subset_taxa(physeq, phylum == "Echinodermata"))

### Create blues palette
library(RColorBrewer)
my_blu = brewer.pal(9,"Blues")[3:8]
my_blu1 = c("#C6DBEF","#4292C6", "#2171B5") #ast, #hol, #oph

### Filter out zero values so they are not plotted

library(tidyverse)

filt <- prune_taxa(taxa_sums(echino) > 0, echino)

otu_df <- as.data.frame(otu_table(filt))
sample_df <- as.data.frame(sample_data(filt))
taxa_df <- as.data.frame(tax_table(filt))

otu_df$uniq <- rownames(otu_df)
taxa_df$uniq <- rownames(taxa_df)
sample_df$sample <- rownames(sample_df)

merged <- otu_df %>%
  pivot_longer(cols = -c(uniq), names_to = "sample", values_to = "counts") %>%
  left_join(taxa_df, by="uniq") %>%
  left_join(sample_df, by="sample")

merged_filt <- merged %>%
  filter(counts > 0)

#### Prepare the data for the stacked barplot
merged_filt$uniq <-NULL

merged_filt1 <- merged_filt %>%
  group_by(class, sample, site.short, year) %>%
  summarize(total_counts=sum(counts))

merged_filt1$total_counts <- as.numeric(merged_filt1$total_counts)
merged_filt1$year <- as.character(merged_filt1$year)

merged_filt1$site.short <- factor(merged_filt1$site.short, levels = c("Fore-Aft", "Bramble","Brittomart",
                                                                      "Otter","Yamacutta","Eddy","Hall-Thompson",
                                                                      "Gibson","Sudbury","Green","Pixie",
                                                                      "Tongue","Undine","Osterlund",
                                                                      "Lizard"))
### Stacked barplot
fig2 <- ggplot(merged_filt1, aes(x=site.short, y=total_counts,fill=class)) +
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Site", y = "OTU counts", color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 5, 1), limits = c(0, 5)) +
  theme_classic() +
  scale_fill_manual(values=my_blu1) +
  scale_color_manual(values=my_blu1) +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 14, face="bold")) 

fig2 <- fig2 + facet_wrap(~year, nrow=3, strip.position = "right") +
  theme(panel.spacing = unit(4, "mm"),
        strip.background = element_rect(fill="lightgrey"),
        strip.placement = "outside",
        strip.text = element_text(size=12, face="bold"))
fig2


fig2$data$site <- factor(fig2$data$site, levels = c("Fore-Aft", "Bramble","Brittomart",
                                                    "Otter","Yamacutta","Eddy","Hall-Thompson",
                                                    "Gibson","Sudbury","Green","Pixie",
                                                    "Tongue","Undine","Osterlund",
                                                    "Lizard"))

fig2$data$class <- factor(fig2$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig2

### Merge figures for publication

#### Load patchwork
library(patchwork)

#### Nest figures
nested <- (fig1|fig2)+
  plot_annotation(tag_levels = 'a') & #add figure labels 
  theme(plot.tag = element_text(face = 'bold', size = 16)) #edit text
nested #view multi-panel figure

#### Save as high-res png
ggsave(# _plot name_ if you have it, else it will save the last one printed 
  filename = "/Users/uqibyrne/Desktop/echinoMS/longshore/figures/long_counts_main_25x20_250823.png", 
  dpi = 300, # for publications  
  units = "cm", 
  width = 25, height = 20,
  bg = "white")

### Tips if you want to subset the data further
#### Plot 2019 only (also remove any sites with zero reads)
(echino19 <- subset_samples(echino, year == "2019"))

#### Define a vector of site names to remove
sites_to_remove <- c("Fore & Aft Reef", "Tongue Reef", "Pixie Reef",
                     "Brittomart Reef", "Yamacutta Reef")

#### Subset the samples in echino19 that are not from the specified sites
echino19b <- subset_samples(echino19, !site %in% sites_to_remove)
