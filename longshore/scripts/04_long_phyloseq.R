
## Load phyloseq
install.packages("phyloseq")
library(phyloseq)

## Plot all samples post-filtering

### Load directory
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/longshore/phyloseq")

### Load files
otumat <- read.csv("otumat/04_long_reads_percentThresh_phyloseq.csv")
taxmat <- read.csv("taxonomy/05_long_subsetTaxa_echinoBOLD_phyloseq.csv")
sampledata <- read.csv("metadata/01_long_sampledata_phyloseq.csv")

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

cols1 <- c("#9E0142", "#C2294A","#F46D43", "#FA9856", "#FDBE6E", "#FEE08B",
           "#F6FBB2", "#E6F598", "#BEE5A0", "#94D4A4", "#66C2A5", "#439BB5",
           "#4075B4", "#5E4FA2")


fig1 <- plot_bar(physeq, x = "tow", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 35000, 5000), limits = c(0, 35000)) +
  theme_classic() +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1$data$tow <- as.character(fig1$data$tow)
fig1$data$tow <- factor(fig1$data$tow, levels = c())

fig1$data$phylum <- factor(fig1$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Echinodermata", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Porifera", "Sipuncula"))
fig1

### Remove the positive controls
subset = subset_samples(physeq, tow != c("EC", "MC", "OC"))

fig1b <- plot_bar(subset, x = "tow", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 27000, 3000), limits = c(0, 27000)) +
  theme_classic() +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1b$data$tow <- as.character(fig1b$data$tow)
fig1b$data$tow <- factor(fig1b$data$tow, levels = c("1623","1624","1629","1630","1632","1633","1635","1636","1638","1639","1641","1642","1644","1645","1650","1651","1656","1657",
                                                    "1662","1663","1668","1669","1674","1675","1677","1678","1689","1690","1755","1756",
                                                    "2662","2663","2647","2648","2644","2645","2641","2642","2638","2639","2635","2636","2632","2633","2626","2627",
                                                    "2620","2621","2614","2615","2608","2609","2602","2603","2599","2600","2587","2588","2572","2573",
                                                    "101341","101343","101320","101321","101295","101297","101283","101286","101276","101279",
                                                    "101268","101270","101262","101264","101246","101248","101236","101238","101226","101228",
                                                    "101221","101222","101200","101202","101176","101178","101169","101171","101156","101157"))

fig1b$data$phylum <- factor(fig1b$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Echinodermata", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Porifera", "Sipuncula"))
fig1b


## Plot reads across unique samples

### Remove the neg control samples
subset1 = subset_samples(subset, tow != c("NCA", "NCB", "NCC"))

### Stacked barplot
fig2 <- plot_bar(subset1, x = "sample.rep", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Site", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30000, 5000), limits = c(0, 30000)) +
  theme_classic() +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 14, face="bold"))

fig2$data$phylum <- factor(fig2$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Porifera", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Echinodermata","Sipuncula"))
fig2

S1F3 <- fig2 + facet_wrap(~year, nrow=3, strip.position = "right") +
  theme(panel.spacing = unit(4, "mm"),
        strip.background = element_rect(fill="lightgrey"),
        strip.placement = "outside",
        strip.text = element_text(size=12, face="bold"))

S1F3

#### Save as png for publication
ggsave(# _plot name_ if you have it, else it will save the last one printed 
  filename = "/Users/uqibyrne/Desktop/echinoMS/S1_Fig3.png", 
  dpi = 300, # for publications  
  units = "cm", 
  width = 25, height = 20,
  bg = "white")


## Plot Echinoderm reads across all samples

(echino <- subset_taxa(subset, phylum == "Echinodermata"))

### Load R Color Brewer
install.packages("RColorBrewer")
library(RColorBrewer)

### Create blues palette
my_blu = brewer.pal(9,"Blues")[3:8]

### Stacked barplot
fig3 <- plot_bar(echino, x = "tow", fill = "class") + 
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

fig3$data$tow <- factor(fig3$data$tow, levels = c("1623","1624","1629","1630","1632","1633","1635","1636","1638","1639","1641","1642","1644","1645","1650","1651","1656","1657",
                                                              "1662","1663","1668","1669","1674","1675","1677","1678","1689","1690","1755","1756",
                                                              "2662","2663","2647","2648","2644","2645","2641","2642","2638","2639","2635","2636","2632","2633","2626","2627",
                                                              "2620","2621","2614","2615","2608","2609","2602","2603","2599","2600","2587","2588","2572","2573",
                                                              "101341","101343","101320","101321","101295","101297","101283","101286","101276","101279",
                                                              "101268","101270","101262","101264","101246","101248","101236","101238","101226","101228",
                                                              "101221","101222","101200","101202","101176","101178","101169","101171","101156","101157"))

fig3$data$class <- factor(fig3$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig3


## Plot Echinoderm reads across unique samples

### Stacked barplot
fig4 <- plot_bar(echino, x = "site", fill = "class", facet_grid = "year") + 
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

fig4$data$site.year <- factor(fig4$data$site.year, levels = c("Fore&Aft Reef", "Bramble Reef","Brittomart Reef",
                                                              "Otter Reef","Yamacutta Reef","Eddy Reef","Hall-Thompson Reef",
                                                              "Gibson Reef","Sudbury Reef","Green Island","Pixie Reef",
                                                              "Tongue Reef","Undine Reef","Osterlund Reef",
                                                              "Lizard Island"))

fig4$data$class <- factor(fig4$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig4

### Plot 2019 only (also remove any sites with zero reads)
(echino19 <- subset_samples(echino, year == "2019"))

fig4b$data$site <- factor(fig4b$data$site, levels = c("Fore&Aft Reef", "Bramble Reef","Brittomart Reef",
                                                      "Otter Reef","Yamacutta Reef","Eddy Reef","Hall-Thompson Reef",
                                                      "Gibson Reef","Sudbury Reef","Green Island","Pixie Reef",
                                                      "Tongue Reef","Undine Reef","Osterlund Reef",
                                                      "Lizard Island"))

### Define a vector of site names to remove
sites_to_remove <- c("Fore&Aft Reef", "Tongue Reef", "Pixie Reef",
                     "Brittomart Reef", "Yamacutta Reef")

### Subset the samples in echino19 that are not from the specified sites
echino19b <- subset_samples(echino19, !site %in% sites_to_remove)

fig4c <- plot_bar(echino19b, x = "site", fill = "class") + 
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


fig4c$data$site <- factor(fig4c$data$site, levels = c("Bramble Reef","Otter Reef","Eddy Reef","Hall-Thompson Reef",
                                                      "Gibson Reef","Sudbury Reef","Green Island",
                                                      "Undine Reef","Osterlund Reef",
                                                      "Lizard Island"))

fig4c$data$class <- factor(fig4c$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig4c
