# Start R session

### Setup workspace
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/longshore")

#### Load phyloseq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)

# Create phyloseq object

## Remove singletons
tab1 <- read.csv("otumat/01_long_reads_phyloseq.csv")
tab1[tab1 < 2] <- 0 # remove hits with <2 reads

## Merge echino BLAST hits
tab2 <- read.csv("long_subsetTaxa_phyloseq.csv")
tab3 <- read.csv("leray_data_echino_ncbi_tax.csv")

install.packages("dplyr")
library(dplyr)

tab4 <- left_join(tab1, tab2, by ='uniq')
write.csv(tab4, "taxonomy/03_long_subsetTaxa_phyloseq.csv")

hits <-read.csv("taxonomy/03_long_subsetTaxa_phyloseq.csv")

hits <- read.csv("taxonomy/05_long_subsetTaxa_echinoBOLD_phyloseq.csv")

## Convert OTU column to row names
otumat <- tab1[,-1]
rownames(otumat) <- tab1[,1]

taxmat1 <- hits[,-1]
rownames(taxmat1) <- hits[,1]
taxmat2 <- as.matrix(taxmat1)

OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat2)

sampledata <- read.csv("metadata/01_long_sampledata_phyloseq.csv")
sampledata1 <- sampledata[,-1]
rownames(sampledata1) <- sampledata[,1]

sampledata2 <- sample_data(sampledata1)

physeq = phyloseq(OTU, TAX, sampledata2)

# Stacked barplots

## Remove control samples
subset = subset_samples(physeq, sample.rep != c("Pos", "Neg"))

## Load required packages
library(ggplot2)

install.packages("devtools")
library(devtools)

devtools::install_github("kwstat/pals")
library(pals)

### Create custom colour palette
cols1 <- c("#9E0142", "#C2294A","#F46D43", "#FA9856", "#FDBE6E", "#FEE08B",
           "#F6FBB2", "#E6F598", "#BEE5A0", "#94D4A4", "#66C2A5", "#439BB5",
           "#4075B4", "#5E4FA2")

## Plot all phyla across biological reps
fig1 <- plot_bar(subset, x = "sample.rep", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30000, 5000), limits = c(0, 30000)) +
  theme_classic() +
  scale_fill_manual(values = cols1) +
  scale_color_manual(values = cols1) +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 15, face="bold"))

fig1$data$sample.rep <- factor(fig1$data$sample.rep, levels = c("Fore-Aft.1","Fore-Aft.2",
                                                              "Bramble.1","Bramble.2",
                                                              "Brittomart.1","Brittomart.2",
                                                              "Otter.1","Otter.2",
                                                              "Yamacutta.1","Yamacutta.2",
                                                              "Eddy.1","Eddy.2",
                                                              "Hall-Thompson.1","Hall-Thompson.2",
                                                              "Gibson.1","Gibson.2",
                                                              "Sudbury.1","Sudbury.2",
                                                              "Green.1","Green.2",
                                                              "Pixie.1","Pixie.2",
                                                              "Tongue.1","Tongue.2",
                                                              "Undine.1","Undine.2",
                                                              "Osterlund.1","Osterlund.2",
                                                              "Lizard.1","Lizard.2"))

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

### Save for publication
ggsave(# _plot name_ if you have it, else it will save the last one printed 
  filename = "/Users/uqibyrne/Desktop/hons_manuscript_2023/longshore/figures/S1Fig3.png", 
  dpi = 300, # for publications  
  units = "cm", 
  width = 25, height = 20,
  bg = "white")
