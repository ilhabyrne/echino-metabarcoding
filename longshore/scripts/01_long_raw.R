# Setup workspace
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/longshore")

# Load pyloseq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)

# First glance at the data 
# Remove singletons
tab1 <- read.csv("long_reads_phyloseq.csv")
tab1[tab1 < 2] <- 0 # remove hits with <2 reads

hits <-read.csv("long_allHits_phyloseq.csv")

# Convert OTU column to row names
otumat <- tab1[,-1]
rownames(otumat) <- tab1[,1]

taxmat1 <- hits[,-1]
rownames(taxmat1) <- hits[,1]
taxmat2 <- as.matrix(taxmat1)

OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat2)

sampledata <- read.csv("long_sampledata_phyloseq.csv")
sampledata1 <- sampledata[,-1]
rownames(sampledata1) <- sampledata[,1]

sampledata2 <- sample_data(sampledata1)

physeq = phyloseq(OTU, TAX, sampledata2)


# Stacked barplot
library(ggplot2)

install.packages("devtools")
library(devtools)

devtools::install_github("kwstat/pals")
library(pals)

fig1 <- plot_bar(physeq, x = "site.year") + 
  geom_bar(stat="identity", color="black", fill="black") +
  labs(x = "Sample", y = "Read counts") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 65000, 5000), limits = c(0, 65000)) +
  theme_classic() +
  theme(
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1$data$site.year <- factor(fig1$data$site.year, levels = c("Fore_2017","Bram_2017","Brit_2017","Otte_2017","Yama_2017","Eddy_2017","Hall_2017","Gibs_2017","Sudb_2017","Gree_2017","Pixi_2017","Tong_2017","Undi_2017","Oste_2017","Liza_2017",
                                                        "Fore_2018","Bram_2018","Brit_2018","Otte_2018","Yama_2018","Eddy_2018","Hall_2018","Gibs_2018","Sudb_2018","Gree_2018","Pixi_2018","Tong_2018","Undi_2018","Oste_2018","Liza_2018",
                                                        "Fore_2019","Bram_2019","Brit_2019","Otte_2019","Yama_2019","Eddy_2019","Hall_2019","Gibs_2019","Sudb_2019","Gree_2019","Pixi_2019","Tong_2019","Undi_2019","Oste_2019","Liza_2019",
                                                        "EC","OC","MC", "NCA", "NCB", "NCC"))

fig1$data$phylum <- factor(fig1$data$phylum, levels = c("Eukaryota","non-Eukaryote","Annelida","Arthropoda","Ascomycota",
                                                        "Bacillariophyta","Chaetognatha","Chlorophyta","Chordata","Cnidaria","Ctenophora",
                                                        "Dinophyceae","Echinodermata","Haptista","Hemichordata","Mollusca","Nemertea",
                                                        "Pelagophyceae","Phoronida","Platyhelminthes","Porifera","Sipuncula"))

fig1
