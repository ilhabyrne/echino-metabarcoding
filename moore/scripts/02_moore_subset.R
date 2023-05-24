# Setup workspace
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/moore")

# Load pyloseq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)

# First glance at the data 

# Merge echino BLAST hits
tab1 <- read.csv("moore_subsetTaxa_phyloseq.csv")
tab2 <- read.csv("moore_echino_ncbi.csv")

install.packages("dplyr")
library(dplyr)

tab3 <- left_join(tab1, tab2, by ='uniq')
write.csv(tab3, "moore_subsetTaxa_phyloseq.csv")

# Remove singletons
tab1 <- read.csv("moore_reads_phyloseq.csv")
tab1[tab1 < 2] <- 0 # remove hits with <2 reads

hits <-read.csv("moore_subsetTaxa_phyloseq.csv")

# Convert OTU column to row names
otumat <- tab1[,-1]
rownames(otumat) <- tab1[,1]

taxmat1 <- hits[,-1]
rownames(taxmat1) <- hits[,1]
taxmat2 <- as.matrix(taxmat1)

OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat2)

sampledata <- read.csv("moore_sampledata_phyloseq.csv")
sampledata1 <- sampledata[,-1]
rownames(sampledata1) <- sampledata[,1]

sampledata2 <- sample_data(sampledata1)

physeq = phyloseq(OTU, TAX, sampledata2)

# Stacked barplot
library(ggplot2)

mycolors <- c("black", "dimgray", "lightgrey", "palegreen4", "yellowgreen", 
              "dodgerblue3","deepskyblue3","lightblue","mediumpurple3", "mediumpurple1", 
              "thistle","lightpink4","lightpink3","lightpink2", "lightpink","beige")

fig1 <- plot_bar(physeq, x = "tow", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 70000, 10000), limits = c(0, 70000)) +
  theme_classic() +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1$data$tow <- factor(fig1$data$tow, levels = c("208", "209", "286", "287", "548", "549",
                                                    "639", "640", "684", "685", "774", "775",
                                                    "952", "953", "1167", "1168", "1446", "1447",
                                                    "1470", "1471", "1530", "1531", "1563", "1564",
                                                    "1984", "1985", "100613", "100614", "100991", "100992",
                                                    "EC","MC", "NCA", "NCB", "NCC", "NCD"))

fig1$data$phylum <- factor(fig1$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                          "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                          "Echinodermata", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                          "Phoronida","Platyhelminthes","Porifera", "Sipuncula"))
fig1
