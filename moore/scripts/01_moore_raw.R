# Setup workspace
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/moore")

# Load pyloseq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)

# First glance at the data 
# Remove singletons
tab1 <- read.csv("moore_reads_phyloseq.csv")
tab1[tab1 < 2] <- 0 # remove hits with <2 reads

hits <-read.csv("moore_allHits_phyloseq.csv")

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


# Stacked barplot for each tow
library(ggplot2)

fig1 <- plot_bar(physeq, x = "tow") + 
  geom_bar(stat="identity", color="black", fill="black") +
  labs(x = "Sample", y = "Read counts") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 85000, 5000), limits = c(0, 85000)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 16))

fig1$data$tow <- factor(fig1$data$tow, levels = c("208", "209", "286", "287", "548", "549",
                                                  "639", "640", "684", "685", "774", "775",
                                                  "952", "953", "1167", "1168", "1446", "1447",
                                                  "1470", "1471", "1530", "1531", "1563", "1564",
                                                  "1984", "1985", "100613", "100614", "100991", "100992",
                                                  "EC","MC", "NCA", "NCB", "NCC", "NCD"))

fig1

# Stacked barplot for each tow
fig1a <- plot_bar(physeq, x = "tow", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Read counts", color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 85000, 5000), limits = c(0, 85000)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 16))

fig1a$data$tow <- factor(fig1a$data$tow, levels = c("208", "209", "286", "287", "548", "549",
                                                  "639", "640", "684", "685", "774", "775",
                                                  "952", "953", "1167", "1168", "1446", "1447",
                                                  "1470", "1471", "1530", "1531", "1563", "1564",
                                                  "1984", "1985", "100613", "100614", "100991", "100992",
                                                  "EC","MC", "NCA", "NCB", "NCC", "NCD"))

fig1a$data$phylum <- factor(fig1a$data$phylum, levels = c("Eukaryota","non-Eukaryote","Annelida","Arthropoda","Ascomycota",
                                                        "Bacillariophyta","Chaetognatha","Chlorophyta","Chordata","Cnidaria","Ctenophora",
                                                      "Echinodermata", "Entoprocta", "Haptista","Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Porifera", "Rhodophyta", "Sipuncula"))

fig1a


# Stacked barplot for each date
fig2 <- plot_bar(physeq, x = "date") + 
  geom_bar(stat="identity", color="black", fill="black") +
  labs(x = "Sample", y = "Read counts") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 150000, 25000), limits = c(0, 150000)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 16))

fig2$data$date <- factor(fig2$data$date, level = c("Dec_2015", "Jan_2016", "Mar_2016", "Jun_2016", "Aug_2016", 
                                                   "Dec_2016", "Jan_2017", "Mar_2017", "Jun_2017", "Aug_2017", 
                                                   "Dec_2017", "Jan_2018", "Mar_2018", "Nov_2019", "Jan_2020",
                                                   "EC","MC", "NCA", "NCB", "NCC"))

fig2










