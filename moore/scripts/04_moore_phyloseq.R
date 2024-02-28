## Load phyloseq
install.packages("phyloseq")
library(phyloseq)

## Plot all samples post-filtering

setwd("/Users/uqibyrne/Desktop/hons_manuscript_2023/moore/phyloseq")

### Load files
otumat <- read.csv("otumat/04b_moore_reads_pcrMerge_collapsed_phyloseq.csv")
taxmat <- read.csv("taxonomy/moore_subsetTaxa_echinoBOLD_phyloseq.csv")
sampledata <- read.csv("metadata/02_moore_sample_phyloseq.csv")

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
library(RColorBrewer)

# Classic palette BuPu, with 4 colors
cols <- brewer.pal(11, "Spectral") 

# Add more colors to this palette :
cols <- colorRampPalette(cols)(16)

levels(df$replicate) <- c("Dec_2015.1", "Dec_2015.2", "Jan_2016.1","Jan_2016.2", 
                          "Mar_2016.1", "Mar_2016.2",
                          "Jun_2016.1", "Jun_2016.2", "Aug_2016.1", "Aug_2016.2",
                          "Dec_2016.1","Dec_2016.2","Jan_2017.1", "Jan_2017.2",
                          "Mar_2017.1","Mar_2017.2","Jun_2017.1","Jun_2017.2",
                          "Aug_2017.1","Aug_2017.2", "Dec_2017.1","Dec_2017.2",
                          "Jan_2018.1","Jan_2018.2", "Mar_2018.1", "Mar_2018.2",
                          "Nov_2019.1","Nov_2019.2","Jan_2020.1", "Jan_2020.2")

fig1 <- plot_bar(physeq, x = "tow", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 35000, 5000), limits = c(0, 35000)) +
  theme_classic() +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1$data$tow <- as.character(fig1$data$tow)
fig1$data$tow <- factor(fig1$data$tow, levels = c("208", "209", "286", "287", "548", "549",
                                                  "639", "640", "684", "685", "774", "775",
                                                  "952", "953", "1167", "1168", "1446", "1447",
                                                  "1470", "1471", "1530", "1531", "1563", "1564",
                                                  "1984", "1985", "100613", "100614", "100991", "100992"))



fig1$data$phylum <- factor(fig1$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Echinodermata", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Porifera", "Sipuncula"))
fig1

## Plot Echinoderm reads across biological replicates

(echino <- subset_taxa(physeq, phylum == "Echinodermata"))

### Load R Color Brewer
install.packages("RColorBrewer")
library(RColorBrewer)

#### Create blues palette
my_blu = brewer.pal(9,"Blues")[3:8]

### Stacked barplot
fig1b <- plot_bar(echino, x = "tow", fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Sample", y = "Reads", color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 6000, 1000), limits = c(0, 6000)) +
  theme_classic() +
  scale_fill_manual(values=my_blu) +
  scale_color_manual(values=my_blu) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1b$data$tow <- as.character(fig1b$data$tow)
fig1b$data$date <- factor(fig1b$data$date, levels = c("208", "209", "286", "287", "548", "549",
                                                      "639", "640", "684", "685", "774", "775",
                                                      "952", "953", "1167", "1168", "1446", "1447",
                                                      "1470", "1471", "1530", "1531", "1563", "1564",
                                                      "1984", "1985", "100613", "100614", "100991", "100992"))

fig1b$data$class <- factor(fig1b$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig1b

### Remove 548 & 549

echino1 = subset_samples(echino, tow != c("548", "549"))

fig1c <- plot_bar(echino1, x = "tow", fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Sample", y = "Reads",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 850, 50), limits = c(0, 850)) +
  theme_classic() +
  scale_fill_manual(values=my_blu) +
  scale_color_manual(values=my_blu) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1c$data$tow <- as.character(fig1c$data$tow)
fig1c$data$date <- factor(fig1c$data$date, levels = c("208", "209", "286", "287",
                                                      "639", "640", "684", "685", "774", "775",
                                                      "952", "953", "1167", "1168", "1446", "1447",
                                                      "1470", "1471", "1530", "1531", "1563", "1564",
                                                      "1984", "1985", "100613", "100614", "100991", "100992"))

fig1c$data$class <- factor(fig1c$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig1c