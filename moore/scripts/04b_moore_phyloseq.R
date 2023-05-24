## Plot counts across all samples

### Load files
otumat <- read.csv("otumat/06a_moore_counts_pcrMerge_phyloseq.csv")
taxmat <- read.csv("taxonomy/moore_subsetTaxa_phyloseq.csv")
sampledata <- read.csv("02_moore_sample_phyloseq.csv")

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
fig1 <- plot_bar(physeq, x = "tow", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 225, 25), limits = c(0, 225)) +
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

## Plot counts across all unique samples

### Load files
otumat <- read.csv("otumat/06b_moore_counts_sampleMerge_phyloseq.csv")
taxmat <- read.csv("taxonomy/moore_subsetTaxa_phyloseq.csv")
sampledata <- read.csv("03_moore_sample_phyloseq.csv")

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
fig2 <- plot_bar(physeq, x = "date", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Date", y = "Counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 300, 25), limits = c(0, 300)) +
  theme_classic() +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig2$data$date <- factor(fig2$data$date, levels = c("Dec_2015","Jan_2016","Mar_2016","Jun_2016",
                                                    "Aug_2016","Dec_2016","Jan_2017","Mar_2017",
                                                    "Jun_2017","Aug_2017","Dec_2017","Jan_2018",
                                                    "Mar_2018","Nov_2019","Jan_2020"))

fig2$data$phylum <- factor(fig2$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Echinodermata", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Porifera", "Sipuncula"))
fig2


## Plot Echinoderms counts across unique samples

(echino <- subset_taxa(physeq, phylum == "Echinodermata"))

### Load R Color Brewer
install.packages("RColorBrewer")
library(RColorBrewer)

## Create blues palette
my_blu = brewer.pal(9,"Blues")[3:8]

fig3 <- plot_bar(echino, x = "date", fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Date", y = "Counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 16, 1), limits = c(0, 16)) +
  theme_classic() +
  scale_fill_manual(values=my_blu) +
  scale_color_manual(values=my_blu) +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig3$data$date <- factor(fig3$data$date, levels = c("Dec_2015","Jan_2016","Mar_2016","Jun_2016",
                                                    "Aug_2016","Dec_2016","Jan_2017","Mar_2017",
                                                    "Jun_2017","Aug_2017","Dec_2017","Jan_2018",
                                                    "Mar_2018","Nov_2019","Jan_2020"))

fig3$data$class <- factor(fig3$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig3
