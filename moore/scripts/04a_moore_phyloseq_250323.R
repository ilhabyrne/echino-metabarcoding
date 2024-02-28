
## Load phyloseq
install.packages("phyloseq")
library(phyloseq)

## Plot reads across all unique samples

### Load files
otumat <- read.csv("otumat/05_moore_reads_sampleMerge_phyloseq.csv")
taxmat <- read.csv("taxonomy/moore_subsetTaxa_echinoBOLD_phyloseq.csv")
sampledata <- read.csv("metadata/03_moore_sample_phyloseq.csv")

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
## Plot read counts across all unique time points

# Load RColorBrewer
library(RColorBrewer)

# Classic palette BuPu, with 4 colors
cols <- brewer.pal(11, "Spectral") 

# Add more colors to this palette :
cols <- colorRampPalette(cols)(16)

# Plot it
pie(rep(1, length(cols)), col = cols , main="") 

fig2 <- plot_bar(physeq, x = "date", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sampling date", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 47000, 5000), limits = c(0, 47000)) +
  theme_classic() +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 14, face="bold"))

fig2$data$date <- factor(fig2$data$date, levels = c("Dec_2015","Jan_2016","Mar_2016","Jun_2016",
                                                    "Aug_2016","Dec_2016","Jan_2017","Mar_2017",
                                                    "Jun_2017","Aug_2017","Dec_2017","Jan_2018",
                                                    "Mar_2018","Nov_2019","Jan_2020"))

fig2$data$phylum <- factor(fig2$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Porifera", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Echinodermata", "Sipuncula"))
fig2

ggsave(# _plot name_ if you have it, else it will save the last one printed 
  filename = "/Users/uqibyrne/Desktop/hons_manuscript_2023/moore/phyloseq/test.png", 
  dpi = 300, # for publications  
  units = "cm", 
  width = 23, height = 13,
  bg = "white")

## Plot Echinoderm reads across unique samples

(echino <- subset_taxa(physeq, phylum == "Echinodermata"))

### Load R Color Brewer
install.packages("RColorBrewer")
library(RColorBrewer)

### Create blues palette
my_blu = brewer.pal(9,"Blues")[3:8]

### Stacked barplot

fig2b <- plot_bar(echino, x = "date", fill = "class") + 
  geom_bar(aes(color = class, fill = class), stat="identity", position="stack") +
  labs(x = "Sampling date", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2050, 100), limits = c(0, 2050)) +
  theme_classic() +
  scale_fill_manual(values=my_blu) +
  scale_color_manual(values=my_blu) +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 14, face="bold"))

fig2b$data$date <- factor(fig2b$data$date, levels = c("Dec_2015","Jan_2016","Mar_2016","Jun_2016",
                                                      "Aug_2016","Dec_2016","Jan_2017","Mar_2017",
                                                      "Jun_2017","Aug_2017","Dec_2017","Jan_2018",
                                                      "Mar_2018","Nov_2019","Jan_2020"))

fig2b$data$class <- factor(fig2b$data$class, levels = c("Asteroidea","Crinoidea","Echinoidea","Holothuroidea","Ophiuroidea","Echinodermata"))

fig2b

# Create panel of figs for supp

install.packages("patchwork")
library(patchwork)

nested <- (fig2/fig2b)+
  plot_annotation(tag_levels = 'a') & #add figure labels 
theme(plot.tag = element_text(face = 'bold', size = 16)) #edit text
nested #view multi-panel figure

ggsave(# _plot name_ if you have it, else it will save the last one printed 
  filename = "/Users/uqibyrne/Desktop/echinoMS/moore/figures/moore_reads_supp_20x25_250823.png", 
  dpi = 300, # for publications  
  units = "cm", 
  width = 20, height = 25,
  bg = "white")
