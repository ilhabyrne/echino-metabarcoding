# Set working directory
setwd("~/Desktop/hons_manuscript_2023/moore/variation")


### PLOT BOXPLOTS ###

# Load dataframe
df <- read.csv("tech-reps.csv")

# Load library reshape & format df (only need to do this at the start when first processing the data)
library(reshape)
df_melt <- melt(df)

write.csv(df_melt, "tech-reps.csv")

# Load data
df <- read.csv("tech-reps.csv")

# Load ggplot2
library(ggplot2)

# Set replicate as a factor
df$sample <- as.factor(df$sample)
df$rep <- as.factor(df$rep)

library(dplyr)
levels(df$replicate) <- c("Dec_2015.1", "Dec_2015.2", "Jan_2016.1","Jan_2016.2", 
                                                       "Mar_2016.1", "Mar_2016.2",
                                                       "Jun_2016.1", "Jun_2016.2", "Aug_2016.1", "Aug_2016.2",
                                                       "Dec_2016.1","Dec_2016.2","Jan_2017.1", "Jan_2017.2",
                                                       "Mar_2017.1","Mar_2017.2","Jun_2017.1","Jun_2017.2",
                                                       "Aug_2017.1","Aug_2017.2", "Dec_2017.1","Dec_2017.2",
                                                       "Jan_2018.1","Jan_2018.2", "Mar_2018.1", "Mar_2018.2",
                                                       "Nov_2019.1","Nov_2019.2","Jan_2020.1", "Jan_2020.2")

# Plot boxplot
fig1 <- ggplot(df, aes(x=sample, y=total, fill=rep)) +
  geom_boxplot(alpha=0.5) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 32000, 4000), limits = c(0, 32000)) +
  theme_classic() +
  scale_fill_manual(values=c("white","darkgray")) +
  labs(x="Sample", y = "Read counts", fill = "Replicate") +
  theme(axis.text.x = element_text(colour="black", size=10, angle = 50, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour="black", size=10),
        legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(axis.title = element_text(size = 14, face="bold"))

fig1$data$sample <- factor(fig1$data$sample, levels = c("Dec_2015.1", "Dec_2015.2", "Jan_2016.1","Jan_2016.2", 
                                                        "Mar_2016.1", "Mar_2016.2",
                                                        "Jun_2016.1", "Jun_2016.2", "Aug_2016.1", "Aug_2016.2",
                                                        "Dec_2016.1","Dec_2016.2","Jan_2017.1", "Jan_2017.2",
                                                        "Mar_2017.1","Mar_2017.2","Jun_2017.1","Jun_2017.2",
                                                        "Aug_2017.1","Aug_2017.2", "Dec_2017.1","Dec_2017.2",
                                                        "Jan_2018.1","Jan_2018.2", "Mar_2018.1", "Mar_2018.2",
                                                        "Nov_2019.1","Nov_2019.2","Jan_2020.1", "Jan_2020.2"))
fig1

### PLOT AS PHYSEQ BARGRAPH ###
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

fig2 <- plot_bar(physeq, x = "sample", fill = "phylum") + 
  geom_bar(aes(color = phylum, fill = phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Read counts",color = "Taxa", fill = "Taxa") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 33000, 1500), limits = c(0, 33000)) +
  theme_classic() +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
theme(axis.text.x = element_text(colour="black", size=10, angle = 50, vjust = 1, hjust = 1),
      axis.text.y = element_text(colour="black", size=10),
      legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(axis.title = element_text(size = 14, face="bold"))

fig2$data$sample <- factor(fig2$data$sample, levels = c("Dec_2015.1", "Dec_2015.2", "Jan_2016.1","Jan_2016.2", 
                                                 "Mar_2016.1", "Mar_2016.2",
                                                 "Jun_2016.1", "Jun_2016.2", "Aug_2016.1", "Aug_2016.2",
                                                 "Dec_2016.1","Dec_2016.2","Jan_2017.1", "Jan_2017.2",
                                                 "Mar_2017.1","Mar_2017.2","Jun_2017.1","Jun_2017.2",
                                                 "Aug_2017.1","Aug_2017.2", "Dec_2017.1","Dec_2017.2",
                                                 "Jan_2018.1","Jan_2018.2", "Mar_2018.1", "Mar_2018.2",
                                                 "Nov_2019.1","Nov_2019.2","Jan_2020.1", "Jan_2020.2"))



fig2$data$phylum <- factor(fig2$data$phylum, levels = c("Annelida","Arthropoda", "Brachiopoda",
                                                        "Chaetognatha","Chordata","Cnidaria","Ctenophora",
                                                        "Porifera", "Entoprocta", "Hemichordata","Mollusca","Nemertea",
                                                        "Phoronida","Platyhelminthes","Echinodermata", "Sipuncula"))
fig2

### PANEL FOR MANUSCRIPT ###

install.packages("patchwork")
library(patchwork)

nested <- (fig2/fig1)+
  plot_annotation(tag_levels = 'a') & #add figure labels 
  theme(plot.tag = element_text(face = 'bold', size = 16)) #edit text
nested #view multi-panel figure

ggsave(# _plot name_ if you have it, else it will save the last one printed 
  filename = "/Users/uqibyrne/Desktop/echinoMS/moore/figures/moore_var_supp_20x25_250823.png", 
  dpi = 300, # for publications  
  units = "cm", 
  width = 20, height = 25,
  bg = "white")
