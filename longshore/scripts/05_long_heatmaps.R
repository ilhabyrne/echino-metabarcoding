## Heatmaps of echino species at each site

### Load packages required
install.packages("reshape")
library(reshape)

### Merge taxonomy and otumat
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/longshore")
tax <- read.csv("long_echino_tax.csv")
otu <- read.csv("06b_long_counts_sampleMerge_phyloseq.csv")

otutax <- merge(tax, otu, by.x = "otuID", by.y = "otuID")

### Save new dataframe as csv
write_csv(otutax, "long_counts_echino_otutax.csv")

### Load species data
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/longshore/heatmap")
df <- read.csv("long_counts_echinoSp_otutax_23-08-23.csv")
df <- as.data.frame(df) 
head(df)

### Format data frame for plot
df_melt <- melt(df)
head(df_melt)

df_melt$variable <- factor(df_melt$variable, 
                           level = c("Fore_2017", "Bram_2017","Brit_2017",
                                     "Otte_2017","Yama_2017","Eddy_2017","Hall_2017",
                                     "Gibs_2017","Sudb_2017","Gree_2017","Pixi_2017",
                                     "Tong_2017","Undi_2017","Oste_2017",
                                     "Liza_2017",
                                     "Fore_2018", "Bram_2018","Brit_2018",
                                     "Otte_2018","Yama_2018","Eddy_2018","Hall_2018",
                                     "Gibs_2018","Sudb_2018","Gree_2018","Pixi_2018",
                                     "Tong_2018","Undi_2018","Oste_2018",
                                     "Liza_2018",
                                     "Fore_2019", "Bram_2019","Brit_2019",
                                     "Otte_2019","Yama_2019","Eddy_2019","Hall_2019",
                                     "Gibs_2019","Sudb_2019","Gree_2019","Pixi_2019",
                                     "Tong_2019","Undi_2019","Oste_2019",
                                     "Liza_2019"))


library(ggplot2)

df_melt$value <- as.factor(df_melt$value) 

fig1 <- ggplot(df_melt, aes(variable, species, fill=value)) +
  geom_tile(color="grey", show.legend = FALSE) + 
  theme_classic() +
  labs(x="Site & year", y="Species") +
  scale_fill_manual(values = c("white","#4075B4","#4075B4","#4075B4")) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10, face="italic"),
        axis.title = element_text(size=14, face="bold")) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1)) +
  annotate("rect", xmin=0.5, xmax=45.4, ymin=0.512, 
           ymax=1.5, colour="black", fill="transparent", size=1) 
fig1

fig1 <- fig1 + theme(axis.ticks = element_blank(),
                     axis.line = element_blank())
fig1

### Load genus data
df1 <- read.csv("long_counts_echinoGenus_otutax_23-08-23.csv")
df1 <- as.data.frame(df1) 
head(df1)

library(plyr)
df1 <- ddply(df1, "genus", numcolwise(sum)) 
head(df1)

### Format data frame for plot
df1_melt <- melt(df1)
head(df1_melt)

df1_melt$variable <- factor(df1_melt$variable, 
                            level = c("Fore_2017", "Bram_2017","Brit_2017",
                                      "Otte_2017","Yama_2017","Eddy_2017","Hall_2017",
                                      "Gibs_2017","Sudb_2017","Gree_2017","Pixi_2017",
                                      "Tong_2017","Undi_2017","Oste_2017",
                                      "Liza_2017",
                                      "Fore_2018", "Bram_2018","Brit_2018",
                                      "Otte_2018","Yama_2018","Eddy_2018","Hall_2018",
                                      "Gibs_2018","Sudb_2018","Gree_2018","Pixi_2018",
                                      "Tong_2018","Undi_2018","Oste_2018",
                                      "Liza_2018",
                                      "Fore_2019", "Bram_2019","Brit_2019",
                                      "Otte_2019","Yama_2019","Eddy_2019","Hall_2019",
                                      "Gibs_2019","Sudb_2019","Gree_2019","Pixi_2019",
                                      "Tong_2019","Undi_2019","Oste_2019",
                                      "Liza_2019"))

df1_melt$value <- as.factor(df1_melt$value) 

fig2 <- ggplot(df1_melt, aes(variable, genus, fill=value)) +
  geom_tile(color="grey",show.legend = FALSE) + 
  labs(x="Site & year", y="Genera") +
  theme_classic() +
  scale_fill_manual(values=c("white","#4075B4","#4075B4","#4075B4","#4075B4",
                             "#4075B4","#4075B4")) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10, face="italic"),
        axis.title = element_text(size=14, face="bold")) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1)) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1)) +
  theme(axis.title = element_text(colour="black", size=14)) +
  annotate("rect", xmin=0.5, xmax=45.4, ymin=0.515, 
           ymax=1.5, colour="black", fill="transparent", size=1)
fig2

fig2 <- fig2 + theme(axis.ticks = element_blank(),
                     axis.line = element_blank())
fig2


### Combined plot

library(patchwork)

#### Nest panels together
nested <- (fig1/fig2)+
  plot_annotation(tag_levels = 'a') & # add figure labels 
  theme(plot.tag = element_text(face = 'bold', size = 16)) # edit text
nested # view multi-panel figure

#### Save at high-res for publication
ggsave(# _plot name_ if you have it, else it will save the last one printed 
  filename = "long_tile_main_25x20.png", 
  dpi = 300, # for publications  
  units = "cm", 
  width = 25, height = 20,
  bg = "white")


### Just plot 2019
echino19 <- subset(df1, select = c("genus", "Fore_2019", "Bram_2019","Brit_2019",
                                                "Otte_2019","Yama_2019","Eddy_2019","Hall_2019",
                                                "Gibs_2019","Sudb_2019","Gree_2019","Pixi_2019",
                                                "Tong_2019","Undi_2019","Oste_2019",
                                                "Liza_2019"))
df2_melt <- melt(echino19)
head(df2_melt)

df2_melt$variable <- factor(df2_melt$variable, 
                            level = c(
                                      "Fore_2019", "Bram_2019","Brit_2019",
                                      "Otte_2019","Yama_2019","Eddy_2019","Hall_2019",
                                      "Gibs_2019","Sudb_2019","Gree_2019","Pixi_2019",
                                      "Tong_2019","Undi_2019","Oste_2019",
                                      "Liza_2019"))

fig3 <- ggplot(df2_melt, aes(variable, genus, fill=value)) +
  geom_tile(color="grey") + 
  labs(x="Site", y="Genera", fill="# OTUs") +
  scale_fill_gradient(low="white", high="blue") +
  theme(axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=12)) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))
fig3
