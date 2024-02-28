# Heatmaps of echino species at each time-point

### Set working directory
setwd("/Users/uqibyrne/Desktop/hons_manuscript_2023/moore/heatmap")

### Load packages required
install.packages("reshape")
library(reshape)
install.packages("dplyr")
library(dplyr)
library(ggplot2)
library(tidyverse)

### Format data for heatmpas (only need to do this once)
#### Merge taxonomy and otumat 
tax <- read.csv("moore_echino_tax.csv")
otu <- read.csv("moore_counts_sampleMerge.csv")

otutax <- merge(tax, otu, by.x = "otuID", by.y = "otuID")
  
#### Save new dataframe as csv
write_csv(otutax, "moore_counts_echino_otutax.csv")

### Load species data
df <- read.csv("moore_counts_echinoSp_otutax_23-08-23.csv")
df <- as.data.frame(df) 
head(df)

### Format data frame for plot
df_melt <- melt(df)
head(df_melt)

df_melt$variable <- factor(df_melt$variable, 
                             level = c("Dec_2015", "Jan_2016", "Mar_2016", 
                                       "Jun_2016", "Aug_2016", "Dec_2016", 
                                       "Jan_2017", "Mar_2017", "Jun_2017", 
                                       "Aug_2017", "Dec_2017", "Jan_2018", 
                                       "Mar_2018", "Nov_2019", "Jan_2020"))

df_melt$value <- as.factor(df_melt$value) 
                           
fig1 <- ggplot(df_melt, aes(variable, BOLDspecies, fill=value)) +
  geom_tile(color="grey", show.legend = FALSE) + 
  theme_classic() +
  labs(x="Sampling date", y="Species") +
  scale_fill_manual(values = c("white","#4075B4","#4075B4","#4075B4")) +
  annotate("rect", xmin=0.5, xmax=15.45, ymin=0.5, 
           ymax=1.5, colour="black", fill="transparent", linewidth=1) +
  theme(axis.text.x = element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10, face="italic"),
        axis.title = element_text(size=14, face="bold")) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))
fig1

fig1 <- fig1 + theme(axis.ticks = element_blank(),
             axis.line = element_blank())
fig1

### Load genus data
df1 <- read.csv("moore_counts_echinoGenus_otutax_23-08-23.csv")
df1 <- as.data.frame(df1) 
head(df1)

### Format data frame for plot
df1_melt <- melt(df1)
head(df1_melt)

df1_melt <- df1_melt %>%
  group_by(BOLDgenus,variable) %>%
  summarize(value=sum(value, na.rm=TRUE))

df1_melt$variable <- factor(df1_melt$variable, 
                           level = c("Dec_2015", "Jan_2016", "Mar_2016", 
                                     "Jun_2016", "Aug_2016", "Dec_2016", 
                                     "Jan_2017", "Mar_2017", "Jun_2017", 
                                     "Aug_2017", "Dec_2017", "Jan_2018", 
                                     "Mar_2018", "Nov_2019", "Jan_2020"))

df1_melt$value <- as.factor(df1_melt$value) 

fig2 <- ggplot(df1_melt, aes(variable, BOLDgenus, fill=value)) +
  geom_tile(color="grey", show.legend = FALSE) + 
  theme_classic() +
  labs(x="Sampling date", y="Genera", fill="OTU counts") +
  scale_fill_manual(values=c("white","#4075B4","#4075B4","#4075B4", "#4075B4")) +
  annotate("rect", xmin=0.5, xmax=15.45, ymin=0.55, 
               ymax=1.5, colour="black", fill="transparent", size=1) +
  annotate("rect", xmin=0.5, xmax=15.45, ymin=8.5, 
           ymax=9.5, colour="black", fill="transparent", size=1) +
  scale_y_discrete(expand=c(0, 0)) +
  theme(axis.text.x = element_text(colour="black", size=10, angle = 50, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour="black", size=10, face="italic"),
        axis.title = element_text(size=14, face="bold"),
        legend.title = element_text(size=14, face="bold"))
fig2

fig2 <- fig2 + theme(axis.ticks = element_blank(),
             axis.line = element_blank())
fig2

## Prep panel for publication 
library(patchwork)

### Nest panels together
nested <- (fig1/fig2)+
  plot_annotation(tag_levels = 'a') & # add figure labels 
  theme(plot.tag = element_text(face = 'bold', size = 16)) # edit text
nested # view multi-panel figure

## Save at high-res for publication
ggsave(# _plot name_ if you have it, else it will save the last one printed 
  filename = "moore_tile_main_20x25.png", 
  dpi = 300, # for publications  
  units = "cm", 
  width = 20, height = 25,
  bg = "white")
