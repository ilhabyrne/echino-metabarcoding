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
df <- read.csv("long_counts_echinoSp_otutax.csv")
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

fig1 <- ggplot(df_melt, aes(variable, species, fill=value)) +
  geom_tile(color="grey") + 
  labs(x="Site", y="Species", fill="# OTUs") +
  scale_fill_gradient(low="white", high="blue") +
  theme(axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=12)) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))
fig1

### Load genus data
df1 <- read.csv("long_counts_echinoGenus_otutax.csv")
df1 <- as.data.frame(df1) 
head(df1)

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

fig2 <- ggplot(df1_melt, aes(variable, genus, fill=value)) +
  geom_tile(color="grey") + 
  labs(x="Site", y="Genera", fill="# OTUs") +
  scale_fill_gradient(low="white", high="blue", breaks=c(0,1)) +
  theme(axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=12)) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1)) +
  theme(axis.title = element_text(colour="black", size=16))
fig2

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
