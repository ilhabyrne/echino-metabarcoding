## Heatmaps of echino species at each time-point

### Load packages required
install.packages("reshape")
library(reshape)

### Merge taxonomy and otumat
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/moore")
tax <- read.csv("moore_echino_tax.csv")
otu <- read.csv("moore_counts_sampleMerge.csv")

otutax <- merge(tax, otu, by.x = "otuID", by.y = "otuID")
  
### Save new datafrmae as csv
write_csv(otutax, "moore_counts_echino_otutax.csv")

### Load species data
df <- read.csv("moore_counts_echinoSp_otutax.csv")
df <- as.data.frame(df) 
head(df)

df <- ddply(df, "species", numcolwise(sum)) 

### Format data frame for plot
df_melt <- melt(df)
head(df_melt)

df_melt$variable <- factor(df_melt$variable, 
                             level = c("Dec_2015", "Jan_2016", "Mar_2016", 
                                       "Jun_2016", "Aug_2016", "Dec_2016", 
                                       "Jan_2017", "Mar_2017", "Jun_2017", 
                                       "Aug_2017", "Dec_2017", "Jan_2018", 
                                       "Mar_2018", "Nov_2019", "Jan_2020"))

fig1 <- ggplot(df_melt, aes(variable, species, fill=value)) +
  geom_tile(color="grey") + 
  labs(x="Date", y="Species", fill="# OTUs") +
  scale_fill_gradient(low="white", high="blue") +
  theme(axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=12)) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))
fig1

### Load genus data
df1 <- read.csv("moore_counts_echinoGenus_otutax.csv")
df1 <- as.data.frame(df1) 
head(df1)

df1 <- ddply(df1, "genus", numcolwise(sum)) 
head(df1)

### Format data frame for plot
df1_melt <- melt(df1)
head(df1_melt)

df1_melt$variable <- factor(df1_melt$variable, 
                           level = c("Dec_2015", "Jan_2016", "Mar_2016", 
                                     "Jun_2016", "Aug_2016", "Dec_2016", 
                                     "Jan_2017", "Mar_2017", "Jun_2017", 
                                     "Aug_2017", "Dec_2017", "Jan_2018", 
                                     "Mar_2018", "Nov_2019", "Jan_2020"))

fig2 <- ggplot(df1_melt, aes(variable, genus, fill=value)) +
  geom_tile(color="grey") + 
  labs(x="Date", y="Genera", fill="# OTUs") +
  scale_fill_gradient(low="white", high="blue") +
  theme(axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=12)) +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))
fig2
