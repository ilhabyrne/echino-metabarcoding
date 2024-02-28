# Set working directory
setwd("~/Desktop/hons_manuscript_2023/moore/variation")

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

# Plot as boxplot
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
fig1
