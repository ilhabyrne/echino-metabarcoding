setwd("~/Desktop/Undergraduate/manuscript/hons_manuscript_2023/moore")

df <- read.csv("tech-reps.csv")

library(reshape)
df_melt <- melt(df)

write.csv(df_melt, "tech-reps.csv")

df <- read.csv("tech-reps.csv")

library(ggplot2)

df$replicate <- as.factor(df$replicate)

ggplot(df, aes(x=replicate, y=total)) +
  geom_boxplot(fill="lightgrey") +
  scale_y_continuous(breaks = seq(0, 32000, 4000), limits = c(0, 32000)) +
  theme_classic() +
  labs(x="Sample", y = "Number of reads") +
  theme(axis.text.x = element_text(colour="black", size=12,angle = 50, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour="black", size=12)) +
  theme(axis.title = element_text(size = 16))

