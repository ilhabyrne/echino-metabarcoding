---
title: "RDA analyses for the Moore Reef metabarcoding dataset"
author: "Ilha Byrne"
date: "2024-02-21"
output: html_document
---
  
## Setup workspace
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(gridExtra)

setwd("/Volumes/Ilha/echinoMS/analyses/moore/rda")

### Load data and format
df <- read.csv("moore_counts_echinoSp_otutax_230823.csv")
df1 <- read.csv("moore_counts_echinoGenus_otutax_230823.csv")

#df <- df %>%
#tibble::column_to_rownames("BOLDspecies") 
#df <- lapply(df, as.numeric)
#df <- as.data.frame(df) 
#head(df)

df1 <- df1 %>%
  group_by(BOLDgenus) %>%
  summarise(across(everything(), sum), .groups = "drop")

df1 <- df1 %>%
  tibble::column_to_rownames("BOLDgenus") 
df1 <- as.data.frame(df1) 
df1 <- t(df1)
head(df1)

#### Load environmental data
env <- read_csv("moore_env.csv")
env <- env %>%
  tibble::column_to_rownames("sample")

env1 <- as.data.frame(env)
head(env1)

#### Subset environmental data
env2 <- read_csv("moore_env_subset.csv")
env2 <- env2 %>%
  tibble::column_to_rownames("sample")

env3 <- as.data.frame(env2)
head(env3)

pa <- decostand(df1, "pa") # standardize for presence absence
decorana(pa) # SD should be < 3
#hel <- decostand(pa, "hellinger") # Hellinger transformation if SD is too high (Legendre and Gallagher 2001)
#decorana(hel)

### Additional filtering - yields better results
sum <- apply(pa, 2, sum) # calculate sum per species
sort(sum)
fin <- pa[, ! sum<2] # remove species that occur in fewer than 2 samples
fin <- fin[rowSums(fin[, -1])>0,] # remove sites with 0 species
sort(apply(fin, 2, max)) # variability in max "abundance"
sort(apply(fin, 2, sd)) # variation is not high
decorana(fin) # SD should be < 3

## RDA analyses

### Unconstrained (PCA)
m1 <- rda(pa ~ 1, data=env1, scale=TRUE) # RDA analysis, scale=TRUE if not hellinger transformed
summary(m1) # % explained by variables in constrained 
plot(m1)
envfit <- envfit(m1, env1)

plot(m1)
plot(envfit, add = TRUE)

### Constrained model
m2 <- rda(pa ~ Year + Season + Temperature + Salinity, data=env1, scale=TRUE) # RDA analysis, scale=TRUE if not hellinger transformed
summary(m2) # % explained by variables in constrained 
plot(m2)

coef(m2) # extract canonical coefficients 
set.seed(111)
anova.cca(m2, STEP=1000) # checks if model is statistically significant 
anova.cca(m2, by="axis", step=10000) # statistical significance of each axis
anova.cca(m2, by="term", step=10000) # statistical significance of each term
vif.cca(m2) # anything above 10/20 should be avoided
RsquareAdj(m2)$adj.r.squared

#### Plot with ggplot
spp <- as.data.frame(scores(m2, display = "species"))
site <- as.data.frame(scores(m2, display = "sites"))

loadings <- data.frame(
  Variables = c("Year", "Season", "Temperature", "Salinity"),
  RDA1 = c(0.8130, 0.2807, 0.3651, -0.3582),
  RDA2 = c(0.4276, -0.9424, 0.7540, 0.1742)
)

m2.plot <- ggplot() +
  geom_point(data = site, aes(x = RDA1, y = RDA2), color = "grey50") +
  geom_text_repel(data = spp, aes(x = RDA1, y = RDA2, label = rownames(spp)),
                  size = 3, box.padding = 0.1, point.padding = 0.1,
                  segment.color = 'grey50', color = "blue", segment.size = 0,
                  max.iter = 5000) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
               color = "grey", linewidth = 0.5) +
  geom_text(data = loadings, aes(x = RDA1, y = RDA2, label = Variables), hjust = 1.2, vjust = 1.2) +
  labs(x = "RDA1", y = "RDA2") +
  theme_bw() +
  coord_fixed()

ggsave("moore-rda-full.png", plot = m2.plot, width = 8, height = 8, dpi = 300)


### Constrained model - version 2
m3 <- rda(fin ~ Year + Season + Temperature + Salinity, data=env3, scale=TRUE) # RDA analysis, scale=TRUE if not hellinger transformed
summary(m3) # % explained by variables in constrained 
plot(m3)

coef(m3) # extract canonical coefficients 
set.seed(111)
anova.cca(m3, STEP=1000) # checks if model is statistically significant 
anova.cca(m3, by="axis", step=10000) # statistical significance of each axis
anova.cca(m3, by="term", step=10000) # statistical significance of each term
vif.cca(m3) # anything above 10/20 should be avoided
RsquareAdj(m3)$adj.r.squared

#### Plot with ggplot
spp <- as.data.frame(scores(m3, display = "species"))
site <- as.data.frame(scores(m3, display = "sites"))

loadings1 <- data.frame(
  Variable = c("Year", "Season", "Temperature", "Salinity"),
  RDA1 = c(-0.1046792, -0.2668817, -0.1678972, -0.0308496),
  RDA2 = c(-0.1824914, 0.3812556, 0.3936271, 0.8709725)
)

m3.plot <- ggplot() +
  geom_point(data = site, aes(x = RDA1, y = RDA2), color = "grey50") +
  geom_text_repel(data = spp, aes(x = RDA1, y = RDA2, label = rownames(spp)),
                  size = 3, box.padding = 0.1, point.padding = 0.1,
                  segment.color = 'grey50', color = "blue", segment.size = 0,
                  max.iter = 5000) +
  geom_segment(data = loadings1, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
               color = "grey", linewidth = 0.5) +
  geom_text_repel(data = loadings1, aes(x = RDA1, y = RDA2, label = Variable),
                  box.padding = 0.1, point.padding = 0.1, segment.size = 0.2,
                  segment.color = "grey50", nudge_x = -0.1, nudge_y = 0.0001) +
  labs(x = "RDA1", y = "RDA2") +
  labs(x = "RDA1", y = "RDA2") +
  theme_bw() +
  coord_fixed()

m3.plot

ggsave("moore-rda-subset.png", plot = m3.plot, width = 8, height = 8, dpi = 300)

### Plot side by side
install.packages("patchwork")
library(patchwork)

combined <- (m2.plot | m3.plot) + 
  plot_annotation(tag_levels = 'a') & 
  theme(
    plot.tag = element_text(face = 'bold', size = 16),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt") # adjust top, right, bottom, left margins
  )

combined # view multi-panel figure

ggsave("moore-rda-combined.png", plot = combined, width = 8, height = 8, dpi = 300)

### Reduced model
#m3 <- rda(pa ~ Salinity + Temperature, data=env1, scale=TRUE) # RDA analysis, scale=TRUE if not hellinger transformed
#summary(m3) # % explained by variables in constrained 
#plot(m3)

#coef(m3) # extract canonical coefficients 
#set.seed(111)
#anova.cca(m3, STEP=1000) # checks if model is statistically significant 
#anova.cca(m3, by="axis", step=10000) # statistical significance of each axis
#anova.cca(m3, by="term", step=10000) # statistical significance of each term
#vif.cca(m3) # anything above 10 should be avoided
#RsquareAdj(m3)$adj.r.squared
