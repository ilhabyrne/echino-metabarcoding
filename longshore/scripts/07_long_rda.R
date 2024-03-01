---
title: "RDA analyses for the longshore metabarcoding dataset"
author: "Ilha Byrne"
date: "2024-02-21"
output: html_document
---
  
## Setup workspace
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)

setwd("/Volumes/Ilha/echinoMS/analyses/moore/rda")

### Load data and format
df <- read.csv("long_counts_echinoGenus_otutax_subset_230823.csv")

df <- df %>%
  group_by(genus) %>%
  summarise(across(everything(), sum), .groups = "drop")

df <- df %>%
  tibble::column_to_rownames("genus") 
df <- as.data.frame(df) 
df <- t(df)
head(df)

#### Load environmental data
env <- read_csv("long_env.csv")
env <- env %>%
  tibble::column_to_rownames("sample")

env1 <- as.data.frame(env)
head(env1)

#### Subset environmental data
env2 <- read_csv("long_env_subset.csv")
env2 <- env2 %>%
  tibble::column_to_rownames("label")

env3 <- as.data.frame(env2)
head(env3)

### Check data
pa <- decostand(df, "pa") # standardize for presence absence
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
m2 <- rda(pa ~ year + site + temperature + salinity + chl, data=env1, scale=TRUE) # RDA analysis, scale=TRUE if not hellinger transformed
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
  Variables = c("Year", "Site", "Temperature", "Salinity", "Chl-a"),
  RDA1 = c(-0.58356, 0.21859, 0.72818, -0.01864, 0.67015),
  RDA2 = c(-0.3407, -0.526, -0.3911, 0.48, 0.576)
)

m2.plot <- ggplot() +
  geom_point(data = site, aes(x = RDA1, y = RDA2), color = "grey50") +
  geom_text_repel(data = spp, aes(x = RDA1, y = RDA2, label = rownames(spp)),
                  size = 3, box.padding = 0.5, point.padding = 0,
                  segment.color = 'grey50', color = "blue", segment.size = 0,
                  max.iter = 5000) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
               color = "grey", linewidth = 0.55) +
  geom_text(data = loadings, aes(x = RDA1, y = RDA2, label = Variables),
            nudge_x = 0.07, nudge_y = 0.1, hjust = 0, vjust = 0, color = "black") +
  labs(x = "RDA1 (13.8%)", y = "RDA2 (8.0%)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.tag = element_text(size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed()

m2.plot

ggsave("long-rda-full.png", plot = m2.plot, width = 8, height = 8, dpi = 300)


### Constrained model - version 2
m3 <- rda(fin ~ year + site + temperature + salinity + chl, data=env3, scale=TRUE) # RDA analysis, scale=TRUE if not hellinger transformed
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
spp1 <- as.data.frame(scores(m3, display = "species"))
site1 <- as.data.frame(scores(m3, display = "sites"))

loadings1 <- data.frame(Variables = c("Year", "Site", "Temperature", "Salinity", "Chl-a"),
                        RDA1 = c(-0.76184, -0.02123, 0.20102, 0.2235, 0.94215),
                        RDA2 = c(-0.19552, 0.68426, 0.43079, -0.32157, -0.01841)
)


m3.plot <- ggplot() +
  geom_point(data = site1, aes(x = RDA1, y = RDA2), color = "grey50") +
  geom_text_repel(data = spp1, aes(x = RDA1, y = RDA2, label = rownames(spp1)),
                  size = 3, box.padding = 0.1, point.padding = 0.1,
                  segment.color = 'grey50', color = "blue", segment.size = 0,
                  max.iter = 5000) +
  geom_segment(data = loadings1, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
               color = "grey", linewidth = 0.5) +
  geom_text(data = loadings1, aes(x = RDA1, y = RDA2, label = Variables), hjust = -0.1, vjust = 1.2) +
  labs(x = "RDA1 (22.7%)", y = "RDA2(12.5%)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        plot.tag = element_text(size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed()

m3.plot

ggsave("long-rda-subset.png", plot = m3.plot, width = 8, height = 8, dpi = 300)

### Plot side by side
library(patchwork)

combined <- m2.plot + m3.plot + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = "") 

combined # view multi-panel figure

ggsave("long-rda-combined.png", plot = combined, width = 8, height = 8, dpi = 300)

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
