# Statistical analyses

## Load packages
library(glmmTMB)
library(DHARMa)
library(multcomp)
library(emmeans)
library(readxl)
library(tidyverse)
library(vegan)

## Merge taxonomy and otumat
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/moore/glmm")

tax <- read.csv("06a_moore_counts_pcrMerge_phyloseq.csv")
otu <- read.csv("moore_subsetTaxa_phyloseq.csv")

otutax <- merge(tax, otu, by.x = "otuID", by.y = "otuID")

### Save new dataframe as csv
write_csv(otutax, "moore_counts_otutax.csv")

## Format otutax for glmm
df <- read.csv("moore_countsPhyla_otutax.csv")
df <- as.data.frame(df) 
head(df)

df <- ddply(df, "phylum", numcolwise(sum)) 

df_melt <- melt(df)
head(df_melt)

df_t <- as.data.frame(t(df))
head(df_t)

### Save new dataframe as csv
write_csv(df_melt, "moore_countsPhyla_list.csv")
write_csv(df_t, "moore_countsPhyla_otutax2.csv")


## GLMM models

### Load data
rich <- read_csv("moore_metadata.csv")
rich$season <- as.factor(rich$season)

### Run models and assess AIC value
m1 <- glmmTMB(echino.rich ~ season, family="poisson", data=rich)
summary(m1) #AIC:142.7

m2 <- glmmTMB(echino.rich ~ season + year, family="poisson", data=rich)
summary(m2) #AIC:143.2; year not significant 

m3 <- glmmTMB(total.rich ~ season, family="poisson", data=rich)
summary(m3) #AIC:516.5

m4 <- glmmTMB(total.rich ~ season + year, family="poisson", data=rich)
summary(m4) #AIC:506

### Seasons means test
testsm1 <- glht(m1, linfct=mcp(season="Tukey"))
summary(testsm1)

meansm1 <- emmeans(m1, specs = "season")
pairs(meansm1)

plot(meansm1, horizontal = FALSE)

### Check the dispersion and homogeneity of best models 
res1 <- simulateResiduals(m1)
plot(res1)

## GLM models

m1 <- glm(echino.rich ~ season, data=rich, family="poisson")
summary(m1) #AIC:142.69

### Seasons means test
testsm1 <- glht(m1, linfct=mcp(season="Tukey"))
summary(testsm1)

meansm1 <- emmeans(m1, specs = "season")
pairs(meansm1)

plot(meansm1, horizontal = FALSE)

### Check the dispersion and homogeneity of best models 
res1 <- simulateResiduals(m1)
plot(res1)
