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
otu <- read.csv("moore_subsetTaxa_echinoBOLD_phyloseq.csv") # 25/08/23

otutax <- merge(tax, otu, by.x = "otuID", by.y = "otuID")

### Save new dataframe as csv
write_csv(otutax, "moore_counts_otutax_250823.csv")

## Format otutax for glmm
df <- read.csv("moore_countsPhyla_otutax_250823.csv")
df <- as.data.frame(df) 
head(df)

library(plyr)
library(reshape2)
df <- ddply(df, "phylum", numcolwise(sum)) 

df_melt <- melt(df)
head(df_melt)

df_t <- as.data.frame(t(df))
head(df_t)
df_t$row_names <- row.names(df_t)
head(df_t)


### Save new dataframe as csv
write_csv(df_melt, "moore_countsPhyla_list_250823.csv")
write_csv(df_t, "moore_countsPhyla_otutax2_250823.csv")


df <- read_csv("moore_countsPhyla_otutax2_250823.csv")
df1 <- as.data.frame(df)
df2 <- t(df1)

write_csv(df1, "moore_glmm_250823.csv")

# GLMM models

## Load data
rich <- read_csv("moore_glmm_250823.csv")
rich$season <- as.factor(rich$season)

## Run models and assess AIC value
m1 <- glmmTMB(echino.rich ~ season, family="poisson", data=rich)
summary(m1) #AIC:168.1

m2 <- glmmTMB(echino.rich ~ season + year, family="poisson", data=rich)
summary(m2) #AIC:159.4; year not significant 

m3 <- glmmTMB(total.rich ~ season, family="poisson", data=rich)
summary(m3) #AIC:508.2

m4 <- glmmTMB(total.rich ~ season + year, family="poisson", data=rich)
summary(m4) #AIC:501.5

### Seasons means test
testsm1 <- glht(m1, linfct=mcp(season="Tukey"))
summary(testsm1)

meansm1 <- emmeans(m1, specs = "season")
pairs(meansm1)

plot(meansm1, horizontal = FALSE)

### Check the dispersion and homogeneity of best models 
res1 <- simulateResiduals(m1)
plot(res1)

# GLM models

## Model 1
m1 <- glm(echino.rich ~ season, data=rich, family="poisson")
summary(m1) #AIC:168.14

### Seasons means test
testsm1 <- glht(m1, linfct=mcp(season="Tukey"))
summary(testsm1)

meansm1 <- emmeans(m1, specs = "season")
pairs(meansm1)

plot(meansm1, horizontal = FALSE)

### Check the dispersion and homogeneity of best models 
res1 <- simulateResiduals(m1)
plot(res1)

## Model 2
m2 <- glm(echino.rich ~ season + year, data=rich, family="poisson")
summary(m2) #AIC:159.38.14

### Seasons means test
testsm2 <- glht(m2, linfct=mcp(season="Tukey"))
summary(testsm2)
