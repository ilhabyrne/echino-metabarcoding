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
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/longshore/glmm")

tax <- read.csv("06a_long_counts_percentThresh_phyloseq.csv")
otu <- read.csv("03_long_subsetTaxa_phyloseq.csv")

otutax <- merge(tax, otu, by.x = "otuID", by.y = "otuID")

### Save new dataframe as csv
write_csv(otutax, "long_counts_otutax.csv")

## Format otutax for glmm
df <- read.csv("long_countsPhyla_otutax.csv")
df <- as.data.frame(df) 
head(df)

df <- ddply(df, "phylum", numcolwise(sum)) 

df_melt <- melt(df)
head(df_melt)

df_t <- as.data.frame(t(df))
df_t$row_names <- row.names((df_t))
head(df_t)

### Save new dataframe as csv
write_csv(df_melt, "long_countsPhyla_list.csv")
write_csv(df_t, "long_countsPhyla_otutax2.csv")


## GLMM models

### Load data
rich <- read_csv("longshore_metadata.csv")
rich$year <- as.factor(rich$year)
rich$site <- as.factor(rich$site)

### Run models and assess AIC value
m1 <- glmmTMB(echino.rich ~ year, family="poisson", data=rich)
summary(m1) #AIC:243.5

m2 <- glmmTMB(echino.rich ~ year + site, family="poisson", data=rich)
summary(m2) #AIC:221.6

m3 <- glmmTMB(echino.rich ~ year + region, family="poisson", data=rich)
summary(m3) #AIC:245.3 

### Years means test
testsm2 <- glht(m2, linfct=mcp(year="Tukey"))
summary(testsm1)

meansm2 <- emmeans(m2, specs = "year")
pairs(meansm2)

### Check the dispersion and homogeneity of best models 
res2 <- simulateResiduals(m2)
plot(res2)

## GLM models

m1 <- glm(echino.rich ~ year, data=rich, family="poisson")
summary(m1) #AIC:243.49

m2 <- glm(echino.rich ~ year + site, family="poisson", data=rich)
summary(m2) #AIC:221.56

### Years means test
testsm2 <- glht(m2, linfct=mcp(year="Tukey"))
summary(testsm2)

meansm2 <- emmeans(m2, specs = "year")
pairs(meansm2)

### Sites means test
testsm2 <- glht(m2, linfct=mcp(site="Tukey"))
summary(testsm2)

### Check the dispersion and homogeneity of best models 
res2 <- simulateResiduals(m2)
plot(res2)
