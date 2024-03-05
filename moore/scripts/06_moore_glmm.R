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

### ~~~~~~~~~~ ECHINODERMS ONLY ~~~~~~~~~~ ###

# GLMM models

## Load data
setwd("~/Documents/GitHub/echino-metabarcoding/moore/glmm")
rich <- read_csv("moore_glmm_250823.csv")
rich$season <- as.factor(rich$season)
rich$year <- as.factor(rich$year)

## Run models and assess AIC value
m1 <- glmmTMB(echino.rich ~ season, family="poisson", data=rich)
summary(m1) #AIC:168.1

m1a <- glmmTMB(echino.rich ~ season + (1 | sample), family="poisson", data=rich)
summary(m1a) #AIC:157.4

m1b <- glmmTMB(echino.rich ~ year + (1 | sample), family="poisson", data=rich)
summary(m1b) #AIC:153.4

m2a <- glmmTMB(echino.rich ~ season + year + (1 | sample), family="poisson", data=rich)
summary(m2a) #AIC:149.8; year not significant 

### Seasons means test
testsm1a <- glht(m1a, linfct=mcp(season="Tukey"))
summary(testsm1a)

testsm1b <- glht(m1b, linfct=mcp(year="Tukey"))
summary(testsm1b)

plot(meansm1, horizontal = FALSE)

### Check the dispersion and homogeneity of model 
res1 <- simulateResiduals(m1)
plot(res1) # bad

res1a <- simulateResiduals(m1a)
plot(res1a) # good

res2 <- simulateResiduals(m2)
plot(res2) # bad

res2a <- simulateResiduals(m2a)
plot(res2a) # good

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

### Check the dispersion and homogeneity of best model
res1 <- simulateResiduals(m1)
plot(res1) # significant dispersion

## Model 2
m2 <- glm(echino.rich ~ season + year, data=rich, family="poisson")
summary(m2) #AIC:149.08

### Check the dispersion and homogeneity of model
res2 <- simulateResiduals(m2)
plot(res2) # significant dispersion

### ~~~~~~~~~~ ALL TAXA ~~~~~~~~~~ ###

## Run models and assess AIC value
m1 <- glmmTMB(total.rich ~ season, family="poisson", data=rich)
summary(m1) #AIC:508.2

m2 <- glmmTMB(total.rich ~ season + year, family="poisson", data=rich)
summary(m2) #AIC:458; year not significant 

m3 <- glmmTMB(total.rich ~ season + (1 | sample), family="poisson", data=rich)
summary(m3) #AIC:309.7

m4 <- glmmTMB(total.rich ~ year + (1 | sample), family="poisson", data=rich)
summary(m4) #AIC:313.9

m5 <- glmmTMB(total.rich ~ season + year + (1 | sample), family="poisson", data=rich)
  # rank-deficient conditional model
summary(m5) #AIC:311.2

### Check the dispersion and homogeneity of models
res3 <- simulateResiduals(m3)
plot(res3) # good

res4 <- simulateResiduals(m4)
plot(res4) # good

res5 <- simulateResiduals(m5)
plot(res5) # good

### Seasons means test
testsm3 <- glht(m3, linfct=mcp(season="Tukey"))
summary(testsm3)

testsm4 <- glht(m4, linfct=mcp(year="Tukey"))
summary(testsm4)

### ~~~~~~~~~~ TAXA relationships ~~~~~~~~~~ ###

m7 <- glmmTMB(echino.rich ~ Cnidaria + (1 | sample), family="poisson", data=rich)
summary(m7)

m7a <- lm(echino.rich ~ Cnidaria, data=rich)
summary(m7a)

res7 <- simulateResiduals(m7)
plot(res7) # good
