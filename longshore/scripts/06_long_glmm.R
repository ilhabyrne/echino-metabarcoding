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
setwd("/Users/ilhabyrne/Documents/GitHub/echino-metabarcoding/longshore/glmm")

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

### Load data
setwd("~/Documents/GitHub/echino-metabarcoding/longshore/glmm")
rich <- read.csv("longshore_metadata.csv")
rich$year <- as.factor(rich$year)
rich$site <- as.factor(rich$site)
rich$region <- as.factor(rich$region)

### ~~~~~~~~~~ ECHINODERMS ONLY ~~~~~~~~~~ ###

## GLMM models

### Run models and assess AIC value
m1 <- glmmTMB(echino.rich ~ year + (1|sample), family="poisson", data=rich)
summary(m1) #AIC:223.6

m2 <- glmmTMB(echino.rich ~ year + site + (1|sample), family="poisson", data=rich)
summary(m2) #AIC:221.3

m3 <- glmmTMB(echino.rich ~ year + region + (1|sample), family="poisson", data=rich)
summary(m3) #AIC:225.6 

m4 <- glmmTMB(echino.rich ~ year + site + region + (1|sample), family="poisson", data=rich)
summary(m4) #AIC:221.3

### Years means test
testsm2 <- glht(m2, linfct=mcp(site="Tukey"))
summary(testsm2)


### Check the dispersion and homogeneity of best models 
res2 <- simulateResiduals(m2)
plot(res2) # good

## GLM models

m1 <- glm(echino.rich ~ year, data=rich, family="poisson")
summary(m1) #AIC:243.49

m2 <- glm(echino.rich ~ year + site, family="poisson", data=rich)
summary(m2) #AIC:221.56

### Years means test
testsm2 <- glht(m2, linfct=mcp(year="Tukey"))
summary(testsm2)

### Sites means test
testsm2 <- glht(m2, linfct=mcp(site="Tukey"))
summary(testsm2)

### Check the dispersion and homogeneity of best models 
res2 <- simulateResiduals(m2)
plot(res2)

### ~~~~~~~~~~ ALL TAXA ~~~~~~~~~~ ###

### Run models and assess AIC value
m1 <- glmmTMB(total.rich ~ year + (1|sample), family=nbinom1, data=rich)
summary(m1) #AIC:637.9

m2 <- glmmTMB(total.rich ~ year + site + (1|sample), family=nbinom1, data=rich)
summary(m2) #AIC:638.8

m3 <- glmmTMB(total.rich ~ year + region + (1|sample), family=nbinom1, data=rich)
summary(m3) #AIC:634.3

m4 <- glmmTMB(total.rich ~ year + site + region + (1|sample), family=nbinom1, data=rich)
summary(m4) #AIC:638.8

### Check the dispersion and homogeneity of best models 
res1 <- simulateResiduals(m1)
plot(res1) # good

res2 <- simulateResiduals(m2)
plot(res2) # bad

res3 <- simulateResiduals(m3)
plot(res3) # good

res4 <- simulateResiduals(m4)
plot(res4) # good

### Years means test
testsm3 <- glht(m3, linfct=mcp(year="Tukey"))
summary(testsm3)

### Regions means test
testsm3a <- glht(m3, linfct=mcp(region="Tukey"))
summary(testsm3a)

### ~~~~~~~~~~ TAXA COMPARISONS ~~~~~~~~~~ ###

m5 <- glmmTMB(echino.rich ~ Cnidaria + (1|sample), family=nbinom1, data=rich)
summary(m5) 

res5 <- simulateResiduals(m5)
plot(res5) # good
