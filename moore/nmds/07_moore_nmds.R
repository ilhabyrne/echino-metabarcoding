## PERMANOVA
  
#Assumes no distribution, allows for differences in between-group variation, is insensitive to multicollinearity, allows for multiple variables and is insensitive to many zeros.

#Assumptions:
##Objects in the data set are exchangeable under the null hypothesis
##Exchangable objects (sites, samples, observations, etc) are independent
##Exchangable objects have similar multivariate disperson (i.e. each group has a similar degree of multivariate scatter)


## Load workspace
library(vegan)
setwd("~/Desktop/Undergraduate/manuscript/hons_manuscript_2023/moore/nmds")

## Format otutab for nmds
df <- read_csv("moore_counts_echinoSp_otutax.csv")
df1 <- ddply(df, "species", numcolwise(sum)) 
df2 <- t(df1)
df2 <- as.data.frame(df2)
df2$row_names <- row.names(df2)

write_csv(df2,"moore_counts_echinoSp_otutax2.csv")


## NMDS analyses

### Load data
data <- read.csv("moore_counts_echinoSp_otutax2.csv") #row sum <0 need to be removed, removed Aug_2017
data <- data %>%
  tibble::column_to_rownames("species") 

data1 <- lapply(data,as.numeric)
data2 <- as.data.frame(data1) 
head(data2)

env <- read_csv("moore_env.csv") #removed Aug_2017
env <- env %>%
  tibble::column_to_rownames("sample") 

env1 <- as.data.frame(env)
head(env1)

### vegdist
dis <- vegdist(data2)
mod <- betadisper(dis, env1$season)
mod
plot(mod)

### test assumptions
t1 <- betadisper(dis, env1$season) #test for homogeneity of group dispersion 
anova(t1) #not significant

permutest(t1, permutations = 999, pairwise=TRUE)
TukeyHSD(t1)

t2 <- betadisper(dis, env1$year) #test for homogeneity of group dispersion 
anova(t2) #not significant

### plots

plot(t1, ellipse=FALSE, hull=FALSE, segments=FALSE, conf=0.28, lwd=2, cex=1)

### final model
adonis2(data2 ~ season*year, data=env1, perm=999, method="bray")

adonis2(data2 ~ season + year, data=env1, perm=999, method="bray")

adonis2(data2 ~ season, data=env1, perm=999, method="bray")

adonis2(data2 ~ year, data=env1, perm=999, method="bray")
