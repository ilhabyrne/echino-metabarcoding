# echino_metabarcoding

Welcome to the GitHub repository associated with my Honours project conducted at the University of Queensland as part of the Bachelor of Science (Honours) program (August 2020 - May 2021). In this repository, you can find the raw data associated with this metabarcoding project and the associated scripts and processed files used to produce the final results.

## Project overview
The main aim of this project was to use metabarcoding as a tool to identify the echinoderm larvae present in plankton samples (generously donated by the Australian Institute of Marine Science) and speak to the temporal and spatial variability in their occurrence.

## Bioinformatics pipeline

## Post-bioinformatics wrangling
Below you will find my ramblings on how I tackled the post-bioinfromatics processing of this dataset. This is just an overview of what I did - refer to my scripts to see exactly what I did. 

Main steps
- [x]  Merge MIDORI sintax & echino NCBI results
- [x]  Remove singletons
- [x]  Format for use with phyloseq
- [x]  Remove anything without hits (but keep both file versions)
- [x]  Merge PCR replicates & save as separate files
- [x]  Convert to count data & save as separate file

### Data wrangling in R

#### Merge MIDORI & NCBI results

Remove singletons or filter based on a minimum number of reads. To begin, I simply just removed singletons. 

```r
# Remove singletons
tab1 <- read.csv(".csv")
tab1[tab1 < 1] <- 0 # remove hits with < 1 reads

# If total reads < 10 then disregard taxon
tab2[tab2 < 10] <- 0
write.csv(tab2, ".csv")
```

#### Format for phyloseq using excel

**OTU matrix:**

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/3dd72f14-5492-4cd2-aba0-158fcf2146b0/Untitled.png)

**Taxa matrix:**

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/a2d2bb02-607b-46ea-9535-174fccc588e2/Untitled.png)

**Metadata:**

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/2bc86109-ad8f-4cfb-82cb-103062c55535/Untitled.png)

#### Import & format phyloseq objects in R

```r
# Load pyloseq
install.packages("phyloseq")
library(phyloseq)

otumat <- read.csv(".csv")
taxmat <- read.csv(".csv")

# Convert OTU column to row names
otumat1 <- otumat[,-1]
rownames(otumat1) <- otumat[,1]

taxmat1 <- taxmat[,-1]
rownames(taxmat1) <- taxmat[,1]
taxmat2 <- as.matrix(taxmat1)

OTU <- otu_table(otumat1, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat2)

physeq = phyloseq(OTU,TAX)

sampledata <- read.csv("sampledata-phyloseq.csv")
sampledata1 <- sampledata[,-1]
rownames(sampledata1) <- sampledata[,1]

sampledata2 <- sample_data(sampledata1)

physeq1 = phyloseq(OTU, TAX, sampledata2)
```

#### Prune dataset

Do this in excel and ensure you save each version with a meaningful file name. 

### Exploratory barplots

#### Stacked barplot of the main phyla

```r
# Stacked barplot
fig1 <- plot_bar(physeq1, fill = "Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  labs(x = "Sample", y = "Read counts") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25000)) +
  theme_classic() +
  theme(legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10), 
        axis.text.x = element_text(size = 12, colour = "black", angle=50, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 16))

fig1$data$Sample <- factor(fig1$data$Sample, levels = c("site1","site2"))

fig1
```

### Sample and taxa specific filtering

#### Maximum contamination filtering

To perform maxCont filtering, you need to format the data frame in excel first (or in R):

- [ ]  make sure you are working on the data frame that has already had singletons removed
- [ ]  calculate maximum NC reads ← NCmax column
- [ ]  remove NC columns

```r
# Example code using a random data frame:

# create a sample data frame
df <- data.frame(species = c("A", "B", "C"), sample1 = c(10, 15, 5), sample2 = c(8, 4, 12), sample3 = c(20, 2, 7), min_values = c(7, 5, 10))

# set the row names as the species names
rownames(df) <- df$species

# select all columns except for the minimum value column
df_subset <- df[, -ncol(df)]

# create a function to filter values based on minimum threshold
threshold_filter <- function(row, threshold) {
  ifelse(row <= threshold, 0, row)
}

# apply the threshold filter to each row using mapply()
df_filtered <- mapply(threshold_filter, df_subset, df$min_values)

# set the row names to NULL to remove the species names
rownames(df_filtered) <- NULL

# view the filtered data frame
df_filtered
```

#### Sample-specific filtering

**Longshore dataset: sample % threshold** 
For the longshore dataset, PCR replicates were not uniquely tagged but rather pooled. As such, we cannot merge/filter based on variability among replicates. Here, I used a sample percentage threshold as per Drake et al. (2022). 

```r

```

**Moore dataset: PCR replicate merging** 

For the Moore reef dataset, merge PCR/technical replicates using a restrictive additive approach (see Alberdi et al. 2018). To do this step, ensure:

- [x]  PCR replicate names are structured for R i.e. sample_1, sample_2, sample_3, sample1_1, …

```r
# Merge PCR replicates
library(plyr)
tab1 <- read.csv(".csv")
tab2 <- ddply(tab1, "sample_uniq", numcolwise(sum))

LP_otu <- read_csv("leray_pilot_count_sxs.csv") # make sure to convert values to numbers first
LP_otutab <- ddply(LP_otu, "taxonomy", numcolwise(sum)) %>% mutate_if(is.numeric, ~1 * (. > 2)) # sum technical reps, then convert to counts if otu present in >2 replicates
write.csv(LP_otutab, "leray_pilot_count_reps_sxs.csv")
```

### Convert to count data

```r
# Convert to count data 
tab1c <- tab1 %>% mutate_if(is.numeric, ~1 * (. > 0))
```

### Biological replicate merging

```r

```

## Final figures

### Heatmaps

#### Format files

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/01e8a0e0-32f2-483c-985d-e6e2484d9c7f/Untitled.png)

## Statistical analyses

### glm(m) model

#### Format files

| sample | month | season | year | Acanthaster | Asteroidea | Crinoidea | Echinodermata | Echinoidea | Holothuroidea | Ophiuroidea | echino.rich |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Aug_2016 | 8 | Winter | 2016 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Aug_2016 | 8 | Winter | 2016 | 1 | 1 | 0 | 0 | 4 | 0 | 0 | 5 |
| Aug_2017 | 8 | Winter | 2017 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Aug_2017 | 8 | Winter | 2017 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### NMDS model

#### Format files

otutab:

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/316fa5af-2373-4bb8-870c-ccb97404bfa8/Untitled.png)

metadata:

## References

### Acknowledgements 
A MASSIVE thank you to Dr Iva Popovic, Dr Dean Brookes and Professor Cynthia Riginos who were instrumental in this study and without whom none of this work would have beeen achievable. I also want to thank Dr Sven Uthicke and his entire team at AIMS for all their assistance in collecting samples, sample processing etc. Dr Karin Zwiep, Dr Simone Blomberg and the Queensland Cyber Infrastructure Foundation (QCIF) also provided useful comments and advice. 
