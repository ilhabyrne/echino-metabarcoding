## Load packages
library(plyr)
library(tidyr)
library(dplyr)

install.packages("tidyverse")
library(tidyverse)

## Load data
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/moore")

df <- read.csv("moore_reads_filtered_phyloseq.csv")

## Max contamination filtering

# set the row names as the otu names
rownames(df) <- df$uniq

# select all columns except for the minimum value column
df_subset <- df[, -ncol(df)]

# create a function to filter values based on minimum threshold
threshold_filter <- function(row, threshold) {
  ifelse(row <= threshold, 0, row)
}

# apply the threshold filter to each row using mapply()
df_filtered <- mapply(threshold_filter, df_subset, df$NCmax)

# set the row names to NULL to remove the species names
rownames(df_filtered) <- NULL

# create a new csv file
write.csv(df_filtered, "moore_reads_maxCont_phyloseq.csv")
                                       
## Merge PCR/technical reps
df <- read.csv("moore_reads_maxCont_phyloseq.csv")

# edit the colnames
colnames(df) <- c("otuID","rep_1","rep_2","rep_3","rep1_1","rep1_2","rep1_3",
                  "rep2_1","rep2_2","rep2_3","rep3_1","rep3_2","rep3_3",
                  "rep4_1","rep4_2","rep4_3","rep5_1","rep5_2","rep5_3",
                  "rep6_1","rep6_2","rep6_3","rep7_1","rep7_2","rep7_3",
                  "rep8_1","rep8_2","rep8_3","rep9_1","rep9_2","rep9_3",
                  "rep10_1","rep10_2","rep10_3","rep11_1","rep11_2","rep11_3",
                  "rep12_1","rep12_2","rep12_3","rep13_1","rep13_2","rep13_3",
                  "rep14_1","rep14_2","rep14_3","rep15_1","rep15_2","rep15_3",
                  "rep16_1","rep16_2","rep16_3","rep17_1","rep17_2","rep17_3",
                  "rep18_1","rep18_2","rep18_3","rep19_1","rep19_2","rep19_3",
                  "rep20_1","rep20_2","rep20_3","rep21_1","rep21_2","rep21_3",
                  "rep22_1","rep22_2","rep22_3","rep23_1","rep23_2","rep23_3",
                  "rep24_1","rep24_2","rep24_3","rep25_1","rep25_2","rep25_3",
                  "rep26_1","rep26_2","rep26_3","rep27_1","rep27_2","rep27_3",
                  "rep28_1","rep28_2","rep28_3","rep29_1","rep29_2","rep29_3",
                  "EC","MC","OC")

# Get the unique column names that represent the same variable
cols <- unique(gsub("\\d+$", "", names(df)[grepl("\\d+$", names(df))]))

# Calculate the sum of replicates for each unique column name
for (col in cols) {
  replicate_cols <- grep(paste0("^", col, "\\d+$"), names(df))
  has_two_positives <- function(x) {
    sum(x > 0) >= 2 & sum(x != 0) >= 1
  }
  sum_replicates <- apply(df[, replicate_cols], 1, function(x) ifelse(has_two_positives(x), sum(x), 0))
  df[paste0("sum_", col)] <- sum_replicates
}

# Print the result
df

# Save result
write.csv(df,"moore_reads_pcrMerge_phyloseq.csv")

## Merge biological replicates
library(dplyr)

df <- read.csv("moore_reads_pcrMerge_collapsed_phyloseq.csv")

# edit the colnames
colnames(df) <- c("otuID","sample1_1","sample1_2","sample2_1","sample2_2","sample3_1","sample3_2",
                  "sample4_1","sample4_2","sample5_1","sample5_2",
                  "sample6_1","sample6_2","sample7_1","sample7_2",
                  "sample8_1","sample8_2","sample9_1","sample9_2",
                  "sample10_1","sample10_2","sample11_1","sample11_2",
                  "sample12_1","sample12_2","sample13_1","sample13_2",
                  "sample14_1","sample14_2","sample15_1","sample15_2")

# Get the unique column names that represent the same variable
cols <- unique(gsub("\\d+$", "", names(df)[grepl("\\d+$", names(df))]))

# Calculate the sum of replicates for each unique column name
for (col in cols) {
  replicate_cols <- grep(paste0("^", col, "\\d+$"), names(df))
  sum_replicates <- apply(df[, replicate_cols], 1, sum)
  df[paste0("sum_", col)] <- sum_replicates
}


# Print the result
df

# Save result
write.csv(df,"moore_reads_sampleMerge_phyloseq.csv")

## Convert to count data

### Load data
df1 <- read.csv("moore_reads_pcrMerge_collapsed_phyloseq.csv")
df2 <- read.csv("moore_reads_sampleMerge_phyloseq.csv")

### Collapse biological replicates & sum values
df1_merge <- df1 %>% mutate_if(is.numeric, ~1 * (. > 1))
df2_merge <- df2 %>% mutate_if(is.numeric, ~1 * (. > 1))

### Save result
write_csv(df1_merge, "moore_counts_pcrMerge_phyloseq.csv")
write_csv(df2_merge, "moore_counts_sampleMerge_phyloseq.csv")

