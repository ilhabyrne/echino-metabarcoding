## Load packages
library(plyr)
library(tidyr)
library(dplyr)

install.packages("tidyverse")
library(tidyverse)

## Max contamination filtering

### Load data
setwd("~/Desktop/Undergraduate/manuscript/Data_2023/longshore/phyloseq/")

df <- read.csv("long_reads_noSingletons_phyloseq.csv")

### set the row names as the otu names
rownames(df) <- df$uniq

### select all columns except for the minimum value column
df_subset <- df[, -ncol(df)]

### create a function to filter values based on minimum threshold
threshold_filter <- function(row, threshold) {
  ifelse(row <= threshold, 0, row)
}

### apply the threshold filter to each row using mapply()
df_filtered <- mapply(threshold_filter, df_subset, df$NCmax)

### set the row names to NULL to remove the species names
rownames(df_filtered) <- NULL

### create a new csv file
write.csv(df_filtered, "long_reads_maxCont_phyloseq.csv")


## Sample % threshold
### Remove any read counts within a sample that are less than a proportion of the total sample read count for that sample

### Load data
df <- read.csv("long_reads_maxCont_phyloseq.csv")

### Get the threshold value from the last row
threshold_row <- df[nrow(df),]

### Loop over the columns
for (col in 1:ncol(df)) {
  # Get the threshold value for the column
  threshold_value <- threshold_row[col]
  # Loop over the rows
  for (row in 1:(nrow(df)-1)) {
    # If the value is less than or equal to the threshold, set it to zero
    if (df[row, col] <= threshold_value) {
      df[row, col] <- 0
    }
  }
}

### Save the modified dataset to a new CSV file
write.csv(df, "long_reads_samplePercent_phyloseq.csv", row.names = FALSE)

## Merge biological replicates
library(dplyr)

df <- read.csv("long_reads_samplePercent_phyloseq.csv")

### edit the colnames
colnames(df) <- c("otuID","sample1_1","sample1_2","sample2_1","sample2_2","sample3_1","sample3_2",
                  "sample4_1","sample4_2","sample5_1","sample5_2",
                  "sample6_1","sample6_2","sample7_1","sample7_2",
                  "sample8_1","sample8_2","sample9_1","sample9_2",
                  "sample10_1","sample10_2","sample11_1","sample11_2",
                  "sample12_1","sample12_2","sample13_1","sample13_2",
                  "sample14_1","sample14_2","sample15_1","sample15_2",
                  "sample16_1","sample16_2","sample17_1","sample17_2",
                  "sample18_1","sample18_2","sample19_1","sample19_2",
                  "sample20_1","sample20_2","sample21_1","sample21_2",
                  "sample22_1","sample22_2","sample23_1","sample23_2",
                  "sample24_1","sample24_2","sample25_1","sample25_2",
                  "sample26_1","sample26_2","sample27_1","sample27_2",
                  "sample28_1","sample28_2","sample29_1","sample29_2",
                  "sample30_1","sample30_2","sample31_1","sample31_2",
                  "sample32_1","sample32_2","sample33_1","sample33_2",
                  "sample34_1","sample34_2",
                  "sample35_1","sample35_2","sample36_1","sample36_2",
                  "sample37_1","sample37_2","sample38_1","sample38_2",
                  "sample39_1","sample39_2","sample40_1","sample40_2",
                  "sample41_1","sample41_2","sample42_1","sample42_2",
                  "sample43_1","sample43_2","sample44_1","sample44_2",
                  "sample45_1","sample45_2","EC","MC","OC")

### Get the unique column names that represent the same variable
cols <- unique(gsub("\\d+$", "", names(df)[grepl("\\d+$", names(df))]))

### Calculate the sum of replicates for each unique column name
for (col in cols) {
  replicate_cols <- grep(paste0("^", col, "\\d+$"), names(df))
  sum_replicates <- apply(df[, replicate_cols], 1, sum)
  df[paste0("sum_", col)] <- sum_replicates
}


### Print the result
df

### Save result
write.csv(df,"long_reads_sampleMerge_phyloseq.csv")


## Convert to count data

### Load data
df1 <- read.csv("long_reads_percentThresh_phyloseq.csv")
df2 <- read.csv("long_reads_sampleMerge_phyloseq.csv")

### Collapse biological replicates & sum values
df1_merge <- df1 %>% mutate_if(is.numeric, ~1 * (. > 1))
df2_merge <- df2 %>% mutate_if(is.numeric, ~1 * (. > 1))

### Save result
write_csv(df1_merge, "long_counts_percentThresh_phyloseq.csv")
write_csv(df2_merge, "long_counts_sampleMerge_phyloseq.csv")
