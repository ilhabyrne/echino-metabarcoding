
## Merge BOLD info to taxmat (23/08/23)

library(readr)
library(tidyverse)

allTax <- read_csv("taxonomy/moore_subsetTaxa_phyloseq.csv")
echinoTax <- read_csv("taxonomy/moore-reef_BOLD_SU.csv")

taxmat <- left_join(allTax, echinoTax,
                    by="otuID")

write_csv(taxmat, "taxonomy/moore_subsetTaxa_echinoBOLD_phyloseq.csv")