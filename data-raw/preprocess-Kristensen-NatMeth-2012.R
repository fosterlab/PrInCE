# R script used to generate the `kristensen` data matrix bundled with the 
# PrInCE package from supporting information files available online at the 
# Nature Methods website.
setwd("~/git/PrInCE-R")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(openxlsx)

# read first replicate
repl = read.xlsx("data-raw/nmeth.2131-S2.xlsx", startRow = 30)

# convert to matrix
mat = repl %>%
  drop_na(Uniprot) %>%
  dplyr::select(Uniprot, starts_with("Ratio.M/L")) %>%
  mutate(Uniprot = gsub(";.*$", "", Uniprot)) %>%
  column_to_rownames('Uniprot') %>%
  as.matrix() %>%
  set_colnames(paste0("SEC_", seq_len(ncol(.))))

# drop proteins never quantified
keep = rowSums(is.finite(mat)) > 0

# rename
kristensen = mat[keep, ]
kristensen = kristensen[order(rownames(kristensen)), ]

# save
usethis::use_data(kristensen, overwrite = TRUE)
