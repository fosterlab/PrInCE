# R script used to generate the data matrix bundled with the PrInCE package from 
# supporting information files available online at the Molecular Systems
# Biology website.
setwd("~/git/PrInCE-R")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(openxlsx)

# read SEC (cytoplasmic) profiles
sec = read.xlsx("data-raw/Table EV11.xlsx", sheet = 1, startRow = 2)

# write heavy condition, replicate #1 
# (the smallest matrix after compression)
repl = sec %>%
  filter(`PCP-SILAC.Replicate` == 1) 
groups = repl$Majority.protein.IDs 
heavy = repl %>%
  select(starts_with("Ratio.H/L")) %>%
  mutate_all(as.numeric) %>%
  as.matrix() %>%
  set_rownames(groups) %>%
  set_colnames(paste0("SEC_", seq_len(ncol(.))))

# keep major protein group
rownames(heavy) = gsub(";.*$", "", rownames(heavy))

# rename
scott = heavy

# save
devtools::use_data(scott)
