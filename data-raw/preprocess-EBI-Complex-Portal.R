# R script used to generate the dataset of known protein complexes
# bundled with the PrInCE package from raw data downloaded frm the 
# EBI Complex Portal website.
setwd("~/git/PrInCE-R")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# read data
dat = read.delim("data-raw/homo_sapiens.tsv.gz")

# keep complex names and components
complexes = dat %>%
  dplyr::select(Recommended.name, 
                Identifiers..and.stoichiometry..of.molecules.in.complex) %>%
  set_colnames(c("name", "subunits")) %>%
  mutate(subunits = strsplit(subunits, "\\|")) %>%
  unnest(subunits) %>%
  mutate(subunits = gsub("\\(.*$", "", subunits))

# rename, remove weird accessions, remove homodimers
gold_standard = complexes %>%
  filter(!grepl("CHEBI|-PRO_|_9606|CPX-", subunits)) %>%
  flavin::as_annotation_list("subunits", "name") %>%
  extract(lengths(.) > 1)

# save
devtools::use_data(gold_standard)
