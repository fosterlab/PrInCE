# R script used to generate the list of fitted Gaussian mixture models 
# corresponding to the `scott` dataset bundled with the PrInCE package.
setwd("~/git/PrInCE-R")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(openxlsx)
library(PrInCE)

# load `scott` dataset
load("data/scott.rda", verbose = T)

# fit Gaussians
scott_gaussians = build_gaussians(scott)

# save
devtools::use_data(scott_gaussians)
