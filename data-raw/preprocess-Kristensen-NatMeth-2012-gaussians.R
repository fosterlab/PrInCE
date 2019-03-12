# R script used to generate the list of fitted Gaussian mixture models 
# corresponding to the `kristensen` dataset bundled with the PrInCE package.
setwd("~/git/PrInCE-R")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(openxlsx)
library(PrInCE)

# load `kristensen` dataset
load("data/kristensen.rda", verbose = TRUE)

# fit Gaussians
kristensen_gaussians = build_gaussians(kristensen)

# save
usethis::use_data(kristensen_gaussians)
