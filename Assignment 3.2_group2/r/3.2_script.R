## Assignment 3-2 ----
## Mixtures of t-distributions   ----

## *** ----
## SETUP ----

{
  rm(list=ls())
  if(dev.cur() != 1){dev.off()}
  
  ## Libraries
  library(microbenchmark)
  library(ggplot2)
  library(profvis)
  library(Matrix)
  library(reshape2)
  #library(plyr)
  library(dplyr)
  library(Rcpp)
  
  # Source functions
  source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/EKSAMEN/2/r/2.2_functions.R")
}