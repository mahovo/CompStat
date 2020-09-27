## Assignment 2-2 ----
## Importance sampling: Univariate Simulation ----

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
}

## Constants

m = 1000 ## Number of paths
n = 100 ## Number of steps in path


## *** ----
## ESTIMATION ====


## p(100) ----

mc_results <- MC_ruin(
  x_mat_1(n, m), 
  Sn_mat_gen_2(x_mat_1(n, m), Sn), 
  default_prob_1
)
cat(paste0(
  "mu_hat: ", mc_results$mu_hat,
  "\nProbability of default in ", n, " steps: ", mc_results$prob
))


## Monte Carlo Integral ----


## Importance Sampling ----

## **************************************

## OOP ====  
## CONSTRUCTORS ----



## VALIDATORS ----



## HELPERS ----




## *** ----
## PROFILING ----


## *** ----
## BENCHMARKING ----

## Tests ----
## 1) Runif ----

## 2) Runtime wrt. n ----