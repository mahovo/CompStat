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

num_paths = 1000 ## Number of paths
num_steps = 100 ## Number of steps in path


## *** ----
## ESTIMATION ====


## > p(100) ----

mc_results <- MCI_results(
  h_mat_gen_1(x_mat_1(num_steps, num_paths), Sn, num_steps, num_paths), 
  default_prob_1,
  seed_switch = TRUE,
  seed = 123
)
cat(paste0(
  "mu_hat: ", mc_results$mu_hat,
  "\nProbability of default in ", num_steps, " steps: ", mc_results$prob
))

mc_results <- MCI(
  Sn, ## h(x) function
  h_mat_gen_2, ## Generate matrix of h(x) values
  num_steps,
  num_paths,
  runif(num_steps*num_paths, -1.9, 2.0),
  seed_switch = TRUE,
  seed = 123
)
cat(paste0(
  "mu_hat: ", mc_results$mu_hat,
  "\nProbability of default in ", num_steps, " steps: ", default_prob_1(mc_results$h_mat)
))




## > Monte Carlo Integral ----


## > Importance Sampling ----

## **************************************

## OOP ====  
## > CONSTRUCTORS ----



## > VALIDATORS ----



## > HELPERS ----




## *** ----
## PROFILING ----


## *** ----
## BENCHMARKING ----

## > TESTS ----
## >> Sn_mat_gen_1 vs Sn_mat_gen_2 ----

## >> Runif ----

## >> Runtime wrt. n ----