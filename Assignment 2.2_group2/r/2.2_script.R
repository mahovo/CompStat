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
## > Monte Carlo Integral ----

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



## > Importance Sampling ----

## Test that g is a density function
x <- x_mat_1(100, 1, rfunc = runif(num_steps*num_paths, -1.9, 2.0))
g_i(x, theta = 1, -1.9, 2)
phi(1, -1.9, 2)
tmp_int <- integrate(function(x){g_i(x, theta = 1, -1.9, 2)}, -1.9, 2.0)
str(tmp_int)
{
  theta_vals <- seq(0.5, 1.5, 0.01)
  a=-1.9
  b=2
  integrate(function(x){g_i(x, theta = 1, a, b)}, a, b)
  theta_test <- sapply(
    tmp_theta,
    function(test) {integrate(function(x){g_i(x, theta = test, a, b)}, a, b)$value}
  )
  plot_optimal_theta_sim = qplot(theta_vals, theta_test, xlab="theta", ylab="integral of g") + 
    geom_vline(xintercept = 1, color = "red", size = 0.5) +
    geom_hline(yintercept = 1, color = "red", size = 0.5) +
    #scale_x_continuous(trans = 'log10')+
    labs(title = "Optimal theta", subtitle = "For density integrating to 1")
  plot_optimal_theta_sim
}

## IS

{
#set.seed(123)
is_results <- IS(
  h = Sn, ## h(x) function
  h_mat_gen = h_mat_gen_2, ## Generate matrix of h(x) values
  num_steps = 100, 
  num_paths = 1000, 
  rfunc = runif(num_steps*num_paths, -1.9, 2.0),
  sigma_switch = FALSE,
  theta = 1,
  a = -1.9,
  b = 2.0
)
is_results$mu_hat
}
{
set.seed(123)
mc_results <- MCI(
  Sn, ## h(x) function
  h_mat_gen_2, ## Generate matrix of h(x) values
  num_steps,
  num_paths,
  runif(num_steps*num_paths, -1.9, 2.0)
)
mc_results$mu_hat
}


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