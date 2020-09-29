
## Monte Carlo functions ====

## Generate matrix of simulated X data, version 1
## Columns correspond to sample paths.
## num_steps is the length n of sample path.
## num_paths is the number of sample paths.
## rfunc is the random generating function, eg. runif, rnorm, etc.
## Sn is our h function. A different h function could be supplied.
x_mat_1 = function(num_steps, num_paths, rfunc = runif(num_steps*num_paths, -1.9, 2.0)){
  matrix(rfunc, num_steps, num_paths)
}

## Define function which generates our simulated sample paths.
## This function is specific to the assignment.
## x is a length n vector (also works for matrices with n elements).
Sn = function(x){30 + cumsum(x)}

## Generate matrix of simulated values of h(x), version 1: for()
## x_mat is the matrix of simulated X-values.
## h is a function, h(x). Specifically in this assignment, h is Sn.
## num_paths is the number of simulated paths.
h_mat_gen_1 <- function(x_mat, h) {
  h_mat <- x_mat ## Place holder
  num_paths <- ncol(x_mat) ## Each column is a sample path
  for(i in 1:num_paths){
    h_mat[ ,i] = h(x_mat[ ,i]) ## For each column, apply h to the n x-values
  }
  h_mat
}

## Generate matrix of simulated values of h(x), version 2: apply()
## x_mat is the matrix of simulated X-values.
## h is a function, h(x). Specifically in this assignment, h(x) is Sn.
## num_paths is the number of simulated paths.
h_mat_gen_2 <- function(x_mat, h) {
  h_mat <- NULL
  ## For each column, apply Sn to the n x-values
  h_mat <- apply(x_mat, 2, h) ## 2 for columns
  h_mat
}

MC_integral <- function() {
  
}


## Calculate probability of default
## Idea: Count number of columns with at least one value =< 0

## Version 1: Ratio of simulations that fail to total num of sims
default_prob_1 = function(h_mat){
  num_paths = ncol(h_mat)
  count = 0
  for(i in 1:num_paths){
    count = count + (min(h_mat[ ,i]) <= 0)
  }
  count/num_paths
}


## List MC results for simulation.
## x_mat, Sn_mat_gen and default_prob are functions.

MCI_results <- function(h, h_mat_gen, num_steps, num_paths, rfunc, default_prob, seed_switch = FALSE, seed = 1) {
  if(seed_switch){set.seed(seed)} ## Set seed if seed_switch=TRUE
  h_mat <- h_mat_gen(x_mat_1(num_steps, num_paths, rfunc), h)
  prob <- default_prob(h_mat)
  list("mu_hat" = mean(h_mat), "h_mat" = h_mat, "prob" = prob)
}




MCI_ruin <- function(n, m, a, b) {
  # n is number of days the company earns/loses money
  # m is the number withdrawings
  # a and b are the upper and lower bound where the distribution takes place
  numberOfDefaults <- 0
  probability_of_default <- 0
  
  for(i in 1:m) {
    startCapital <- 30
    earnings <- runif(n, a, b)
    for(j in seq_along(earnings)) {
      startCapital <- startCapital + earnings[j]
      if(startCapital <= 0) {
        numberOfDefaults <- numberOfDefaults + 1
        break
      }
    }
  }
  probabilityOfDefault <- numberOfDefaults / m
  probabilityOfDefault
}


MCI_ruin_vectorized <- function(n, m, a, b) {
  # n is number of days the company earns/loses money
  # m is the number withdrawings
  # a and b are the upper and lower bound where the distribution takes place
  numberOfDefaults <- 0
  probabilityOfDefault <- 0
  
  for(i in 1:m) {
    startCapital <- 30
    earnings <- runif(n, a, b)
    capital <- startCapital + cumsum(earnings)
    if( min(capital) <= 0) {
      numberOfDefaults <- numberOfDefaults + 1
      next
    }
  }
  probabilityOfDefault <- numberOfDefaults / m
  probabilityOfDefault
}


## Importance Sampling functions ====

## g for single sample path, version 1: power
g_1 <- function(smp_path, theta, n) {
  exp(theta * sum(smp_path)) / phi(theta)^n
}

## g for single sample path, version 2: for()
g_2 <- function(smp_path, theta, n) {
  g = 1
  for(i in seq_along(smp_path)) {
    g = g *theta * smp_path[i] / phi(theta)
  }
  g
}


## phi
phi <- function(theta) {
  integrate(function(z){exp(theta*z)}, -1.9, 2)
}
