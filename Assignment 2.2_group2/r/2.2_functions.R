
## Monte Carlo functions ====

## Generate matrix of simulated X data, version 1
## Columns correspond to sample paths.
## Each sample path has length n.
## Number of sample paths is m.
## rfunc is the random generating function, eg. runif, rnorm, etc.
## Sn is our h function. A different h function could be supplied.
x_mat_1 = function(n, m, rfunc = runif(n*m, -1.9, 2.0)){
  matrix(rfunc, n, m)
}

## Define function which generates our simulated sample paths.
## This function is specific to the assignment.
## x is a length n vector (also works for matrices with n elements).
Sn = function(x){30 + cumsum(x)}

## Generate matrix of simulated values of Sn, version 1: for()
## x_mat is the matrix of simulated X-values.
## Sn is a function, specific to this assignment.
Sn_mat_gen_1 <- function(x_mat, Sn, n, m) {
  Sn_mat <- x_mat ## Place holder
  for(i in 1:m){
    Sn_mat[ ,i] = Sn(x_mat[ ,i]) ## For each column, apply Sn to the n x-values
  }
  Sn_mat
}

## Generate matrix of simulated values of Sn, version 2: apply()
## x_mat is the matrix of simulated X-values.
## Sn is a function, specific to this assignment.
Sn_mat_gen_2 <- function(x_mat, Sn) {
  Sn_mat <- NULL
  ## For each column, apply Sn to the n x-values
  Sn_mat <- apply(x_mat, 2, Sn) ## 2 for columns
  Sn_mat
}

MC_integral <- function() {
  
}


## Calculate probability of default
## Idea: Count number of columns with at least one value =< 0

## Version 1: Ratio of simulations that fail to total num of sims
default_prob_1 = function(Sn_mat){
  m = ncol(Sn_mat)
  count = 0
  for(i in 1:m){
    count = count + (min(Sn_mat[ ,i]) <= 0)
  }
  count/m
}


## List MC results for simulation.
## x_mat, Sn_mat_gen and default_prob are functions.

MCI_results <- function(Sn_mat_gen, default_prob, seed_switch = FALSE, seed = 1) {
  if(seed_switch){set.seed(seed)} ## Set seed if seed_switch=TRUE
  Sn_mat <- Sn_mat_gen
  prob <- default_prob(Sn_mat)
  list("mu_hat" = mean(Sn_mat), "Sn_mat" = Sn_mat, "prob" = prob)
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
