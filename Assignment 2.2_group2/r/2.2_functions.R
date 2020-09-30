
## Monte Carlo functions ====

## Generate matrix of simulated X data, version 1
## Columns correspond to sample paths.
## num_steps is the length n of sample path.
## num_paths is the number of sample paths.
## rfunc is the random generating function, eg. runif, rnorm, etc.
## Sn is our h function. A different h function could be supplied.
x_mat_1 = function(num_steps, num_paths, rfunc = runif(num_steps*num_paths, -1.9, 2.0)){
  matrix(rfunc, num_steps, num_paths, byrow = FALSE)
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


## MC integration.
## Outputs:
##      mu_hat, the value of the integral.
##      h-matrix, which can be used as input for default_prob.
## h is a function to integrate.
## h_mat_gen is a function that generates a matrix of h(x) values.
## rfunc is a function that generates random x-values from the distribution of x.
## Seed can be switched on or off (TRUE/FALSE).

MCI <- function(h, h_mat_gen, num_steps, num_paths, rfunc, seed_switch = FALSE, seed = 1) {
  if(seed_switch){set.seed(seed)} ## Set seed if seed_switch=TRUE
  h_mat <- h_mat_gen(x_mat_1(num_steps, num_paths, rfunc), h)
  list("mu_hat" = mean(h_mat), "h_mat" = h_mat)
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

## IS function

IS <- function(h, h_mat_gen, num_steps, num_paths, rfunc, seed_switch = FALSE, seed = 1, sigma_switch = TRUE, ...) {
  if(seed_switch){set.seed(seed)} ## Set seed if seed_switch=TRUE
  f_mat <- x_mat_1(num_steps, num_paths, rfunc)
  ## g_mat is g(x) matrix for num_paths paths
  #g_mat <- h_mat_gen_2(x_mat, function(x_mat){g(x_mat, theta, a, b)})
  g_mat <- apply(f_mat, 2, function(f_mat){g(f_mat, ...)}) ## 2 for columns
  h_mat <- h_mat_gen(g_mat, h)
  w_star <- f_mat/g_mat
  
  mu_hat <- mean(h_mat*w_star) ## Element wise matrix multiplication
  if(sigma_switch){
    sigma_hat <- 1
  } else {sigma_hat <- NA}
  list("mu_hat" = mu_hat, "h_mat" = h_mat, "sigma_hat" = sigma_hat)
}


## g for x-vector, version 1: power
gn_1 <- function(smp_path, theta, a, b) {
  n = length(smp_path)
  exp(theta * sum(smp_path)) / phi(theta, a, b)^n
}

## g for x-vector, version 2: for()
gn_2 <- function(smp_path, theta, a, b) {
  g = 1
  for(i in seq_along(smp_path)) {
    g = g * exp(theta * smp_path[i]) / phi(theta, a, b)
  }
  g
}

## g function calls the phi function
gn_3 <- function(x, theta = 0.4, phi_func = phi) {
  n <- length(x)
  res <- 1 / phi_func(theta, -1.9, 2)**n * exp(theta * sum(x))
  
  return(res)
}

## g for single x_ik
## Index i: Path
## Index k: Step
g <- function(x_ki, theta, a, b) {
  exp(theta * x_ki) / phi(theta, a, b)
}
str(g)

## phi
phi <- function(theta, a, b) {
  #integrate(function(z){exp(theta*z)}, a, b)
  exp(b*theta) - exp(a*theta)
}
