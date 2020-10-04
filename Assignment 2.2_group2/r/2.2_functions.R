
## Monte Carlo functions ====

## Generate matrix of simulated X data, version 1
## Columns correspond to sample paths.
## num_steps is the length n of sample path.
## num_paths is the number of sample paths.
## rfunc is the random generating function, eg. runif, rnorm, etc.
## Different h functions can be supplied, such as Sn or default.
x_mat_1 = function(num_steps, num_paths, rfunc = runif(num_steps*num_paths, -1.9, 2.0)){
  matrix(rfunc, num_steps, num_paths, byrow = FALSE)
}

## Calculate value of Sn for a single path
## Define function which generates our simulated sample paths.
## This function is specific to the assignment.
## x is a length n vector (also works for matrices with n elements).
Sn = function(x){
  #n = length(x)
  #Sn <- 30 + cumsum(x)
  #list("Sn" = Sn[n], "payoff_vect" = Sn)
  start_cap + sum(x)
}

## h function for calculating probability of default
## Returns 1 if default occurs in a given path, 0 otherwise
default = function(x){
  any(30 + cumsum(x) <= 0)
}

default_if <- function(x)
{
  default <- 0
  path <- 30 + cumsum(x)
  if (any(path <= 0)) default <- 1
  return(default)
}



## Generate matrix of simulated values of h(x), version 1: for()
## x_mat is the matrix of simulated X-values.
## h is a function, h(x). Specifically in this assignment, h is Sn.
## num_paths is the number of simulated paths.
h_vect_gen_1 <- function(x_mat, h) {
  h_vect <- NULL
  num_paths <- ncol(x_mat) ## Each column is a sample path
  for(i in 1:num_paths){
    h_vect[i] = h(x_mat[ ,i]) ## For each column, apply h to the n x-values
  }
  h_vect
}

## Generate matrix of simulated values of h(x), version 2: apply()
## x_mat is the matrix of simulated X-values.
## h is a function, h(x). Specifically in this assignment, h(x) is Sn.
## num_paths is the number of simulated paths.
h_vect_gen_2 <- function(x_mat, h) {
  ## For each column, apply h to the n x-values
  h_vect <- apply(x_mat, 2, h) ## 2 for columns
  h_vect
}

## for() loop with pre-assigned vector
h_vect_gen_3 <- function(x_mat, h) {
  num_paths <- ncol(x_mat) ## Each column is a sample path
  h_vect <- numeric(num_paths)
  for(i in 1:num_paths){
    h_vect[i] = h(x_mat[ ,i]) ## For each column, apply h to the n x-values
  }
  h_vect
}

## MC integration.
## Outputs:
##      mu_hat, the value of the integral.
##      h-vector, vector of values of h for each path.
## h is a function to integrate.
## h_vect_gen is a function that generates a vector of h(x) values.
## rfunc is a function that generates random x-values from the distribution of x.
## Seed can be switched on or off (TRUE/FALSE).

MCI <- function(h, h_vect_gen, num_steps, num_paths, rfunc, sigma_switch = TRUE) {
  h_vect <- h_vect_gen(x_mat_1(num_steps, num_paths, rfunc), h)
  mu_hat <- cumsum(h_vect) / 1:num_paths
  if(sigma_switch){
    sigma_hat <- numeric(num_paths)
    dev <- numeric(num_paths)
    ci_l <- numeric(num_paths) ## c.i. lower
    ci_u <- numeric(num_paths) ## c.i. upper
    for(i in 1:num_paths){
      sigma_hat[i] <- sd(h_vect[1:i]) ## sd for paths 1 thru i
      dev[i] <- 1.96 * sigma_hat[i] / sqrt(i)
    }
    ci_l <- mu_hat - dev
    ci_u <- mu_hat + dev
  } else {sigma_hat <- NA; ci_l = NA; ci_u = NA}
  list("mu_hat" = mu_hat, "h_vect" = h_vect, "sigma_hat" = sigma_hat, "CI_lower" = ci_l, "CI_upper" = ci_u)
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

IS <- function(h, h_vect_gen, num_steps, num_paths, theta, a, b, sigma_switch = TRUE) {
  ## p_vals is a matrix of random probabilities
  p_vals <- matrix(
    #runif(num_steps * num_paths, 0.02066011, 1.0), num_steps, num_paths, byrow = FALSE ## 0.02066011 < p < 1, lower values of p produce out of range x values
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  xg_mat <- xg_gen(p_vals, theta, a, b) ## Quantile function for g density
  h_vect <- h_vect_gen(xg_mat, h)
  g_dens <- apply(xg_mat, 2, function(x) {gn_1(x, theta, a, b)})
  #g_dens <- apply(xg_mat, 2, function(x) {g(x, theta, a, b)})
  f_dens <- replicate(num_paths, (1/3.9)^num_steps)
  #f_dens <- replicate(num_paths, (1/3.9))
  w_star <-  f_dens / g_dens
  mu_hat <- cumsum(h_vect*w_star) / 1:num_paths ## Element wise matrix multiplication
  if(sigma_switch){
    sigma_hat <- numeric(num_paths)
    dev <- numeric(num_paths)
    ci_l <- numeric(num_paths) ## c.i. lower
    ci_u <- numeric(num_paths) ## c.i. upper
    for(i in 1:num_paths){
      sigma_hat[i] <- sd(h_vect[1:i]*w_star[1:i]) ## sd for paths 1 thru i
      dev[i] <- 1.96 * sigma_hat[i] / sqrt(i)
    }
    ci_l <- mu_hat - dev
    ci_u <- mu_hat + dev
  } else {sigma_hat <- NA; ci_l = NA; ci_u = NA}
  list("mu_hat" = mu_hat, "h_vect" = h_vect, "gx_mat" = xg_mat, "sigma_hat" = sigma_hat, "CI_lower" = ci_l, "CI_upper" = ci_u)
}



## phi
phi <- function(theta, a, b) {
  #integrate(function(z){exp(theta*z)}, a, b)
  (exp(b*theta) - exp(a*theta))/theta
}



## g for single x_ik
## Index i: Path
## Index k: Step
g <- function(x_ki, theta, a, b) {
  exp(theta * x_ki) / phi(theta, a, b)
}

## Generate random x matrix from g distribution
## p_vect is a vector of random probabilities (vals in [0,1])
# xg_gen <- function(p_vect, theta, a, b) {
#   log(-(exp(a*theta) - exp(b*theta)) * p_vect * theta)/theta ## Inverse distribution (quantile function)
# }


## Inverse distribution (quantile function) for marginal g distribution function
xg_gen <- function(p_mat, theta, a = -1.9, b = 2.0) {
  if(theta != 0.0) {
    (
      log(
        p_mat*(
          exp(b * theta)-exp(a * theta)
        ) + exp(a * theta)
      )
  ) / theta
  } else {matrix(runif(length(p_mat), a, b), dim(p_mat))} ## Exception if theta=0: Uniform distribution
}


## gn for x-vector, version 1: power
gn_1 <- function(smp_path, theta, a, b) {
  n = length(smp_path)
  phi <- (exp(b*theta) - exp(a*theta))/(theta)
  exp(theta * sum(smp_path)) / phi(theta, a, b)^n
}

## gn for x-vector, version 2: for()
gn_2 <- function(smp_path, theta, a, b) {
  g = 1
  for(i in seq_along(smp_path)) {
    g = g * exp(theta * smp_path[i]) / phi(theta, a, b)
  }
  g
}

## gn function calls the phi function
gn_3 <- function(x, theta = 0.4, phi_func = phi) {
  n <- length(x)
  res <- 1 / phi_func(theta, -1.9, 2)**n * exp(theta * sum(x))
  return(res)
}


