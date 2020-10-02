
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
  min(30 + cumsum(x)) <= 0 
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


## Calculate probability of default
## Idea: Count number of columns with at least one value =< 0

## Version 1: Ratio of simulations that fail to total num of sims
# default_prob_1 = function(h_mat){
#   num_paths = ncol(h_mat)
#   count = 0
#   for(i in 1:num_paths){
#     count = count + (min(h_mat[ ,i]) <= 0)
#   }
#   count/num_paths
# }


## MC integration.
## Outputs:
##      mu_hat, the value of the integral.
##      h-matrix, which can be used as input for default_prob.
## h is a function to integrate.
## h_mat_gen is a function that generates a matrix of h(x) values.
## rfunc is a function that generates random x-values from the distribution of x.
## Seed can be switched on or off (TRUE/FALSE).

MCI <- function(h, h_vect_gen, num_steps, num_paths, rfunc) {
  h_vect <- h_vect_gen(x_mat_1(num_steps, num_paths, rfunc), h)
  list("mu_hat" = mean(h_vect), "h_vect" = h_vect)
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
  #x_mat <- x_mat_1(num_steps, num_paths, runif(num_steps*num_paths, a, b))
  ## g_mat is g(x) matrix for num_paths paths
  #g_mat <- apply(x_mat, 2, function(x_mat){g(x_mat, theta, a, b)}) ## 2 for columns
  ## p_vals is a matrix of random probabilities
  ## 0.02066011 < p < 1
  ## Lower values of p produce out of range x values
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.02066011, 1.0), num_steps, num_paths, byrow = FALSE
  )
  xg_mat <- xg_gen(p_vals, theta, a, b) ## Quantile function for g density
  #g_mat <- apply(xg_mat, 2, function(xg_mat){g(xg_mat, theta, a, b)}) ## 2 for columns
  h_vect <- h_vect_gen(xg_mat, h)
  
  g_dens <- apply(xg_mat, 2, function(x) {gn_1(x, theta, a, b)})
  
  f_dens <- replicate(num_paths, (1/3.9)^num_steps)
  
  
  #w_star <- dunif(xg_mat, a, b)/g_mat ## dunif(x, -1.9, 2.0) = 0.2564103 for all x
  w_star <-  f_dens / g_dens
  
  mu_hat <- cumsum(h_vect*w_star) / 1:num_paths ## Element wise matrix multiplication
  if(sigma_switch){
    # *** missing ***
  } else {sigma_hat <- NA}
  list("mu_hat" = mu_hat, "h_vect" = h_vect, "gx_mat" = xg_mat, "sigma_hat" = sigma_hat, test = 1)
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


## Inverse distribution (quantile function)
xg_gen <- function(p_mat, theta, a = -1.9, b = 2.0) {
  (
    log(
      p_mat*(
        exp(b * theta)-exp(a * theta)
      ) + exp(a * theta)
    )
  ) / theta
}


## g for x-vector, version 1: power
gn_1 <- function(smp_path, theta, a, b) {
  n = length(smp_path)
  phi <- (exp(b*theta) - exp(a*theta))/(theta)
  exp(theta * sum(smp_path)) / phi(theta, a, b)^n
}

## g function calls the phi function
#g <- function(x, theta = 0.4, phi_func = phi) {
#  n <- length(x)
#  res <- 1 / phi_func(theta, -1.9, 2)**n * exp(theta * sum(x))
#  return(res)
#}


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


