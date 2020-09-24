# Monte Carlo St Version 1 ====
# Takes inputs:
# Sn (Sn function)
# x (function that generates vector of random data for single simulation)
# n (numeric, number of simulations)
mc_St_1 = function(St, x, n){
  if(!is.function(St)){stop("Sn is not a function.")}
  if(!is.function(x)){stop("x is not a function.")}
  if(!is.numeric(n)){stop("n is not a number.")}
  sim_mat = replicate(n, St(x())) # Matrix of sim data, cols are sample paths
  list("mu_hat" = mean(sim_mat), "sim_mat" = sim_mat)
}


# Monte Carlo St Version 2 ====
# Takes inputs:
# sim_mat (matrix, columns are sample paths)
# Sn (Sn function)

mc_St_2 = function(St, sim_mat){
  n = ncol(sim_mat) # number of sample paths
  for(i in 1:n){
    sim_mat[ ,i] = St(sim_mat[ ,i])
  }
  list("mu_hat" = mean(sim_mat), "sim_mat" = sim_mat)
}


# Monte Carlo Integral Version 1 ====

# h (h function, the function to integrate, outputs a sample path vector)
# x (vector of sampled values)

# h() is specific to the task, outputs a vector x = x_1, x_2, ..., x_n
h = function(x){
  n = ncol(x)
  St = function(x){30+cumsum(x)}
  sim_mat = apply(x, 2, St)
  evals = NULL # Vector of evaluated values
  for(i in 1:n){
    evals[i] = (min(sim_mat[ ,i]) <= 0)
  }
  evals
}
mc_integral_2 = function(x){
  nn = seq_along(x)
  muhat = cumsum(x) / nn;
  sigmahat = sd(x)
  list("mu_hat_vect" = muhat, "mu_hat" = muhat[length(muhat)], "sigma_hat" = sigmahat)
}



# Monte Carlo Integral Version 2 ====

# h (h function, the function to integrate, outputs a sample path vector)
# n (number of simulations)

# h() is specific to the task, outputs a vector x = x_1, x_2, ..., x_n
h = function(n){
  T_ = 100 # To enable mc_integral to take generic h-function, h doesn't
  # take T_ as input; only n.
  sim_Xn = function(T_, n){runif(T_*n, -1.9, 2)} # Vector of rand vals
  sim_mat = matrix(sim_Xn(T_, n), nrow = T_, byrow = T_) # As matrix
  St = function(sim_Xn){30+cumsum(sim_Xn)}
  for(i in 1:n){
    sim_mat[ ,i] = St(sim_mat[ ,i])
  }
  evals = NULL # Vector of evaluated values
  for(i in 1:n){
    evals[i] = (min(sim_mat[ ,i]) <= 0)
  }
  evals
}
mc_integral_1 = function(x){
  nn = seq_along(x)
  muhat = cumsum(x) / nn;
  sigmahat = sd(x)
  list("mu_hat" = muhat, "sigma_hat" = sigmahat)
}




# Probability ====
# Idea: Count number of columns with at least one value =< 0

# Ver 1: Ratio of simulations that fail to total num of sims
mc_prob_1 = function(spl_paths){
  n = ncol(spl_paths)
  count = 0
  for(i in 1:n){
    count = count + (min(spl_paths[ ,i]) <= 0)
  }
  count/n
}

# Ver 2: Ratio of num of elem <= 0 to tot num elem
mc_prob_2 = function(spl_paths){
  num_elem = length(spl_paths)
  count = sum(spl_paths <= 0)
  count/num_elem
}


# Helper functions (CSwR) ====
# Gamma simulation ----

# -> 4.3.1

rng_stream <- function(m, rng, ...) {
  args <- list(...)
  cache <- do.call(rng, c(m, args))
  j <- 0
  fact <- 1
  next_rn <- function(r = m) {
    j <<- j + 1
    if(j > m) {
      if(fact == 1 && r < m) fact <<- m / (m - r)
      m <<- floor(fact * (r + 1))
      cache <<- do.call(rng, c(m, args))
      j <<- 1
    }
    cache[j] 
  }
  next_rn
}

# -> 4.3.2

## r >= 1 
tfun <- function(y, a) {
  b <- 1 / (3 * sqrt(a))
  (y > -1/b) * a * (1 + b * y)^3  ## 0 when y <= -1/b
}

qfun <- function(y, r) {
  a <- r - 1/3
  tval <- tfun(y, a)
  exp(a * log(tval / a) - tval + a)
}

gammasim <- function(n, r, trace = FALSE) {
  count <- 0
  y <- numeric(n)
  y0 <- rng_stream(n, rnorm)
  u <- rng_stream(n, runif)
  for(i in 1:n) {
    reject <- TRUE
    while(reject) {
      count <- count + 1
      z <- y0(n - i)
      reject <- u(n - i) > qfun(z, r) * exp(z^2/2)
    }
    y[i] <- z
  }
  if(trace)
    cat("r =", r, ":", (count - n)/ count, "\n")  ## Rejection frequency
  tfun(y, r - 1/3)
}



