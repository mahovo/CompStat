## EM ----
## https://cswr.nrhstat.org/8-1-basic-properties.html






## EM
{
  EM <- function(
    par,
    x,
    d = 0.8, 
    c = 0.1, 
    gamma0 = 0.01,
    maxit = 1000,
    epsilon = 1e-6,
    cb = NULL) {
    
    nu <- 5
    n <- length(x)
    const <- gamma((nu + 1)/2)/(sqrt(nu * pi) * gamma(nu/2))
    
    ## t-distribution density function
    t_dens <- function(x_i, mu, sigma) {
      if(sigma < 0)
        warning("sigma < 0")
      (const * 
         (1 + (
           (x_i - mu)^2)/(nu * sigma^2)
         )^(-(nu + 1)/2)
      )/sigma
    }
    
    ## Q
    Q <- function(par, p_tilde){
      p <- par[1]
      mu1 <- par[2]
      mu2 <- par[3]
      sigma1 <- par[4]
      sigma2 <- par[5]
      nu <- nu1 <- nu2 <- 5
      if(p < 0)
        #warning("p < 0")
        print("p < 0")
      if(sigma1 < 0 )
        #warning("sigma1 < 0")
        print("sigma1 < 0")
      if(sigma2 < 0 )
        #warning("sigma2 < 0")
        print("sigma2 < 0")
      sum = 0
      for(i in seq_along(x)){
        sum = sum + p_tilde[i] * (
          log(p * t_dens(x[i], mu1, sigma1))
          + (1 - p_tilde[i]) * (
            log((1 - p) * t_dens(x[i], mu2, sigma2))
          )
        )
      }
      - sum ## Negative log-likelihood
    }
    
    grad_Q <- function(par, p_tilde){
      p <- par[1]
      mu1 <- par[2]
      mu2 <- par[3]
      sigma1 <- par[4]
      sigma2 <- par[5]
      nu <- nu1 <- nu2 <- 5
      
      const <- (nu + 1) / 2
      
      dp <- -sum((p_tilde - p) /(p * (1-0)))
      dmu1 <- -sum((2 * p_tilde * const)/(x - mu1))
      dmu2 <- -sum((2 * (1 - p_tilde) * const)/(x - mu2))
      dsigma1 <- -sum(p_tilde * ((2 * const - 1) / sigma1))
      dsigma2 <- -sum((1 - p_tilde) * ((2 * const - 1) / sigma2))
      
      c(dp, dmu1, dmu2, dsigma1, dsigma2) ## Negative log likelihood
    }
    
    ## Calculate conditional expectations
    EStep <- function(par) {
      p <- par[1]
      mu1 <- par[2]
      mu2 <- par[3]
      sigma1 <- par[4]
      sigma2 <- par[5]
      a <- numeric(n)
      b <- numeric(n)
      for(i in seq_along(x)) {
        a[i] <- p * t_dens(x[i], mu1, sigma1)
        b[i] <- (1 - p)* t_dens(x[i], mu2, sigma2)
      }
      print(paste("EStep sigma1 = ", sigma1))
      print(paste("EStep sigma2 = ", sigma2))
      b / (a + b)
    }
    ## Maximize Q
    MStep <- function(par, p_tilde) {
      value <- Q(par, p_tilde)
      print(paste("MStep value = ", value))
      grad <- grad_Q(par, p_tilde)
      print(paste("min(grad) = ", min(grad)))
      print(paste("max(grad) = ", max(grad)))
      h_prime <- sum(grad^2)
      print(paste("h_prime = ", h_prime))
      ## Convergence criterion based on gradient norm
      if(h_prime <= epsilon) break
      gamma <- gamma0
      ## Proposed descent step
      par1 <- par - gamma * grad
      print(paste("par1 = ", par1))
      ## Backtracking while descent is insufficient
      while(Q(par1, p_tilde) > value - c * gamma * h_prime) {
        gamma <- d * gamma
        par1 <- par - gamma * grad
      }
      par <- par1
    }
    for(i in 1:maxit) {
      par0 <- par
      p_tilde <- EStep(par)
      par <- MStep(par, p_tilde)
      diff <- par - par0
      sum_sq_diff <- sum((par - par0)^2)
      if(!is.null(cb)) cb()
      if(sum_sq_diff <= epsilon * (sum(par^2) + epsilon))
        break
    }
    print(paste("Number of iterations: ", i))
    if(i == maxit)
      warning("Maximal number, ", maxit, ", of iterations reached")
    par
  }
}


## Gradient Descent ----
## https://cswr.nrhstat.org/7-2-descent-direction-algorithms.html
# GD <- function(
#   par, 
#   H,
#   gradient,
#   d = 0.8, 
#   c = 0.1, 
#   gamma0 = 0.01, 
#   epsilon = 1e-4, 
#   maxit = 1000,
#   cb = NULL
# ) {
#   for(i in 1:maxit) {
#     value <- H(par)
#     grad <- gradient(par)
#     h_prime <- sum(grad^2)
#     if(!is.null(cb)) cb()
#     # Convergence criterion based on gradient norm
#     if(h_prime <= epsilon) break
#     gamma <- gamma0
#     # Proposed descent step
#     par1 <- par - gamma * grad
#     # Backtracking while descent is insufficient
#     while(H(par1) > value - c * gamma * h_prime) {
#       gamma <- d * gamma
#       par1 <- par - gamma * grad
#     }
#     par <- par1
#   }
#   if(i == maxit)
#     warning("Maximal number, ", maxit, ", of iterations reached")
#   par
# }


