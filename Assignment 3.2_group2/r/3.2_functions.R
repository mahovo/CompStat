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
    const1 <- gamma((nu + 1)/2)/(sqrt(nu * pi) * gamma(nu/2))
    
    ## t-distribution density function
    t_dens <- function(x_i, mu, sigma) {
      if(sigma <= 0)
        warning("sigma <= 0")
      (const1 * 
          (1 + (x_i - mu)^2/(nu * sigma^2))^(-(nu + 1)/2)
      )/sigma
    }
    
    ## Q
    Q <- function(par, p_hat){
      p <- par[1]
      mu1 <- par[2]
      mu2 <- par[3]
      sigma1 <- par[4]
      sigma2 <- par[5]
      nu <- nu1 <- nu2 <- 5
      if(p <= 0)
        warning("p <= 0")
        #print("p <= 0")
      if(sigma1 <= 0 )
        warning("sigma1 <= 0")
        #print("sigma1 <= 0")
      if(sigma2 <= 0 )
        warning("sigma2 <= 0")
        #print("sigma2 <= 0")
      sum = 0
      
      sum(p_hat * (
          log(p * t_dens(x, mu1, sigma1))
          + (1 - p_hat) * (
            log((1 - p) * t_dens(x, mu2, sigma2))
          )
        )
      )
      
      
      # for(i in seq_along(x)){
      #   sum = sum + p_hat[i] * (
      #     log(p * t_dens(x[i], mu1, sigma1))
      #     + (1 - p_hat[i]) * (
      #       log((1 - p) * t_dens(x[i], mu2, sigma2))
      #     )
      #   )
      # }
      # - sum ## Negative log-likelihood
    }
    
    grad_Q <- function(par, p_hat){
      p <- par[1]
      mu1 <- par[2]
      mu2 <- par[3]
      sigma1 <- par[4]
      sigma2 <- par[5]
      nu <- nu1 <- nu2 <- 5
      
      const2 <- (nu + 1) / 2
      
      dp <- sum((p_hat - p) /(p * (1 - p)))
      dmu1 <- sum((2 * p_hat * const2 * (x - mu1))/(nu * sigma1^2 + (x - mu1)^2))
      dmu2 <- sum((2 * (1 - p_hat) * const2 * (x - mu2))/(nu * sigma2^2 + (x - mu2)^2))
      dsigma1 <- sum(p_hat * (
        (
          2 * const2 * (x - mu1)^2) / (nu * sigma1^3 + sigma1 * (x - mu1)^2
        )
        - (1/sigma1)
      ))
      dsigma2 <- sum((1 - p_hat) * (
        (
          2 * const2 * (x - mu1)^2) / (nu * sigma1^3 + sigma1 * (x - mu1)^2
          )
        - (1/sigma1)
      ))
      
      c(dp, dmu1, dmu2, dsigma1, dsigma2)
    }
    
    ## Calculate conditional expectations
    EStep <- function(par) {
      mu1 <- par[2]
      mu2 <- par[3]
      sigma1 <- par[4]
      sigma2 <- par[5]
      a <- numeric(n)
      b <- numeric(n)
      # for(i in seq_along(x)) {
      #   a[i] <- p * t_dens(x[i], mu1, sigma1)
      #   b[i] <- (1 - p)* t_dens(x[i], mu2, sigma2)
      # }
      a <- p * t_dens(x, mu1, sigma1)
      b <- (1 - p) * t_dens(x, mu2, sigma2)
      #print(paste("EStep sigma1 = ", sigma1))
      #print(paste("EStep sigma2 = ", sigma2))
      a / (a + b)
    }
    ## Maximize Q
    MStep <- function(par, p_hat) {
      # N1 <- sum(p_hat)
      # N2 <- n - N1
      # c(N1 / n, sum(p_hat * x) / N1, sum((1 - p_hat) * x) / N2, 1, 2) ## Fixed sigma. 
      
      H_value <- -Q(par, p_hat)
      #print(paste("MStep value = ", H_value))
      grad_H <- - grad_Q(par, p_hat)
      #print(paste("min(grad) = ", min(grad_H)))
      #print(paste("max(grad) = ", max(grad_H)))
      h_prime <- sum(grad_H^2)
      #print(paste("h_prime = ", h_prime))
      ## Convergence criterion based on gradient norm
      if(h_prime <= epsilon) break
      gamma <- gamma0
      ## Proposed descent step
      par1 <- par - gamma * grad_H
      if (par1[1] <= 0 || par1[1] >= 1) {
        warning("p out of range")
        #print("p out of range")
      }
      if (par1[4] <= 0) {
        warning("sigma1 out of range")
      }
      if (par1[5] <= 0) {
        warning("sigma2 out of range")
      }
      #print(paste("par1 = ", par1))
      ## Backtracking while descent is insufficient
      while(Q(par1, p_hat) > H_value - c * gamma * h_prime) {
        gamma <- d * gamma
        par1 <- par - gamma * grad_H
      }
      par <- par1
    }
    for(i in 1:maxit) {
      par0 <- par
      p_hat <- EStep(par)
      par <- MStep(par, p_hat)
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

### Extracting function from EM algorithm
EStep <- environment(EM)$EStep ## E step
MStep <- environment(EM)$MStep ## M step



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


