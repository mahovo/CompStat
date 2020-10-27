
## EM function factory ----
## https://cswr.nrhstat.org/8-1-basic-properties.html
{
  EM <- function(
    par,
    x_sim,
    d = 0.8, 
    c = 0.1, 
    gamma0 = 0.01,
    nu = 5,
    maxit = 100,
    maxit_gd = 100,
    maxit_range_check = 20, ## Number of times to adjust step size while checking parameter ranges
    epsilon = 1e-6,
    optimizer = gd,
    eps_gd = 1e-16,
    cb = NULL) {
    
    n <- length(x_sim)

    ## t_dens() **********************************************
    ## t-distribution density function
    const1 <- gamma((nu + 1)/2)/(sqrt(nu * pi) * gamma(nu/2))
    const2 <- (nu + 1) / 2
    t_dens <- function(x_i, mu, sigma) {
      if(sigma <= 0)
        warning("sigma <= 0")
      (const1 * 
          (1 + (x_i - mu)^2/(nu * sigma^2))^(-const2)
      )/sigma
    }
    
    ## Q() ***************************************************
    Q <- function(par, p_hat, xx){
      p <- par[1]
      mu1 <- par[2]
      mu2 <- par[3]
      sigma1 <- par[4]
      sigma2 <- par[5]
      nu <- nu1 <- nu2 <- 5
      
      ## Check that params are in range
      if(p <= 0)
        warning("p <= 0 in Q()")
      if(p < 0 || p > 1)
        return(Inf)
      if(sigma1 <= 0 )
        warning("sigma1 <= 0  in Q()")
      if(sigma2 <= 0 )
        warning("sigma2 <= 0  in Q()")
      
      sum(p_hat * log(p * t_dens(xx, mu1, sigma1))
          + (1 - p_hat) * log((1 - p) * t_dens(xx, mu2, sigma2))
      )
    }
    
    ## grad_Q() ***********************************************
    ## Gradient of Q
    grad_Q <- function(par, p_hat, x){
      p <- par[1]
      mu1 <- par[2]
      mu2 <- par[3]
      sigma1 <- par[4]
      sigma2 <- par[5]
  
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
          2 * const2 * (x - mu2)^2) / (nu * sigma2^3 + sigma2 * (x - mu2)^2
          )
        - (1/sigma2)
      ))
      
      c(dp, dmu1, dmu2, dsigma1, dsigma2)
    }
    
    ## GD() *********************************************************
    ## Gradient descent
    GD <- function (par, p_hat, eps_gd) {
      H <- function(...) {-Q(...)}
      grad_H <- function(...) {-grad_Q(...)}
      for(i in 1:maxit_gd) {
        it <- i
        H_value <- H(par, p_hat, xx = x_sim)
        if(!is.null(cb)) cb()
        grad_H <- grad_H(par = par, p_hat = p_hat, x = x_sim)
        #if(!is.null(cb)) cb()
        h_prime <- sum(grad_H^2)
        #if(!is.null(cb)) cb()
        
        ## Convergence criterion based on gradient norm
        if(h_prime <= eps_gd) break
        gamma <- gamma0
        
        ## Proposed descent step
        par1 <- par - gamma * grad_H

        ## Check that p is in range.
        ## Else reduce stepsize until p is in range.
        if(par1[1] <= 0 || par1[1] >= 1) {
          for(i in 1:maxit_range_check) {
            it = i
            if (par1[1] <= 0 || par1[1] >= 1) {
              warning("p out of range in GD()")
              gamma = gamma/2
              par1 <- par - gamma * grad_H
            } else {
              break
            }
            if(it == maxit_range_check ) warning("Range check of p exceeded max number of iterations.")
          }
        }
 
        ## Backtracking while descent is insufficient
        while(H(par1, p_hat, xx = x_sim) > H_value - c * gamma * h_prime) {
          gamma <- d * gamma
          par1 <- par - gamma * grad_H
        }
      }
      if(it == maxit_gd)
        warning("Number of iterations > maxit_gd = ", maxit_gd)
      par <- par1
    }
    
    ## E_step() **********************************************
    ## Calculate conditional expectations
    E_step <- function(par) {
      mu1 <- par[2]
      mu2 <- par[3]
      sigma1 <- par[4]
      sigma2 <- par[5]
      a <- numeric(n)
      b <- numeric(n)
      a <- p * t_dens(x_sim, mu1, sigma1)
      b <- (1 - p) * t_dens(x_sim, mu2, sigma2)
      a / (a + b)
    }
    
    ## M_step() **********************************************
    ## Maximize Q
    M_step <- function(par, p_hat) {
      switch(
        optimizer,
        gd = GD(par, p_hat, eps_gd),
        sg = SG(par, p_hat, gamma = 1e-5, num_epochs = 40)
      )
    }
    for(i in 1:maxit) {
      it <- i
      par0 <- par
      p_hat <- E_step(par)
      par <- M_step(par, p_hat)
      sum_sq_diff <- sum((par - par0)^2)
      if(!is.null(cb)) cb()
      if(sum_sq_diff <= epsilon * (sum(par^2) + epsilon))
        break
    }
    print(paste("Number of iterations: ", it))
    if(it == maxit)
      warning("Number of iterations > maxit = ", maxit)
    par
  }
}


