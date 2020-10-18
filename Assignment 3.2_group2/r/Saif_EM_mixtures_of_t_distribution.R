
rm(list = ls())

## tracer function from Niels' package
#install.packages("remotes")  ## if 'remotes' is not installed
#library(remotes)
#install_github("nielsrhansen/CSwR/CSwR_package")


library(CSwR)
library(ggplot2)

# Setting starting variables for simulating variables
sigma1 <- 1
sigma2 <- 2
mu1 <- -4
mu2 <- 4
nu1 = nu2 <- 5 # Shape parameter
p <- 0.5
n <- 5000
initialPoints <- c(0.9, -0.9, 18)
set.seed(42)

# Simulating categorical variables
z = sample( c(TRUE, FALSE), n, replace = TRUE, prob = c(p, 1 - p) )

# Vector needed for simulated data from distribution
x <- numeric(n)
n1 <- sum(z)
x[z] <- rt(n1, nu1) * sigma1 + mu1
x[!z] <- rt(n - n1, nu2) * sigma2 + mu2


# Histogram
xDf <- as.data.frame(x)

p0 <- ggplot(data = xDf, aes(x = x)) + 
  geom_histogram(aes(y = ..density..), alpha = 0.5, color = "black") + 
  geom_density(color = 'red')

p0 + ggtitle("Histogram of simulated data")

#######################################################################################
########                        Function implementations                       ########
#######################################################################################

### t dens function ###
tdens <- function(x, nu, mu = 0, sigma = 1) {
  tmp1 <- (gamma((nu + 1) / 2)) / (sqrt(pi * nu * sigma^2) * gamma(nu / 2))
  tmp2 <- (1 + ((x - mu)^2 / (nu * sigma**2)))**(-(nu + 1) / 2)
  res <- tmp1 * tmp2
  return(res)
}


### Log likelihood for t-distribution ###
loglik <- function (par, x) {
  p <- par [1]
  
  if ( p < 0 || p > 1) {
    return( Inf )
  }
  mu1 <- par[2]
  mu2 <- par[3]
  
  -sum(log( 
            p * tdens(x, nu = nu1, mu = mu1, sigma = sigma1) + 
            (1-p) * tdens(x, nu = nu2, mu = mu2, sigma = sigma2) 
            )
       )
}


### EM implementation ###
EM_t_mix <- function(x) {
  n <- length(x)
  
  EStep <- function(par) {
    p <- par[1]
    mu1 <- par[2]
    mu2 <- par[3]
    
    a <- p * tdens(x, nu1, mu = mu1, sigma = sigma1)
    b <- (1-p) * tdens(x, nu2, mu = mu2, sigma = sigma2)
    a / (a + b)
  }
  
  MStep <- function(p_hat) {
    N1 <- sum(p_hat)
    N2 <- n - N1
    c(N1 / n, sum(p_hat * x) / N1, sum((1 - p_hat) * x) / N2)
  }
  
  function(par, epsilon = 1e-12, maxit = 50, cb = NULL) {
    for(i in 1:maxit) {
      par0 <- par
      par <- MStep(EStep(par))
      if(!is.null(cb)) cb()
      if(sum((par - par0)^2) <= epsilon * (sum(par^2) + epsilon)) break
    } 
    par 
  }
}


### Extracting function from EM algorithm
EStep <- environment(EM)$EStep ## E step
MStep <- environment(EM)$MStep ## M step

Q <- function(par, par_prime) {
  phat <- EStep(par_prime)
  p <- par[1]
  mu1 <- par[2]
  mu2 <- par[3]
  sum(phat * log(p * tdens(x, nu = nu1, mu = mu1, sigma = sigma1))  + 
        (1 - phat) * (log((1 - p) * tdens(x, nu = nu2, mu = mu2, sigma = sigma2))))
}

### Gradient descent ###
GD <- function(
  par, 
  H,
  gr,
  d = 0.8, # Constant used to find optimal step length, gamma, by backtracking
  c = 0.1, # Curvature condition that ensures that we move away from zero sufficiently
  # Also called Sufficient descent condition
  gamma0 = 0.01, # Step length
  epsilon = 1e-4, 
  maxit = 10000,
  cb = NULL
) {
  for(i in 1:maxit) {
    value <- H(par)
    grad <- gr(par)
    h_prime <- sum(grad^2)
    if(!is.null(cb)) cb()
    # Convergence criterion based on gradient norm
    if(h_prime <= epsilon) break
    gamma <- gamma0
    # Proposed descent step
    par1 <- par - gamma * grad
    # Backtracking while descent is insufficient
    while(H(par1) > value - c * gamma * h_prime) {
      gamma <- d * gamma
      par1 <- par - gamma * grad
    }
    par <- par1
  }
  if(i == maxit)
    warning("Maximal number, ", maxit, ", of iterations reached")
  cat("Converged after", i, "iterations! \n")
  par
}

### Newton implementation ###
newton <- function(
  par, 
  H,
  gr,
  hess,
  d = 0.8, 
  c = 0.1, 
  gamma0 = 1, 
  epsilon = 1e-10, 
  maxit = 1000,
  cb = NULL
) {
  for(i in 1:maxit) {
    value <- H(par)
    grad <- gr(par)
    if(!is.null(cb)) cb()
    if(sum(grad^2) <= epsilon) break
    Hessian <- hess(par)
    rho <- - drop(solve(Hessian, grad)) 
    gamma <- gamma0
    par1 <- par + gamma * rho
    h_prime <- t(grad) %*% rho 
    while(H(par1) > value +  c * gamma * h_prime) { 
      gamma <- d * gamma 
      par1 <- par + gamma * rho
    }
    par <- par1 
  }
  if(i == maxit)
    warning("Maximal number, ", maxit, ", of iterations reached")
  cat("Convergence after", i, "iterations! \n")
  par
}

## Tracer function
tracerFunction <- function(traceVariables) {
  givenTracer <- tracer(
    traceVariables, 
    N = 0, 
    expr = quote(
      {
        loglik <- loglik(par, x) 
        h_prime <- sum(gradQ(par)^2)
      }
    )
  )
}



### Numerical computation of gradient and hessian
library(numDeriv)
gradLike <- function(par) grad(function(par) loglik(par, x = x), par) ## Gradient for log likelihood
hessLike <- function(par) hessian(function(par) loglik(par, x = x), par) ## Hessian for log likelihood
gradQ <- function(par) -grad(Q, par, par_prime = par) ## Gradient for Q function


#######################################################################################
########                                EM Algorithm                           ########
#######################################################################################

## Testing log likelihood and EM algorithm
thetaHatOptim <- as.vector(unlist(optim(initialPoints, loglik, x = x)[1])) ## Convert to a vector

## Tracer and EM algorithm
EM <- EM_t_mix(x)
EM_tracer <- tracerFunction(c("par0", "par", "loglik", "h_prime"))
thetaHatEM <- EM(initialPoints, cb = EM_tracer$tracer)
EM_trace <- summary(EM_tracer)

## Plot of convergence
autoplot(EM_trace, y = loglik - min(loglik))
autoplot(EM_trace, y = h_prime)

# h_prime is 2norm of the gradient - information for myself #
## Print values
thetaHatOptim
thetaHatEM

## The gradients using Q function and the numDeriv package yields the more or less the same result
sum(gradLike(thetaHatEM)^2)
gradQ(thetaHatEM)

## Gradients for the two sets of optimal thetas
gradQ(thetaHatOptim)
gradQ(thetaHatEM)



#######################################################################################
########                              Gradient Descent                         ########
#######################################################################################

## log likelihood function which only takes parameters as argument used in GD
LL <- function(par) {
  loglik(par, x)
}

## Gradient descent and tracer
GD_tracer <- tracerFunction(c("par1", "par", "loglik", "h_prime"))
thetaHatGD <- GD(initialPoints, LL, gradLike, cb = GD_tracer$tracer)
GD_trace <- summary(GD_tracer)

autoplot(GD_trace, y = loglik - min(loglik))
autoplot(GD_trace, y = h_prime) # h_prime is 2norm of the gradient

#######################################################################################
########            Comparison of likelihood value for GD and EM theta         ########
#######################################################################################

# Compare log likelihood values for GD and EM
loglik(thetaHatGD, x = x) ## GD
loglik(thetaHatEM, x = x) ## Q function


gradQ(thetaHatGD) ## Gradient for GD parameters
gradQ(thetaHatEM) ## Gradient for EM parameters


testInitials <- c(0.8, 11, 19)


#######################################################################################
########                            Fisher information                         ########
#######################################################################################

# Q function which only takes the parameters as argument
Q2 <- function(par) Q(par, pHatEM)

## deriving hessian of Q using numDerive
iY <- hessian(Q2, pHatEM)

## Deriving the jacobian
Phi <- function(pp) MStep(EStep(pp))
DPhi <- jacobian(Phi, pHatEM)  ## Using numDeriv function 'jacobian'
iX <- (diag(1, 3) - t(DPhi)) %*% iY

## Deriving the inverse of the fisher information
iYinv <- solve(iY)
iYinv + iYinv %*% t(solve(diag(1, 3) - DPhi, DPhi))



#######################################################################################
########                               Newtons method                          ########
#######################################################################################

# Hessian of the log likelihood
hessLike <- function(par) hessian(function(par) loglik(par, x = x), par)

init <- c(0.9, -2, 4)

# Newton tracer
newton_tracer <- tracerFunction(c("par1", "par", "loglik", "h_prime"))
thetaHatNewton <- newton(init, LL, gradLike, hessLike, cb = newton_tracer$tracer)
newton_trace <- summary(newton_tracer)

autoplot(newton_trace, y = loglik - min(loglik))
autoplot(newton_trace, y = h_prime) # h_prime is 2norm of the gradient





newton(init, LL, gradLike, hessLike)
GD(init, LL, gradLike)
EM(init)




