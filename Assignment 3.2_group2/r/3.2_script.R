## Assignment 3-2 ----
## Mixtures of t-distributions   ----

## *** ----
## SETUP ----

{
  rm(list=ls())
  if(dev.cur() != 1){dev.off()}
  
  ## Libraries
  library(microbenchmark)
  library(ggplot2)
  library(profvis)
  library(Matrix)
  library(reshape2)
  #library(plyr)
  library(dplyr)
  library(Rcpp)
  
  #install.packages("remotes")  ## if 'remotes' is not installed
  #library(remotes)
  #install_github("nielsrhansen/CSwR/CSwR_package")
  library(CSwR)
  
  # Source functions
  source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/EKSAMEN/2/r/2.2_functions.R")
}


# *** ----
## DATA ----
## Simulate data ----

## Sample from t-distribution
## CSwR 4-2-1
{
  n <- 10000
  nu <- 5
  mu <- runif(1, -50, 50)
  sigma <- runif(1, 1, 20)
  #set.seed(2603)
  W <- rgamma(n, shape = nu/2, scale = 2) ## chi^2 distribution
  Z <- rnorm(n, 0, 1)
  T <- Z * sqrt(nu/W)
  X <- T * sigma + mu
  
  f_h <- function(x, sigma, mu, nu) {
    out <- numeric(length(x))
    const <- gamma((nu+1)/2)/(sqrt(nu*pi)*gamma(nu/2))
    for (i in seq_along(x)) {
      out[i] <- (const * (1 + ((x[i]-mu)^2)/(nu*sigma^2))^(-(nu+1)/2))/sigma
    }
    out
  }
  
  ## Check if it seems to match rt()s
  ## X <- T * sigma + mu
  par(mfrow=c(1,2))
  hist(X, prob=TRUE, breaks = 40)
  rug(X)
  x <- seq(min(X), max(X), length.out = 200)
  curve(f_h(x, sigma, mu, nu), add = TRUE, col = "red")
  ## x = rt(n, df=nu)*sigma + mu
  x = rt(n, df=nu)*sigma + mu
  hist(x, prob=TRUE, breaks = 40)
  rug(x)
  x <- seq(min(x), max(x), length.out = 200)
  curve(f_h(x, sigma, mu, nu), add = TRUE, col = "red")
}
par(mfrow=c(1,1))

## > t-distribution mixture ----
{
  p <- 0.25
  mu1 <- -3.5
  mu2 <- 4
  sigma1 <- 1
  sigma2 <- 2
  nu <- 5
  n <- 10000
  z <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(p, 1 - p))
  ## Conditional simulation from the mixture components
  y <- numeric(n)
  n1 <- sum(z)
  y[z] <- rt(n1, df=5)*sigma1 + mu1 #rnorm(n1, mu1, sigma1)
  y[!z] <- rt(n - n1, df=5)*sigma2 + mu2 #rnorm(n - n1, mu2, sigma2)
  
  f_h <- function(x, sigma1, sigma2, mu1, mu2, p, nu) {
    out <- numeric(length(x))
    const <- gamma((nu+1)/2)/(sqrt(nu*pi)*gamma(nu/2))
    for (i in seq_along(x)) {
      out[i] <- (p * const * (1 + ((x[i]-mu1)^2)/(nu*sigma1^2))^(-(nu+1)/2))/sigma1 + ((1-p) * const * (1 + ((x[i]-mu2)^2)/(nu*sigma2^2))^(-(nu+1)/2))/sigma2
    }
    out
  }
  
  hist(y, prob=TRUE, breaks=40)
  rug(y)
  x <- seq(min(y), max(y), length.out = 200)
  curve(f_h(x, sigma1, sigma2, mu1, mu2, p, nu), add = TRUE, col = "red")
}






{
  #set.seed(2603)
  par <- c(runif(5, 0.1, 1))
  EM_tracer = tracer(c('sum_sq_diff'), N = 1)
  EM(
    par,
    x,
    d = 0.8, 
    c = 0.1, 
    gamma0 = 0.01,
    maxit = 1000,
    epsilon = 1e-6,
    #cb = NULL
    cb = EM_tracer$trace
  )
  #summary(EM_tracer)
}





