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
  
  library(plotly)
  
  # Source functions
  source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 3.2_group2/r/3.2_functions.R")
}


# *** ----
## DATA ----
## > Simulate data ----

## >> Sample from t-distribution ----
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
  set.seed(42)
  p <- 0.5
  mu1 <- -10
  mu2 <- 10
  sigma1 <- 1
  sigma2 <- 2
  nu <- 5
  n <- 5000
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



## ESTIMATION ====

## > Simulate data ----

{
  set.seed(3630)
  true_params <- c(0.3, -4, 5, 1, 3)
  p <- true_params[1]
  mu1 <- true_params[2]
  mu2 <- true_params[3]
  sigma1 <- true_params[4]
  sigma2 <- true_params[5]
  nu <- 5
  n <- 5000
  z <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(p, 1 - p))
  ## Conditional simulation from the mixture components
  x_sim <- numeric(n)
  n1 <- sum(z)
  x_sim[z] <- rt(n1, df=5)*sigma1 + mu1 #rnorm(n1, mu1, sigma1)
  x_sim[!z] <- rt(n - n1, df=5)*sigma2 + mu2 #rnorm(n - n1, mu2, sigma2)
}

## Estimate ----
{
  set.seed(3630)
  
  ## Initial params: mu1 and mu2 should have same signs as real params (eg. guess from density plot) for convergence
  par <- c(runif(1, 0.1, 0.9), runif(1, -0.9, -0.1), runif(1, 0.1, 0.9), runif(1, 0.5, 2.0), runif(1, 0.5, 2.0))
  
  ## Tracer
  #EM_tracer = tracer(c("sum_sq_diff", "H_value", "h_prime", "gamma"), N = 1)
  #EM_tracer = tracer(c("sum_sq_diff"), N = 1)
  #EM_tracer = tracer(c("H_value"), N = 1)
  #EM_tracer = tracer(c("h_prime"), N = 1)
  #EM_tracer = tracer(c("grad_H"), N = 1)
  #EM_tracer = tracer(c("par1"), N = 1)
  EM_tracer = tracer(c("par", "sum_sq_diff", "H_value"), N = 1)
  
  estimate <- EM(
    par,
    x_sim,
    d = 0.7, 
    c = 0.3, 
    gamma0 = 0.01,
    nu = 5, 
    maxit = 1000,
    maxit_gd = 1, ## No need to optimize too much in GD. Even one step is enough.
    maxit_range_check = 1000, ## Drive p into range if needed
    epsilon = 1e-10,
    optimizer = "gd",
    eps_gd = 1e-16,
    #cb = NULL
    cb = EM_tracer$trace
  )
  #summary(EM_tracer)
  #head(summary(EM_tracer))
  #tail(summary(EM_tracer))
  print("Estimate: "); estimate
  print("True params: "); true_params
  print("Error: "); true_params - estimate
}

## > Plot estimated distribution ----
{
  hist(x_sim, prob=TRUE, breaks=40)
  rug(x_sim)
  x <- seq(min(y), max(y), length.out = 200)
  curve(f_h(x, estimate[4], estimate[5], estimate[2], estimate[3], estimate[1], nu), add = TRUE, col = "red")
}

## > Plot convergence ----

bind_rows(tracer = summary(EM_tracer)) %>% 
  autoplot(y = sum_sq_diff)

bind_rows(tracer = summary(EM_tracer)) %>% 
  autoplot(y = H_value)

bind_rows(tracer = summary(EM_tracer)) %>% 
  autoplot(y = H_value - min(H_value, na.rm = TRUE)) + scale_y_log10()

## same
# ggplot(summary(EM_tracer), aes(.time, H_value - min(H_value, na.rm = TRUE))) +
#   geom_point() + 
#   geom_line() + 
#   scale_y_log10()


## > Contour ----
## Note: The fixed parameter values are the theoretical values, not the empirical!

## >> Contour wrt p, mu1 ----
{
  mu2 <- true_params[3]
  sigma1 <- true_params[4]
  sigma2 <- true_params[5]
  x_contour = seq(0.0001, 0.999, length.out = 101)
  y_contour =  seq(-10, 10, length.out = 101)
  ## Define function to compute in each grid point
  ll_contour <- function(x, y) {
    negll(c(x, y, mu2, sigma1, sigma2), x_sim)
  }
  z_contour <- outer(x_contour, y_contour, Vectorize(ll_contour))
  
  fig3 <- plot_ly(
    x = x_contour, 
    y = y_contour, 
    z = t(z_contour), ## z matrix must be transposed for correct contour plot! (Don't know why)
    type = "contour",
    contours = list(
      colorscale='heatmap',
      showlabels = TRUE),
    ncontours = 30
  )
  x_label <- list(
    title = "p"
  )
  y_label <- list(
    title = "mu 1"
  )
  fig3 <- fig3 %>% 
    add_trace(type = "scatter", x = summary(EM_tracer)$par.1, y = summary(EM_tracer)$par.2, mode = "lines+markers") %>% 
    layout(xaxis = x_label, yaxis = y_label) %>% 
    layout(title=paste("mu2 = ", mu2, ", sigma1 = ", sigma1, ", sigma2 = ", sigma2))
  fig3
}


## >> Contour wrt sigma1, sigma2 ----
## Why does the algorithm move away from the minimum at first?
## The other params decreas sufficiently to drive the square d difference down.

{
  p <- true_params[1]
  mu1 <- true_params[2]
  mu2 <- true_params[3]
  x_contour = seq(0.1, 6.1, length.out = 101)
  y_contour =  seq(0.1, 6.1, length.out = 101)
  ## Define function to compute in each grid point
  ll_contour <- function(x, y) {
    negll(c(p, mu1, mu2, x, y), x_sim)
  }
  z_contour <- outer(x_contour, y_contour, Vectorize(ll_contour))
  
  fig4 <- plot_ly(
    x = x_contour, 
    y = y_contour, 
    z = t(z_contour), ## z matrix must be transposed for correct contour plot! (Don't know why)
    type = "contour",
    contours = list(
      showlabels = TRUE),
    ncontours = 300
  )
  x_label <- list(
    title = "sigma 1"
  )
  y_label <- list(
    title = "sigma 2"
  )
  fig4 <- fig4 %>% 
    add_trace(x = summary(EM_tracer)$par.4, y = summary(EM_tracer)$par.5, type = "scatter", mode = "lines+markers") %>% 
    layout(xaxis = x_label, yaxis = y_label) %>% 
    layout(title=paste("p2 = ", p, ", mu1 = ", mu1, ", mu2 = ", mu2))
  fig4
}

## z matrix must be transposed for correct contour plot! (Don't know why)
## Test:
{
  xtmp <- 1:5
  ytmp <- 1:5 * 10
  ztmp <- outer(xtmp, ytmp, function(x,y) {y - x})
  
  ## The contour line should be y - x, like here:
  plot_ly(
    x = xtmp, 
    y = ytmp, 
    z = t(ztmp), ## z matrix must be transposed for correct contour plot! (Don't know why)
    type = "contour",
    contours = list(
      showlabels = TRUE)
  )
}


