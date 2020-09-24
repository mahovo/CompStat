## Assignment 2-2 ----
## Importance sampling: Univariate Simulation ----

## *** ----
## SETUP ----

{
  rm(list=ls())
  if(dev.cur() != 1){dev.off()}
  
  ## Work. dir. and libraries
  library(microbenchmark)
  library(Rcpp)
  library(ggplot2)
  library(profvis)
  library(splines)
  library(glmnet)
  library(MASS)
  library(tidyr)
  library(dplyr)
  library(numDeriv)
  library(Matrix)
  library(tibble)
  library(reshape2)
  
  ## Source functions
  source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/Assignment 2.2_group2/r/2.2_functions_v01.R")
}
{
## save the workspace to the file .RData in the cwd 
# save.image()

## save the workspace to the file "workspace_HjOpg2.RData" in the cwd 
# save.image("workspace_HjOpg2.RData")

## Load a workspace into the current session.
## If you don't specify the path, the cwd is assumed.
# load("workspace_HjOpg2.RData") 
}




## *** ----
## ESTIMATION ----

## Test Monte Carlo average on gamma distributed data
## -> Simulated data (from CSwR 5.1.1) ----
x <- gammasim(1000, 8) 
nn <- seq_along(x)
muhat <- cumsum(x) / nn; sigmahat <- sd(x)
qplot(nn, muhat) + 
  geom_ribbon(aes(
    ymin = muhat - 1.96 * sigmahat / sqrt(nn), 
    ymax = muhat + 1.96 * sigmahat / sqrt(nn)), 
    fill = "gray") + 
  geom_line() + 
  geom_point()

## -> Gamma tails (from CSwR 5.1.3) ----
lambda <- 10

curve(pgamma(x, lambda, lower.tail = FALSE), 
      lambda, lambda + 30,
      ylab = "probability",
      main = "Gamma tail",
      ylim = c(0, 1))
curve(exp(lambda - x)*(x/lambda)^lambda, lambda, lambda + 30,
      add = TRUE, col = "red")
curve(lambda/(x-lambda)^2, lambda, lambda + 30,
      add = TRUE, col = "blue")

curve(pgamma(x, lambda, lower.tail = FALSE, log.p = TRUE), 
      lambda, lambda + 30,
      ylab = "log-probability",
      main = "Logarithm of Gamma tail",
      ylim = c(-20, 0))
curve(lambda - x + lambda*log(x/lambda), lambda, lambda + 30,
      add = TRUE, col = "red")
curve(-2 * log(x-lambda) + log(lambda), lambda, lambda + 30,
      add = TRUE, col = "blue")


## St Version 1 ====

tmax = 100
n = 1e5
x = function(){runif(tmax, -1.9, 2)} # Generate vector of simulated data for single sim path
St = function(x){30 + cumsum(x())}
mc_1 = mc_St_1(St, x, n)
mc_1$mu_hat
mc_prob_1(mc_1$sim_mat)
mc_prob_2(mc_1$sim_mat) # Don't use

## Plot St Version 1 simulations

tmax = 100
n = 1e3
x = function(){runif(T, -1.9, 2)} # Generate vector of simulated data for single sim path
St = function(x){30 + cumsum(x())}
mc_1 = mc_St_1(Sn, x, n)
mc_1$mu_hat

plot_data_1 = melt(mc_1$sim_mat) # From reshape2 library
plot_St_1_1000sim = ggplot() +
  geom_line(data = plot_data_1, aes(x = Var1, y = value, group = Var2, colour=factor(Var2))) +
  geom_hline(yintercept = mu_hat$mu_hat, color = "black") +
  labs(title = "St", subtitle = substitute(paste(n, " simulations, T = ", tmax, ". Black line is mu_hat."), list(n = n)), x = "t", y = "St") +
  theme(legend.position = "none")
plot_St_1_1000sim
{
  png('plot_St_1_1000sim.png', width=1800, height=1200, res=300)
  plot_St_1_1000sim
  dev.off()
}

## St Version 2 ====

tmax = 100
n = 1e5
sim_Xn = function(tmax, n){runif(T*n, -1.9, 2)} # Vector of rand vals
sim_mat = matrix(sim_Xn(tmax, n), nrow = tmax, byrow = T) # As matrix
St = function(x){30+cumsum(x)}

mc_2 = mc_Sn_2(Sn, sim_mat, n)
mc_2$mu_hat
mc_prob_1(mc_2$sim_mat)

## St Version 2 simulations

tmax = 100
n = 1e3
sim_Xn = function(tmax, n){runif(tmax*n, -1.9, 2)} # Vector of rand vals
sim_mat = matrix(sim_Xn(tmax, n), nrow = tmax, byrow = T) # As matrix
St = function(x){30+cumsum(x)}

mc_2 = mc_St_2(St, sim_mat)
mc_2$mu_hat
mc_prob_1(mc_2$sim_mat)

plot_data_2 = melt(mc_2$sim_mat) # From reshape2 library
plot_St_2_1000sim = ggplot() +
  geom_line(data = plot_data_2, aes(x = Var1, y = value, group = Var2, colour=factor(Var2))) +
  geom_hline(yintercept = mu_hat$mu_hat, color = "black") +
  labs(title = "St", subtitle = substitute(paste(n, " simulations, T = ", tmax, ". Black line is mu_hat."), list(n = n)), x = "t", y = "St") +
  theme(legend.position = "none")
plot_St_2_1000sim
{
  png('plot_St_1_1000sim.png', width=1800, height=1200, res=300)
  plot_St_2_1000sim
  dev.off()
}




## Monte Carlo Integral Version 1====

n = 1e5
mc_est_1 = mc_integral_1(h(n))
muhat = mc_est_1$mu_hat[length(mc_est_1$mu_hat)]; muhat
sigmahat = mc_est_1$sigma_hat; sigmahat
nn = seq_along(mc_est_1$mu_hat)

## Version 1, Plot 1
plot_mc_est_1 = qplot(seq_along(mc_est_1$mu_hat), mc_est_1$mu_hat) +
  # geom_ribbon(aes(
  #   ymin = mc_est_1$mu_hat - 1.96 * sigmahat / sqrt(nn),
  #   ymax = mc_est_1$mu_hat + 1.96 * sigmahat / sqrt(nn)),
  #   fill = "green") +
  geom_line() + 
  geom_point() +
  labs(title = "Monte Carlo estimate of p(n)", subtitle = paste("Max n = ", n), x = "n", y = "mu_hat")
plot_mc_est_1
{
  png('plot_mc_est_1.png', width=1800, height=1200, res=300)
  plot_mc_est_1
  dev.off()
}

## Version 1, Plot 2
## Running average
{
  n = 1e5
  cut = 1e4
  mc_est_2 = mc_integral_1(h(n))
  muhat = mc_est_2$mu_hat[length(mc_est_2$mu_hat)]; muhat
  sigmahat = mc_est_2$sigma_hat; sigmahat
  nn = seq_along(mc_est_2$mu_hat)[-(1:cut)]
  plot_mc_est_2 = qplot(seq_along(mc_est_2$mu_hat)[-(1:cut)], mc_est_2$mu_hat[-(1:cut)]) +
    geom_ribbon(aes(
      ymin = mc_est_2$mu_hat[-(1:cut)] - 1.96 * sigmahat / sqrt(nn),
      ymax = mc_est_2$mu_hat[-(1:cut)] + 1.96 * sigmahat / sqrt(nn)),
      fill = "green") +
    geom_line() + 
    geom_point() +
    labs(title = "Monte Carlo estimate of p(n)", subtitle = paste("Max n = ", n, ". First ", cut, " estimates removed."), x = "n", y = "mu_hat")
  plot_mc_est_2
}
{
  png('plot_mc_est_2.png', width=1800, height=1200, res=300)
  plot_mc_est_2
  dev.off()
}

## Monte Carlo Integral version 2 ====

## Version 2, Plot 1
{
  n = 1e5
  tmax = 100
  x = matrix(runif(tmax*n, -1.9, 2), tmax, n)
  mc_est_3 = mc_integral_2(h(x))
  muhat_3 = mc_est_3$mu_hat; muhat_3
  sigmahat_3 = mc_est_3$sigma_hat; sigmahat_3
  nn = seq_along(mc_est_3$mu_hat_vect)
  plot_mc_est_3a = qplot(seq_along(mc_est_3$mu_hat_vect), mc_est_3$mu_hat_vect) +
    ## geom_ribbon(aes(
    ##   ymin = mc_est_1$mu_hat - 1.96 * sigmahat / sqrt(nn),
    ##   ymax = mc_est_1$mu_hat + 1.96 * sigmahat / sqrt(nn)),
    ##   fill = "green") +
    geom_line() + 
    geom_point() +
    labs(title = "Monte Carlo estimate of p(t)", subtitle = paste("Max n = ", n), x = "n", y = "mu_hat")
  plot_mc_est_3a
}
{
  png('plot_mc_est_3a.png', width=1800, height=1200, res=300)
  plot_mc_est_3a
  dev.off()
}

## Version 2, Plot 2
{
  nn = seq_along(mc_est_3$mu_hat_vect)[-(1:cut)]
  plot_mc_est_3b = qplot(seq_along(mc_est_3$mu_hat_vect)[-(1:cut)], mc_est_3$mu_hat_vect[-(1:cut)]) +
    geom_ribbon(aes(
      ymin = mc_est_3$mu_hat_vect[-(1:cut)] - 1.96 * sigmahat / sqrt(nn),
      ymax = mc_est_3$mu_hat_vect[-(1:cut)] + 1.96 * sigmahat / sqrt(nn)),
      fill = "green") +
    geom_line() + 
    geom_point() +
    labs(title = "Monte Carlo estimate of p(t)", subtitle = paste("Max n = ", n, ". First ", cut, " estimates removed."), x = "n", y = "mu_hat")
  plot_mc_est_3b
}
{
  png('plot_mc_est_3b.png', width=1800, height=1200, res=300)
  plot_mc_est_3b
  dev.off()
}

## ====

## Get columns with elements <= 0
zero_col_IDs = which(apply(mc_1$sim_mat, 2, function(col) any(col <= 0)))
zero_cols = mc_1$sim_mat[ , zero_col_IDs]

## Count number of columns with an element <= 0 in column
sum(apply(mc_1$sim_mat, 2, function(col) any(col <= 0)))

## Count number of elements <= 0 in matrix
sum(mc_1$sim_mat <= 0)



## Importance Sampling ====





## *** ----
## FUNCTIONS in development ----
## See CS2_functions



## **************************************


## *** ----
## CONSTRUCTORS ----


## *** ----
## VALIDATORS ----


## *** ----
## HELPERS ----




## *** ----
## PROFILING ----


## *** ----
## BENCHMARKING ----

## Tests:
## 1) Runif
tmax = 100
n = 1e5
bench_runif = microbenchmark(
  matrix(runif(n*tmax, -1.9, 2), tmax, n),
  replicate(n, runif(tmax, -1.9, 2))
)
levels(bench_runif$expr) = c("B", "A")
plot_bench_runif = autoplot(bench_runif) +
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  labs(title = "Benchmark generation of sample matrix") +
  theme(legend.position = "none")
plot_bench_runif
{
  png('plot_bench_runif.png', width=1800, height=1200, res=300)
  plot_bench_runif
  dev.off()
}

## 2) Runtime wrt. n






n_bench = microbenchmark(
  {
    T = 100
    n = 1e1
    x = function(){runif(T, -1.9, 2)} # Generate vector of simulated data for single sim path
    h = function(x){30 + cumsum(x())}
    mc_1 = mc_Sn_1(h, x, n)
  },
  {
    T = 100
    n = 1e2
    sim_Xn = function(T, n){runif(T*n, -1.9, 2)} # Vector of rand vals
    spl_paths = matrix(sim_Xn(T, n), nrow = T, byrow = T) # As matrix
    h = function(x){30+cumsum(x)}
    mc_2 = mc_Sn_2(h, spl_paths, n)
  },
  {
    T = 100
    n = 1e3
    x = function(){runif(T, -1.9, 2)} # Generate vector of simulated data for single sim path
    h = function(x){30 + cumsum(x())}
    mc_1 = mc_Sn_1(h, x, n)
  },
  {
    T = 100
    n = 1e4
    x = function(){runif(T, -1.9, 2)} # Generate vector of simulated data for single sim path
    h = function(x){30 + cumsum(x())}
    mc_1 = mc_Sn_1(h, x, n)
  },
  {
    T = 100
    n = 1e5
    x = function(){runif(T, -1.9, 2)} # Generate vector of simulated data for single sim path
    h = function(x){30 + cumsum(x())}
    mc_1 = mc_Sn_1(h, x, n)
  }
)
levels(n_bench$expr) = c("v1", "v2", "v3", "v4", "v5")
{
  n_bench_df = data.frame(
    num_iter = 10^(1:5),
    runtime = summary(n_bench)$median
  )
}
{
  plot_n_bench<- ggplot(data = n_bench_df, aes(x = num_iter, y = runtime)) +
    geom_point(aes(x = num_iter, y = runtime), size = 3, col = "blue") + 
    geom_line(aes(x = num_iter, y = runtime), col = "blue") +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    labs(x = "Number of samples", y = "Median runtime (microseconds)") +
    labs(title = 'St() benchmark', subtitle = 'Runtime of St() wrt. number of samples')
  plot_n_bench
}
{
  png('plot_n_bench.png', width=1800, height=1200, res=300)
  plot_n_bench
  dev.off()
}
