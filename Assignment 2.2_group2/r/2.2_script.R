## Assignment 2-2 ----
## Importance sampling: Univariate Simulation ----

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
  
  # Source functions
  source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/r/2.2_functions.R")
}

## Constants

num_paths = 1000 ## Number of paths
num_steps = 100 ## Number of steps in path


## *** ----
## ESTIMATION ====
## > Monte Carlo Integral ----

## > p(100) ----

mc_results <- MCI_results(
  h_mat_gen_1(x_mat_1(num_steps, num_paths), Sn, num_steps, num_paths), 
  default_prob_1,
  seed_switch = TRUE,
  seed = 123
)
cat(paste0(
  "mu_hat: ", mc_results$mu_hat,
  "\nProbability of default in ", num_steps, " steps: ", mc_results$prob
))

mc_results <- MCI(
  Sn, ## h(x) function
  h_mat_gen_2, ## Generate matrix of h(x) values
  num_steps,
  num_paths,
  runif(num_steps*num_paths, -1.9, 2.0),
  seed_switch = TRUE,
  seed = 123
)
cat(paste0(
  "mu_hat: ", mc_results$mu_hat,
  "\nProbability of default in ", num_steps, " steps: ", default_prob_1(mc_results$h_mat)
))



## > Importance Sampling ----

## Test that g is a density function
x <- x_mat_1(100, 1, rfunc = runif(num_steps*num_paths, -1.9, 2.0))
g_i(x, theta = 1, -1.9, 2)
phi(1, -1.9, 2)
tmp_int <- integrate(function(x){g_i(x, theta = 1, -1.9, 2)}, -1.9, 2.0)
str(tmp_int)
{
  theta_vals <- seq(0.5, 1.5, 0.01)
  a=-1.9
  b=2
  integrate(function(x){g_i(x, theta = 1, a, b)}, a, b)
  theta_test <- sapply(
    tmp_theta,
    function(test) {integrate(function(x){g_i(x, theta = test, a, b)}, a, b)$value}
  )
  plot_optimal_theta_sim = qplot(theta_vals, theta_test, xlab="theta", ylab="integral of g") + 
    geom_vline(xintercept = 1, color = "red", size = 0.5) +
    geom_hline(yintercept = 1, color = "red", size = 0.5) +
    #scale_x_continuous(trans = 'log10')+
    labs(title = "Optimal theta", subtitle = "For density integrating to 1")
  plot_optimal_theta_sim
}

## IS

## Simulating x values with g:

## Test that x values are in [-1.9, 2.0]
{
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  range(xg_gen(p_vals, theta=1, a=-1.9, b=2.0))
}

## g is a density (integrates to 1, except for theta = 0, where g is undefined)


## Compare IS and MC for Sn-function
num_paths = 1000 ## Number of paths
num_steps = 100 ## Number of steps in path
{
set.seed(123)
is_results <- IS(
  h = Sn, ## h(x) function
  #h = default,
  h_vect_gen = h_vect_gen_2, ## Generate matrix of h(x) values
  num_steps, 
  num_paths,
  theta = -0.2,
  a = -1.9,
  b = 2.0,
  sigma_switch = FALSE
)
is_results$mu_hat
}
{
## mu_hat is expected value of Sn
#set.seed(123)
mc_results_1 <- MCI(
  Sn, ## h(x) function
  h_vect_gen_2, ## Generate vector of h(x) values
  num_steps,
  num_paths,
  runif(num_steps*num_paths, -1.9, 2.0)
)
mc_results_1$mu_hat
}

## Compare IS and MC for default function
num_paths = 1000 ## Number of paths
num_steps = 100 ## Number of steps in path
{
  set.seed(42)
  is_results <- IS(
    h = default,
    h_vect_gen = h_vect_gen_2, ## Generate vector of h(x) values
    num_steps, 
    num_paths,
    theta = -0.2,
    a = -1.9,
    b = 2.0,
    sigma_switch = TRUE
  )
  cat("mu_hat = ", paste(is_results$mu_hat[num_paths], "\n"))
  cat("sigma_hat = ", is_results$sigma_hat[num_paths], "\n")
  cat("c.i. = ", is_results$CI_lower[num_paths], "\n")
  cat("c.i. = ", is_results$CI_upper[num_paths], "\n")
}
{
  qplot(1:num_paths, is_results$mu_hat) + 
    geom_ribbon(
      mapping = aes(
        ymin = is_results$CI_lower, 
        ymax = is_results$CI_upper
      ), fill = "gray") +
    #coord_cartesian(ylim = c(min(is_results$CI_lower, na.rm = TRUE), max(is_results$CI_upper, na.rm = TRUE))) +
    coord_cartesian(ylim = c(-0.001, 0.006)) +
    geom_line() + 
    #geom_point() +
    labs(title = "Importance sampling", subtitle = paste(num_steps, " steps", num_paths, " paths"), x = "number of paths", y = "mu_hat")
}




{
  ## mu_hat = expected probability of default
  set.seed(42)
  mc_results_2 <- MCI(
    default, ## h(x) function
    h_vect_gen_2, ## Generate vector of h(x) values
    num_steps,
    num_paths,
    runif(num_steps*num_paths, -1.9, 2.0),
    sigma_switch = TRUE
  )
  mc_results_2$mu_hat
}
{
  qplot(1:num_paths, mc_results_2$mu_hat) + 
    geom_ribbon(
      mapping = aes(
        ymin = mc_results_2$CI_lower, 
        ymax = mc_results_2$CI_upper
      ), fill = "gray") +
    #coord_cartesian(ylim = c(min(is_results$CI_lower, na.rm = TRUE), max(is_results$CI_upper, na.rm = TRUE))) +
    coord_cartesian(ylim = c(-0.003, 0.009)) +
    geom_line() + 
    #geom_point() +
    labs(title = "Monte Carlo", subtitle = paste(num_steps, " steps", num_paths, " paths"), x = "number of paths", y = "mu_hat")
}





## **************************************

## OOP ====  
## > CONSTRUCTORS ----



## > VALIDATORS ----



## > HELPERS ----




## *** ----
## PROFILING ----

source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/r/2.2_functions.R")

profvis(IS(
  h = Sn, ## h(x) function
  #h = default,
  h_mat_gen = h_mat_gen_2, ## Generate matrix of h(x) values
  num_steps = 100, 
  num_paths = 1000,
  theta = -0.2,
  a = -1.9,
  b = 2.0,
  sigma_switch = FALSE
))


## *** ----
## BENCHMARKING ----

## > TESTS ----
## >> Sn_mat_gen_1 vs Sn_mat_gen_2 ----

## >> Runif ----

## >> Runtime wrt. n ----



## *** ----
## PLOT

## Plot gn against theta
{
  num_curves = 1000
  theta_min <- -1
  theta_max <- -0.01
  tmp_vals <- matrix(numeric(1000), 100, num_curves)
  for(i in 1:num_curves) {
    tmp_vals[, i] <- sapply(seq(theta_min, theta_max, 0.01), function(theta) {gn_1(xg_mat[,i], theta, a, b)})
  }
  #plot(seq(theta_min, theta_max, 0.01), tmp_vals, pch=16, cex = 0.3, xlab = "theta", ylab = "gn", log="y", main="lin-log")
  #plot(seq(theta_min, theta_max, 0.01), tmp_vals[, 1], pch=16, cex = 0.3, xlab = "theta", ylab = "gn", main="lin-lin")
  plot_data <- melt(tmp_vals) ## Create data frame with columns as groups
  plot_data$Var1 <- rep(seq(theta_min, theta_max, 0.01), num_curves)
  theta_test_plot = ggplot() +
    geom_line(data = plot_data, aes(x = Var1, y = value, group = Var2, colour=factor(Var2))) +
    #scale_y_continuous(trans='log10') +
    labs(title = "lin-log", subtitle = substitute(paste("gn wrt. theta for ",  num_curves, " simulated paths"), list(num_curves = num_curves)), x = "theta", y = "gn") +
    theme(legend.position = "none")
  theta_test_plot
}




