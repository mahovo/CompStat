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
  #library(plyr)
  library(dplyr)
  library(Rcpp)

  # Source functions
  source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/r/2.2_functions.R")
}

## Constants

num_paths = 1000 ## Number of paths
num_steps = 100 ## Number of steps in path


## *** ----
## DENSITIES ====

## Test that g is a density function
g(x, theta = 1, -1.9, 2)
phi(1, -1.9, 2)
tmp_int <- integrate(function(x){g(x, theta = 1, -1.9, 2)}, -1.9, 2.0)
str(tmp_int)
## g is a density (integrates to 1, except for theta = 0, where g is undefined)

## > Plot histograms ----
{
  num_steps = 100
  num_paths = 1e5
  x <- x_mat_1(100, 1, rfunc = runif(num_steps*num_paths, -1.9, 2.0))
  
  df = data.frame(x_mat_1(1e5, 1, rfunc = runif(num_steps*num_paths, -1.9, 2.0)))
  names(df)="X1"
  dens_plot_01 <- ggplot(df, aes(x=X1)) + 
    xlab("x") +
    geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 10*log(100, 10))+
    geom_density(alpha=.2, fill="#FF6666") +
    labs(title = "Density of g", subtitle = paste("theta = ", 0))
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/dens_plot_01.png', width=1800, height=1200, res=300)
# dens_plot_01
# dev.off()

{
  num_steps = 100
  num_paths = 1e5
  theta = 1.0
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  xg_mat <- xg_gen(p_vals, theta = theta, a = -1.9, b = 2)
  
  ## Histogram of num_steps*num_paths values drawn from g
  df = data.frame(matrix(xg_mat[, 1:1000], num_steps*1000, 1))
  names(df)="X1"
  dens_plot_02 <- ggplot(df, aes(x=X1)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 10*log(1e5, 10))+
    geom_density(alpha=.2, fill="#FF6666") +
    xlab("x") +
    labs(title = "Density of g", subtitle = paste("theta = ", theta))
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/dens_plot_02.png', width=1800, height=1200, res=300)
# dens_plot_02
# dev.off()

{
  num_steps = 100
  num_paths = 1e5
  theta = -1
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  xg_mat <- xg_gen(p_vals, theta = theta, a = -1.9, b = 2)
  
  ## Histogram of num_steps*num_paths values drawn from g
  df = data.frame(matrix(xg_mat[, 1:1000], num_steps*1000, 1))
  names(df)="X1"
  dens_plot_03 <- ggplot(df, aes(x=X1)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 10*log(1e5, 10))+
    geom_density(alpha=.2, fill="#FF6666") +
    xlab("x") +
    labs(title = "Density of g", subtitle = paste("theta = ", theta))
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/dens_plot_03.png', width=1800, height=1200, res=300)
# dens_plot_03
# dev.off()

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

## This test is not relevant. g integrates to 1 for all theta, except theta==
# {
#   
#   num_steps = 100
#   num_paths = 1e4
#   a=-1.9
#   b=2
# 
#   integrate(function(x){g(x, theta = 1, a, b)}, a, b)
#   
#   theta_test <- sapply(
#     theta_vals,
#     function(test) {integrate(function(x){g(x, theta = test, a, b)}, a, b)$value}
#   )
#   plot_optimal_theta_sim = qplot(theta_vals, theta_test, xlab="theta", ylab="integral of g") + 
#     geom_vline(xintercept = 1, color = "red", size = 0.5) +
#     geom_hline(yintercept = 1, color = "red", size = 0.5) +
#     labs(title = "Optimal theta", subtitle = "For density integrating to 1")
#   plot_optimal_theta_sim
# }

## IS

## Simulating x values with g:

## Test that x values are in [-1.9, 2.0]
{
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  range(xg_gen(p_vals, theta=1, a=-1.9, b=2.0))
}




## Compare IS and MC for Sn-function
num_paths = 1000 ## Number of paths
num_steps = 100 ## Number of steps in path
{
set.seed(123)
is_results <- IS(
  h = Sn, ## h(x) function
  #h = default,
  h_vect_gen = h_vect_gen_3, ## Generate matrix of h(x) values
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
  h_vect_gen_3, ## Generate vector of h(x) values
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
    h_vect_gen = h_vect_gen_3, ## Generate vector of h(x) values
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
  IS_plot_01 <- qplot(1:num_paths, is_results$mu_hat) + 
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
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/IS_plot_01.png', width=1800, height=1200, res=300)
# IS_plot_01
# dev.off()



{
  ## mu_hat = expected probability of default
  set.seed(42)
  mc_results_2 <- MCI(
    default, ## h(x) function
    h_vect_gen_3, ## Generate vector of h(x) values
    num_steps,
    num_paths,
    runif(num_steps*num_paths, -1.9, 2.0),
    sigma_switch = TRUE
  )
  mc_results_2$mu_hat
}
{
  IS_plot_02 <- qplot(1:num_paths, mc_results_2$mu_hat) + 
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
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/IS_plot_02.png', width=1800, height=1200, res=300)
# IS_plot_02
# dev.off()


## ***
## OPTIMAL THETA ====
## > mu_hat wrt num_steps ----
{
  set.seed(42)
  num_steps = 100
  num_paths = 1e5
  theta_test_vector <- seq(-0.9, 0.9, 0.1) #c(-0.6, -0.2, 0.2, 0.6)
  num_test_runs <- length(theta_test_vector)
  is_res <- matrix(numeric(num_test_runs*(num_paths)), num_paths, num_test_runs)
  
  
  for (i in 1:num_test_runs) {
    is_res_tmp <- IS(
      h = default,
      h_vect_gen = h_vect_gen_3, ## Generate vector of h(x) values
      num_steps, 
      num_paths,
      theta_test_vector[i],
      a = -1.9,
      b = 2.0,
      sigma_switch = TRUE
    )
    is_res[, i] <- is_res_tmp$mu_hat
    # print("")
    # print(range(is_res_tmp$w_star))
    # print(range(is_res_tmp$h_vect))
    # print(range(is_res_tmp$w_star * is_res_tmp$h_vect))
  }
  
  plot_data <- melt(is_res) ## Create data frame with columns as groups
  plot_data$Var1 <- rep(1:num_paths, num_test_runs)
  theta_test_plot_1 = ggplot() +
    geom_line(data = plot_data, aes(x = Var1, y = value, group = Var2, colour=factor(Var2, labels = theta_test_vector))) +
    #scale_y_continuous(trans='log10') +
    labs(title = "theta", subtitle = "", x = "num_steps", y = "mu_hat", color = "theta") 
  #+
  #theme(legend.position = "none")
}
#png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_test_plot_1.png', width=1800, height=1200, res=300)
#theta_test_plot_1
#dev.off()

## > mu_hat wrt num_steps with cut----
{
  set.seed(42)
  num_steps = 100
  num_paths = 1e5
  cut = 1e4 ## omit the left side of graph
  theta_test_vector <- seq(-0.9, 0.9, 0.1) #c(-0.6, -0.2, 0.2, 0.6)
  num_test_runs <- length(theta_test_vector)
  is_res <- matrix(numeric(num_test_runs*(num_paths-cut)), num_paths-cut, num_test_runs)
  
  for (i in 1:num_test_runs) {
    is_res_tmp <- IS(
      h = default,
      h_vect_gen = h_vect_gen_3, ## Generate vector of h(x) values
      num_steps, 
      num_paths,
      theta_test_vector[i],
      a = -1.9,
      b = 2.0,
      sigma_switch = TRUE
    )
    is_res[, i] <- is_res_tmp$mu_hat[-(1:cut)]
  }
  
  plot_data <- melt(is_res) ## Create data frame with columns as groups
  plot_data$Var1 <- rep(1:num_paths, num_test_runs)
  theta_test_plot_1_cut = ggplot() +
    geom_line(data = plot_data, aes(x = Var1, y = value, group = Var2, colour=factor(Var2, labels = theta_test_vector))) +
    #scale_y_continuous(trans='log10') +
    labs(title = "theta", subtitle = "", x = "num_steps", y = "mu_hat", color = "theta") 
  #+
  #theme(legend.position = "none")
}
#png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_test_plot_1.png', width=1800, height=1200, res=300)
#theta_test_plot_1_cut
#dev.off()

## ***
## > Ginv wrt theta
{
  a = -1.9
  b = 2.0
  p_range = c(0, 1)
  
  theta_vals <- seq(-0.9, 0.9, 0.02)
  #test_seq <- test_seq[which(test_seq == 0)] ## remove 0 from list
  # p_vals <- matrix(
  #   runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  # )
  xg_mat <- xg_gen(p_vals, theta = theta_seq[1], a = -1.9, b = 2)
  
  theta_test_a <- sapply(
    theta_vals,
    function(theta) {
      if(theta != 0.0) {
        (
          log(
            p_range[1]*(
              exp(b * theta)-exp(a * theta)
            ) + exp(a * theta)
          )
        ) / theta
      } else {
        a = -1.9
        b = 2.0
        (p_range[1] * b) + ((1 - p_range[1]) * a)}
    }
  )
  theta_test_b <- sapply(
    theta_vals,
    function(theta) {
      if(theta != 0.0) {
        (
          log(
            p_range[2]*(
              exp(b * theta)-exp(a * theta)
            ) + exp(a * theta)
          )
        ) / theta
      } else {
        a = -1.9
        b = 2.0
        (p_range[2] * b) + ((1 - p_range[2]) * a)}
    }
  )
  
  theta_test_plot_2 <- qplot(theta_vals, theta_test_a, xlab="theta", ylab="Ginv(x)", colour = I("red")) + 
    geom_point(aes(y = theta_test_b), colour = I("blue")) +
    ylim(-2, 10) +
    labs(title = "Optimal theta", subtitle = "x-values from g-dristribution")
}
#png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_test_plot_2.png', width=1800, height=1200, res=300)
#theta_test_plot_2
#dev.off()


## ***
## > Ginv wrt p ----
{
  num_steps = 100
  num_paths = 1e3
  a = -1.9
  b = 2.0
  
  theta_vals <- c(-1.0, -0.6, -0.2, 0.2, 0.6, 1.0)
  p_vals <- seq(0.0, 1.0, 0.01)
  #test_seq <- test_seq[which(test_seq == 0)] ## remove 0 from list
  # p_vals <- matrix(
  #   runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  # )
  
  p_test_zero <- sapply(
    p_vals,
    function(p) {
      a = -1.9
      b = 2.0
      (p * b) + ((1 - p) * a)
    }
  )
  p_test_mat <- matrix(numeric(length(theta_vals ) * length(p_vals)), length(p_vals), length(theta_vals ))
  for (i in 1:length(theta_vals)) {
    p_test_mat[, i] <- sapply(
      p_vals,
      function(p) {
        theta=theta_vals[i]
        (
          log(
            p*(
              exp(b * theta)-exp(a * theta)
            ) + exp(a * theta)
          )
        ) / theta
      }
    )
  }
  
  theta_test_plot_3 <- qplot(p_vals, p_test_zero, xlab="p", ylab="Ginv(x)", colour = I("red")) + 
    geom_point(aes(y = p_test_mat[, 1]), colour = I("green")) +
    geom_point(aes(y = p_test_mat[, 2]), colour = I("orange")) +
    geom_point(aes(y = p_test_mat[, 3]), colour = I("blue")) +
    geom_point(aes(y = p_test_mat[, 4]), colour = I("black")) +
    geom_point(aes(y = p_test_mat[, 5]), colour = I("cyan")) +
    geom_point(aes(y = p_test_mat[, 6]), colour = I("magenta")) +
    ylim(-2, 2.1) +
    labs(title = "x-values from g-distribution wrt p", subtitle = "red: theta = 0, green: theta = -1, magenta: theta = 1")
  theta_test_plot_3
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_test_plot_3.png', width=1800, height=1200, res=300)
# theta_test_plot_3
# dev.off()

## > sd(h(x)) wrt theta ----
{
  num_steps = 100
  num_paths = 1e6
  theta_vals = seq(-0.9, 0.9, 0.01)
  #theta_vals = c(-0.5, 0.5)
  num_test_runs = length(theta_vals)
  
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  # xg_mat <- xg_gen(p_vals, theta = theta, a = -1.9, b = 2)
  # h_vect_gen_3(xg_mat, default)
  
  sd_vect <- numeric(num_test_runs)
  for (i in 1:num_test_runs) {
    xg_mat <- xg_gen(p_vals, theta = theta_vals[i], a = -1.9, b = 2)
    sd_vect_tmp <- sd(h_vect_gen_3(xg_mat, default))
    sd_vect[i] <- sd_vect_tmp
  }
  range(sd_vect)
  
  ## Which theta value gives highest variance?
  theta_opt = theta_vals[which(sd_vect==max(sd_vect))]
  ## theta = -0.27
  sd_test_plot_1 <- qplot(theta_vals, sd_vect, xlab="theta", ylab="sd(h(x))") +
    geom_vline(xintercept = theta_opt, color = "red", size = 0.5) +
    labs(title = "sd(h(x)) for x-values from g-distribution wrt theta", subtitle = paste("Optimal theta: ", theta_opt, ", Number of paths: ", num_paths))
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/sd_test_plot_1.png', width=1800, height=1200, res=300)
# sd_test_plot_1
# dev.off()


## ***
## > sd(h(x)*w_star) wrt theta ----
{
  num_steps = 100
  num_paths = 1e2
  theta_vals = seq(-1.0, 1.0, 0.1)
  num_test_runs = length(theta_vals)
  
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  # xg_mat <- xg_gen(p_vals, theta = theta, a = -1.9, b = 2)
  # h_vect_gen_3(xg_mat, default)
  
  
  ## For each theta value:
  ## - Generate num_paths paths of num_steps steps
  ## - Calcuate a vector gn-density density values for that theta value
  sd_vect <- numeric(num_test_runs)
  g_dens_mat <- matrix(numeric(num_test_runs * num_paths), num_paths, num_test_runs)
  w_star_mat <- matrix(numeric(num_test_runs * num_paths), num_paths, num_test_runs)
  for (i in 1:num_test_runs) {
    xg_mat <- xg_gen(p_vals, theta = theta_vals[i], a = -1.9, b = 2)
    g_dens <- apply(xg_mat, 2, function(x) {gn_1(x, theta = theta_vals[i], a = -1.9, b = 2)})
    f_dens <- replicate(num_paths, (1/3.9)^num_steps)
    #f_dens <- replicate(num_paths, (1/3.9))
    g_dens_mat[, i] <- g_dens
    w_star_mat[, i] <-  f_dens / g_dens
    sd_vect_tmp <- sd(h_vect_gen_3(xg_mat, default) * w_star_mat[, i])
    sd_vect[i] <- sd_vect_tmp
  }
  # range(sd_vect)
  # range(w_star_mat)
  # range(g_dens_mat)
  # sum(g_dens_mat == Inf) ## 1000 g_dens vals are 0. I.e. those for theta=0
  # ((1/3.9)^num_steps) / min(g_dens_mat) ## max w_star

  
  g_dens_mat_plot_data <- melt(g_dens_mat) ## Create data frame with columns as groups
  g_dens_mat_plot_data$Var2 <- rep(theta_vals, each = num_paths)
  
  theta_g_dens_plot <-  ggplot() +
    geom_point(data = g_dens_mat_plot_data, aes(x = Var2, y = value, group = Var1, colour=factor(Var1)), size = 0.2) +
    geom_line(data = g_dens_mat_plot_data, aes(x = Var2, y = value, group = Var1, colour=factor(Var1)), size = 0.2) +
    scale_y_continuous(trans='log10') +
    geom_hline(yintercept=7.82599e-60) +
    geom_text(data=data.frame(x=0, y=7.82599e-60), aes(x, y), label="fn", vjust=-0.5, hjust=17.5) +
    labs(title = "lin-log", subtitle = paste("gn wrt. theta for",  num_paths, "simulations at each theta"), x = "theta", y = "gn") +
    theme(legend.position = "none")

  
  w_star_mat_plot_data <- melt(w_star_mat) ## Create data frame with columns as groups
  w_star_mat_plot_data$Var2 <- rep(theta_vals, each = num_paths)
  
  theta_w_star_plot <- ggplot() +
    geom_point(data = w_star_mat_plot_data, aes(x = Var2, y = value, group = Var1, colour=factor(Var1)), size = 0.2) +
    geom_line(data = w_star_mat_plot_data, aes(x = Var2, y = value, group = Var1, colour=factor(Var1)), size = 0.2) +
    scale_y_continuous(trans='log10') +
    labs(title = "lin-log", subtitle = paste("w_star wrt. theta for",  num_paths, "simulations at each theta"), x = "theta", y = "w_star") +
    theme(legend.position = "none")

}
{
theta_g_dens_plot
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_g_dens_plot.png', width=1800, height=1200, res=300)
# theta_g_dens_plot
# dev.off()

theta_w_star_plot
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_w_star_plot.png', width=1800, height=1200, res=300)
# theta_w_star_plot
# dev.off()
}
{ 
  ## Which theta value gives highest variance?
  theta_opt = theta_vals[which(sd_vect==max(sd_vect))]

  sd_test_plot_2_linlog <- qplot(theta_vals, sd_vect, xlab="theta", ylab="sd(h(x) * w_star)") +
    geom_vline(xintercept = theta_opt, color = "red", size = 0.5) +
    scale_y_continuous(trans='log10') +
    labs(title = "sd(h(x) * w_star) for x-values from g-distribution wrt theta (lin-log)", subtitle = paste("Optimal theta: ", theta_opt, ", Number of paths: ", num_paths))

  sd_test_plot_2_linlin <- qplot(theta_vals, sd_vect, xlab="theta", ylab="sd(h(x) * w_star)") +
    geom_vline(xintercept = theta_opt, color = "red", size = 0.5) +
    #scale_y_continuous(trans='log10') +
    labs(title = "sd(h(x) * w_star) for x-values from g-distribution wrt theta (lin-lin)", subtitle = paste("Optimal theta: ", theta_opt, ", Number of paths: ", num_paths))

  sd_test_plot_2_linlin_cut <- qplot(theta_vals[-(1:95)], sd_vect[-(1:95)], xlab="theta", ylab="sd(h(x) * w_star)") +
    geom_vline(xintercept = theta_opt, color = "red", size = 0.5) +
    xlim(0, 1) +
    #scale_y_continuous(trans='log10') +
    labs(title = "sd(h(x) * w_star) for x-values from g-distribution wrt theta (lin-lin)", subtitle = paste("Optimal theta: ", theta_opt, ", Number of paths: ", num_paths))
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/sd_test_plot_2_linlog.png', width=1800, height=1200, res=300)
# sd_test_plot_2_linlog
# dev.off()
# 
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/sd_test_plot_2_linlin.png', width=1800, height=1200, res=300)
# sd_test_plot_2_linlin
# dev.off()
# 
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/sd_test_plot_2_linlin_cut.png', width=1800, height=1200, res=300)
# sd_test_plot_2_linlin_cut
# dev.off()

## ***
## > sd(w_star) wrt theta ----
{
  num_steps = 100
  num_paths = 1e4
  theta_vals = seq(-0.9, 0.9, 0.01)
  num_test_runs = length(theta_vals)

  # xg_mat <- xg_gen(p_vals, theta = theta, a = -1.9, b = 2)
  # h_vect_gen_3(xg_mat, default)
  
  sd_vect <- numeric(num_test_runs)
  for (i in 1:num_test_runs) {
    p_vals <- matrix(
      runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
    )
    xg_mat <- xg_gen(p_vals, theta = theta_vals[i], a = -1.9, b = 2)
    g_dens <- apply(xg_mat, 2, function(x) {gn_1(x, theta = theta_vals[i], a = -1.9, b = 2)})
    f_dens <- replicate(num_paths, (1/3.9)^num_steps)
    #f_dens <- replicate(num_paths, (1/3.9))
    w_star <-  f_dens / g_dens
    sd_vect_tmp <- sd(w_star)
    sd_vect[i] <- sd_vect_tmp
    print(" - - - - - ")
    print(paste("w_star: ", range(w_star)))
    print(paste("g_dens: ", range(g_dens)))
  }
  #range(sd_vect)
  
  
  ## Which theta value gives highest variance?
  theta_opt = theta_vals[which(sd_vect==max(sd_vect))]
  
  sd_test_plot_3_linlog <- qplot(theta_vals, sd_vect, xlab="theta", ylab="sd(w_star)") +
    geom_vline(xintercept = theta_opt, color = "red", size = 0.5) +
    scale_y_continuous(trans='log10') +
    labs(title = "sd(w_star) for x-values from g-distribution wrt theta (lin-log)", subtitle = paste("Optimal theta: ", theta_opt, ", Number of paths: ", num_paths))
  
  sd_test_plot_3_linlin <- qplot(theta_vals, sd_vect, xlab="theta", ylab="sd(w_star)") +
    geom_vline(xintercept = theta_opt, color = "red", size = 0.5) +
    #scale_y_continuous(trans='log10') +
    labs(title = "sd(w_star) for x-values from g-distribution wrt theta (lin-lin)", subtitle = paste("Optimal theta: ", theta_opt, ", Number of paths: ", num_paths))
  
  sd_test_plot_3_linlin_cut <- qplot(theta_vals[-(1:95)], sd_vect[-(1:95)], xlab="theta", ylab="sd(w_star)") +
    geom_vline(xintercept = theta_opt, color = "red", size = 0.5) +
    xlim(0, 1) +
    #scale_y_continuous(trans='log10') +
    labs(title = "sd(w_star) for x-values from g-distribution wrt theta (lin-lin)", subtitle = paste("Optimal theta: ", theta_opt, ", Number of paths: ", num_paths))
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/sd_test_plot_3_linlog.png', width=1800, height=1200, res=300)
# sd_test_plot_3_linlog
# dev.off()

# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/sd_test_plot_3_linlin.png', width=1800, height=1200, res=300)
# sd_test_plot_3_linlin
# dev.off()

# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/sd_test_plot_3_linlin_cut.png', width=1800, height=1200, res=300)
# sd_test_plot_3_linlin_cut
# dev.off()


## ***
## > sd mu_hat wrt theta ----
{
  num_steps = 100
  num_paths = 1e5
  theta_vals = seq(-1.0, 1.0, 0.1)
  #theta_vals = c(-0.5, 0.5)
  num_test_runs = length(theta_vals)
  
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  
  sd_vect <- numeric(num_test_runs)
  for (i in 1:num_test_runs) {
    xg_mat <- xg_gen(p_vals, theta = theta_vals[i], a = -1.9, b = 2)
    is_result <- IS(
      default,
      h_vect_gen = h_vect_gen_4, ## Generate vector of h(x) values
      num_steps, 
      num_paths,
      theta = theta_vals[i],
      a = -1.9,
      b = 2.0,
      sigma_switch = TRUE
    )
    sd_vect[i] <- is_result$sigma[num_paths]
  }
  print(paste("Range of sd(mu_hat): [",range(sd_vect)[1], ", ", range(sd_vect)[2], "]"))
  
  
  ## Which theta value gives lowest variance?
  theta_opt = theta_vals[which(sd_vect==min(sd_vect))]
  ## theta = -0.27
  sd_test_plot_4 <- qplot(theta_vals, sd_vect, xlab="theta", ylab="sd(mu_hat)") +
    geom_vline(xintercept = theta_opt, color = "red", size = 0.5) +
    scale_y_continuous(trans='log10') +
    labs(title = "sd(mu_hat) for x-values from g-distribution wrt theta", subtitle = paste("Optimal theta: ", theta_opt, ", Number of paths: ", num_paths))
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/sd_test_plot_4.png', width=1800, height=1200, res=300)
# sd_test_plot_4
# dev.off()



## > gn wrt. theta ----
## lin-log
{
  num_curves = 1e1
  theta_min <- -1
  theta_max <- 1
  theta_vals <- seq(theta_min, theta_max, 0.02)
  #theta_vals <- theta_vals[-which(theta_vals == 0)] # Remove theta = 0, see phi()
  
  num_steps <- 100
  num_paths <- 1
  a = -1.9
  b = 2.0
  
  f_dens <- (1/3.9)^num_steps
  
  tmp_vals <- matrix(length(theta_vals) * num_curves, length(theta_vals), num_curves)
  ## Do 100 times:
  ## For each theta value, create a sample single paths of x-values from the g-distribution, 
  ## and calculate a gn value
  
  for(i in 1:num_curves) {
    tmp_vals[, i] <- sapply(theta_vals, function(theta) {
      p_vals <- matrix(
        runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
      )
      xg_mat <- xg_gen(p_vals, theta, a = -1.9, b = 2)
      gn_1(xg_mat, theta, a, b) ## A number
      #apply(xg_mat, 2, function(x) {gn_1(x, theta, a = -1.9, b = 2)})
    })
  }
  plot_data <- melt(tmp_vals) ## Create data frame with columns as groups
  plot_data$Var1 <- rep(theta_vals, num_curves)
  theta_test_plot_linlog = ggplot() +
    geom_point(data = plot_data, aes(x = Var1, y = value, group = Var2, colour=factor(Var2)), size = 0.2) +
    geom_line(data = plot_data, aes(x = Var1, y = value, group = Var2, colour=factor(Var2)), size = 0.2) +
    scale_y_continuous(trans='log10') +
    geom_hline(yintercept=f_dens) +
    geom_text(data=data.frame(x=0, y=f_dens), aes(x, y), label="fn", vjust=-0.5, hjust=17.5) +
    labs(title = "lin-log", subtitle = substitute(paste("gn wrt. theta for ",  num_curves, " simulations"), list(num_curves = num_curves)), x = "theta", y = "gn") +
    theme(legend.position = "none")
  theta_test_plot_linlog
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_test_plot_linlog.png', width=1800, height=1200, res=300)
# theta_test_plot_linlog
# dev.off()

## Why spike at theta=0?
## theta =0 -> phi = 0 -> gn = inf -> w_star = 0
## For theta = 0, w_star should be 1, so mu_hat_IS = mu_hat_MC.
## Instead for theta = 0, all weights are 0.
## However, close to theta=0, mu_hat_IS  will be close to mu_hat_MC
range(tmp_vals)

## lin-log cut
{
  num_curves = 1e4
  theta_min <- -1
  theta_max <- -0.01
  theta_vals <- seq(theta_min, theta_max, 0.01)
  #theta_vals <- theta_vals[-which(theta_vals == 0)] # Remove theta = 0, see phi()
  
  tmp_vals <- matrix(length(theta_vals) * num_curves, length(theta_vals), num_curves)
  for(i in 1:num_curves) {
    tmp_vals[, i] <- sapply(theta_vals, function(theta) {gn_1(xg_mat[,i], theta, a, b)})
  }
  plot_data <- melt(tmp_vals) ## Create data frame with columns as groups
  plot_data$Var1 <- rep(theta_vals, num_curves)
  theta_test_plot_linlog_cut = ggplot() +
    geom_line(data = plot_data, aes(x = Var1, y = value, group = Var2, colour=factor(Var2))) +
    scale_y_continuous(trans='log10') +
    labs(title = "lin-log", subtitle = substitute(paste("gn wrt. theta for ",  num_curves, " simulated paths"), list(num_curves = num_curves)), x = "theta", y = "gn") +
    theme(legend.position = "none")
  #theta_test_plot_linlog_cut
}
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_test_plot_linlog_cut.png', width=1800, height=1200, res=300)
# theta_test_plot_linlog_cut
# dev.off()




## **************************************

## OOP ====  
## > CONSTRUCTORS ----



## > VALIDATORS ----



## > HELPERS ----




## *** ----
## PROFILING ----

source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/r/2.2_functions.R")

{
  set.seed(42)
  profvis(IS(
      h = default,
      h_vect_gen = h_vect_gen_3, ## Generate vector of h(x) values
      num_steps, 
      num_paths,
      theta = -0.27,
      a = -1.9,
      b = 2.0,
      sigma_switch = TRUE
    )
  )
}


## *** ----
## BENCHMARKING ----

## ***
## > IS w/ h_vect_gen_1 vs IS w/ h_vect_gen_2 ----
{
  bench_01 <- microbenchmark(
    IS(
      h = default,
      h_vect_gen = h_vect_gen_1, ## Generate vector of h(x) values
      num_steps, 
      num_paths,
      theta = -0.2,
      a = -1.9,
      b = 2.0,
      sigma_switch = TRUE
    ),
    IS(
      h = default,
      h_vect_gen = h_vect_gen_2, ## Generate vector of h(x) values
      num_steps, 
      num_paths,
      theta = -0.2,
      a = -1.9,
      b = 2.0,
      sigma_switch = TRUE
    ),
    times = 20L
  )
}
levels(bench_01$expr) <- c("loop", "apply")

bench_plot_01 <- autoplot(bench_01) +
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  labs(title = "Benchmark of IS") +
  theme(legend.position = "none")
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/bench_plot_01.png', width=1800, height=1200, res=300)
# bench_plot_01
# dev.off()


## ***
## > h_vect_gen: if() vs apply() vs if() w/ preallocated vector ----
{
  num_steps = 100
  num_paths = 1e5
  theta = -0.27
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  xg_mat <- xg_gen(p_vals, theta = theta, a = -1.9, b = 2)
}

bench_02 <- microbenchmark(h_vect_gen_1(xg_mat, default), h_vect_gen_2(xg_mat, default), h_vect_gen_3(xg_mat, default), h_vect_gen_4(xg_mat))
levels(bench_02b$expr) <- c("apply()", "if()", "if() prealloc", "rcpp")
bench_02
summary(bench_02)

bench_plot_02 <- autoplot(bench_02b) +
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  labs(title = "Benchmark of h_vect_gen") +
  theme(legend.position = "none")
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/bench_plot_02.png', width=1800, height=1200, res=300)
# bench_plot_02
# dev.off()

## Conclusion: pre-alloc if() is faster
## ToDo: Plot runtimes of if() and apply() and if() w/ preallocated vector wrt. num_paths


## ***
## > default(): any() vs if() ----
num_steps = 1e5
num_paths = 1 ## Only look at one path (vector)
p_vals <- matrix(
  runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
)
xg_mat <- xg_gen(p_vals, theta = -0.2, a = -1.9, b = 2)

bench_03 <- microbenchmark(default(x = xg_mat), default_if(x = xg_mat), times = 100L)
levels(bench_03$expr) <- c("default: any())", "default: if()")
summary(bench_03)

bench_plot_03 <- autoplot(bench_03) +
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  labs(title = "Benchmark of default()", subtitle = paste(num_paths, " paths")) +
  theme(legend.position = "none")
# png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/bench_plot_03.png', width=1800, height=1200, res=300)
# bench_plot_03
# dev.off()
## ToDo: Plot runtimes of default() and default_if() wrt. num_steps

## ***
## > Sn_mat_gen_1 vs Sn_mat_gen_2 ----

## ***
## > Runif ----

## ***
## > Runtime wrt. n ----

{
  minn = 1
  maxn = 5
  num_paths = 10^maxn
  num_steps = 100
  
  conf <- expand.grid(
    fun = c("h_vect_gen_1", "h_vect_gen_2", "h_vect_gen_3", "h_vect_gen_4"),
    n = 10^(minn:maxn)
  )
  
  set.seed(42)
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  xg_mat <- xg_gen(p_vals, theta = -0.27, a = -1.9, b = 2)
  
  calls <- paste0(conf[, 1], "(xg_mat[, 1:", conf[, 2], "])")
  expr_list <- lapply(calls, function(x) parse(text = x)[[1]])
  kern_benchmarks <- microbenchmark(list = expr_list, times = 40L)
  
  class(kern_benchmarks) <- "data.frame"
  kern_benchmarks <- dplyr::bind_cols(conf, expr = calls) %>% 
    dplyr::left_join(kern_benchmarks, .)
  
  kern_benchmarks$time = kern_benchmarks$time/1000 ## Convert from nanoseconds to milliseconds
  
  
  runtime_h_vect_gen_03 <- ggplot(kern_benchmarks, aes(x = n, y = time, color = fun)) + 
    geom_abline(intercept = 0, slope = 1, color = "gray", linetype = 2) +
    stat_summary(fun = "median", geom = "line") + 
    stat_summary(fun = "median", geom = "point") + 
    #facet_wrap(~ fun) + 
    scale_x_continuous(trans = "log10") + 
    scale_y_continuous("Time (ms)", trans = "log10") +
    #scale_y_continuous("Time (ms)", trans = "log10", 
    #                    breaks = c(1e5, 1e6, 1e7, 1e8), 
    #                    labels = c("0.1", "1", "10", "100")) +
    scale_color_discrete("Function:", labels = c("apply", "f (a)", "f (b)", "rcpp")) + 
    theme(legend.position="right")
  
  runtime_h_vect_gen_03
  # png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/runtime_h_vect_gen_03.png', width=1800, height=1200, res=300)
  # runtime_h_vect_gen_03
  # dev.off()
}