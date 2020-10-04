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
## DENSITIES ====

## Test that g is a density function
x <- x_mat_1(100, 1, rfunc = runif(num_steps*num_paths, -1.9, 2.0))
g_i(x, theta = 1, -1.9, 2)
phi(1, -1.9, 2)
tmp_int <- integrate(function(x){g_i(x, theta = 1, -1.9, 2)}, -1.9, 2.0)
str(tmp_int)
## g is a density (integrates to 1, except for theta = 0, where g is undefined)

df = data.frame(x_mat_1(1e5, 1, rfunc = runif(num_steps*num_paths, -1.9, 2.0)))
names(df)="X1"
ggplot(df, aes(x=X1)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 10*log(100, 10))+
  geom_density(alpha=.2, fill="#FF6666") 


## gn
{
  num_steps = 100
  num_paths = 1e5
  theta = 0.5
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  xg_mat <- xg_gen(p_vals, theta = theta, a = -1.9, b = 2)
  
  ## Histogram of num_steps*num_paths values drawn from g
  df = data.frame(matrix(xg_mat[, 1:1000], num_steps*1000, 1))
  names(df)="X1"
  ggplot(df, aes(x=X1)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 10*log(1e5, 10))+
    geom_density(alpha=.2, fill="#FF6666") +
    labs(title = paste("theta = ", theta))
}

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


## ***
## OPTIMAL THETA ====
{
  set.seed(42)
  num_steps = 100
  num_paths = 10000
  cut = 1e4 ## omit the left side of graph
  theta_test_vector <- seq(-0.9, 0.9, 0.1) #c(-0.6, -0.2, 0.2, 0.6)
  num_test_runs <- length(theta_test_vector)
  is_res <- matrix(numeric(num_test_runs*(num_paths-cut)), num_paths-cut, num_test_runs)
  
  
  for (i in 1:num_test_runs) {
    is_res_tmp <- IS(
      h = default,
      h_vect_gen = h_vect_gen_2, ## Generate vector of h(x) values
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
  theta_test_plot_1 = ggplot() +
    geom_line(data = plot_data, aes(x = Var1, y = value, group = Var2, colour=factor(Var2, labels = theta_test_vector))) +
    #scale_y_continuous(trans='log10') +
    labs(title = "theta", subtitle = "", x = "num_steps", y = "mu_hat", color = "theta") 
  #+
  #theme(legend.position = "none")
  
  #png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_test_plot_1.png', width=1800, height=1200, res=300)
  #theta_test_plot
  #dev.off()
}

## ***
## xg ====
## Ginv wrt theta
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
## Ginv wrt p
{
  num_steps = 100
  num_paths = 1e3
  a = -1.9
  b = 2.0
  theta_vals = c(0, -0.9)
  
  p_vals <- seq(0.0, 1.0, 0.01)
  #test_seq <- test_seq[which(test_seq == 0)] ## remove 0 from list
  # p_vals <- matrix(
  #   runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  # )
  
  p_test_a <- sapply(
    p_vals,
    function(p) {
      a = -1.9
      b = 2.0
      (p * b) + ((1 - p) * a)
    }
  )
  p_test_b <- sapply(
    p_vals,
    function(p) {
      theta=-1
      (
        log(
          p*(
            exp(b * theta)-exp(a * theta)
          ) + exp(a * theta)
        )
      ) / theta
    }
  )
  p_test_c <- sapply(
    p_vals,
    function(p) {
      theta=-0.6
      (
        log(
          p*(
            exp(b * theta)-exp(a * theta)
          ) + exp(a * theta)
        )
      ) / theta
    }
  )
  p_test_d<- sapply(
    p_vals,
    function(p) {
      theta=-0.2
      (
        log(
          p*(
            exp(b * theta)-exp(a * theta)
          ) + exp(a * theta)
        )
      ) / theta
    }
  )
  p_test_e<- sapply(
    p_vals,
    function(p) {
      theta=0.2
      (
        log(
          p*(
            exp(b * theta)-exp(a * theta)
          ) + exp(a * theta)
        )
      ) / theta
    }
  )
  p_test_f<- sapply(
    p_vals,
    function(p) {
      theta=0.6
      (
        log(
          p*(
            exp(b * theta)-exp(a * theta)
          ) + exp(a * theta)
        )
      ) / theta
    }
  )
  p_test_g<- sapply(
    p_vals,
    function(p) {
      theta=1
      (
        log(
          p*(
            exp(b * theta)-exp(a * theta)
          ) + exp(a * theta)
        )
      ) / theta
    }
  )
  
  
  theta_test_plot_3 <- qplot(p_vals, p_test_a, xlab="p", ylab="Ginv(x)", colour = I("red")) + 
    geom_point(aes(y = p_test_b), colour = I("green")) +
    geom_point(aes(y = p_test_c), colour = I("orange")) +
    geom_point(aes(y = p_test_d), colour = I("blue")) +
    geom_point(aes(y = p_test_e), colour = I("black")) +
    geom_point(aes(y = p_test_f), colour = I("cyan")) +
    geom_point(aes(y = p_test_g), colour = I("magenta")) +
    ylim(-2, 2.1) +
    labs(title = "x-values from g-dristribution wrt p", subtitle = "red: theta = 0, green: theta = -1, magenta: theta = 1")
}
png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/git/CompStat/Assignment 2.2_group2/images/theta_test_plot_3.png', width=1800, height=1200, res=300)
theta_test_plot_3
dev.off()

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
      h_vect_gen = h_vect_gen_2, ## Generate vector of h(x) values
      num_steps, 
      num_paths,
      theta = -0.2,
      a = -1.9,
      b = 2.0,
      sigma_switch = TRUE
    )
  )
}


## *** ----
## BENCHMARKING ----

## ***
## IS w/ h_vect_gen_1 vs IS w/ h_vect_gen_2
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

autoplot(bench_01) +
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  labs(title = "Benchmark of IS") +
  theme(legend.position = "none")



## ***
## h_vect_gen: if() vs apply() vs if() w/ preallocated vector
{
  num_steps = 100
  num_paths = 1e5
  theta = 0.5
  p_vals <- matrix(
    runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
  )
  xg_mat <- xg_gen(p_vals, theta = theta, a = -1.9, b = 2)
}

bench_02 <- microbenchmark(h_vect_gen_1(xg_mat, default), h_vect_gen_2(xg_mat, default), h_vect_gen_3(xg_mat, default))
levels(bench_02$expr) <- c("if()", "apply()", "if() prealloc")
summary(bench_02)

autoplot(bench_02) +
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  labs(title = "Benchmark of h_vect_gen") +
  theme(legend.position = "none")
## Conclusion: pre-alloc if() is faster
## ToDo: Plot runtimes of if() and apply() and if() w/ preallocated vector wrt. num_paths


## ***
## default(): any() vs if()
num_steps = 1e5
num_paths = 1 ## Only look at one path (vector)
p_vals <- matrix(
  runif(num_steps * num_paths, 0.0, 1.0), num_steps, num_paths, byrow = FALSE
)
xg_mat <- xg_gen(p_vals, theta = -0.2, a = -1.9, b = 2)

bench_03 <- microbenchmark(default(x = xg_mat), default_if(x = xg_mat), times = 100L)
levels(bench_03$expr) <- c("default: any())", "default: if()")
summary(bench_03)

autoplot(bench_03) +
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  labs(title = "Benchmark of default()", subtitle = paste(num_paths, " paths")) +
  theme(legend.position = "none")
## ToDo: Plot runtimes of default() and default_if() wrt. num_steps


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




