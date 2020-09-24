###########################################
## Computational Statistics
## KU Science 2020/Block 1
## Martin Hoshi Vognsen
###########################################


## Assignment 1_1: Density estimation ----

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
  source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/Assignment 1.1_group2/r/1.1_v02_functions.R")
  
  infrared <- read.table("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Exercises/data/infrared.txt", header = TRUE)
  F12 <- data.frame(infrared$F12)
  logF12 <- data.frame(log(infrared$F12))
}


## -> Empirical data ----
## logF12

## Martin's h
x = logF12[, 1]

{
  n = length(x)
  r <- r_hat(x, n)
  
  wiggle <- wiggle_gauss_vec(x, r)
  
  hn_hat_opt <- hn_hat(wiggle, n)
  hn_hat_opt
  
  h = hn_hat_opt
  f_hat <- kern_dens(x, h*sqrt(5)) ## Multiply by sqrt(5) to compare with density()
  f_hat_dens <- density(
    x, h, kernel = "epanechnikov", cut = 3*sqrt(5) ## We multiplied h by sqrt(5) in kern_dens()
  )
  hist(x, prob=TRUE, xlim = c(min(x), max(x)))
  lines(f_hat, type = "l", lwd = 4)
  lines(f_hat_dens, col = "red", lwd = 2)
}

plot(f_hat$x, f_hat$y - f_hat_dens$y, 
     type = "l", 
     lwd = 2,
     xlab = "x", 
     ylab = "Difference"
)

## KEDD's h
{
  n = length(x)
  r <- r_hat(x, n)
  
  wiggle <- wiggle_gauss_vec(x, r)
  
  hn_hat_opt <- hn_hat(wiggle, n)
  hn_hat_opt
  
  h = kedd::h.amise(x, kernel = "epanechnikov")$h
  f_hat <- kern_dens(x, h*sqrt(5)) ## Multiply by sqrt(5) to compare with density()
  f_hat_dens <- density(x, h, kernel = "epanechnikov", cut = 3*sqrt(5)) ## We multiplied h by sqrt(5) in kern_dens()
  hist(x, prob=TRUE, xlim = c(min(x), max(x)))
  lines(f_hat, type = "l", lwd = 4)
  lines(f_hat_dens, col = "red", lwd = 2)
}



## -> Simulated data ----
## Standard normal, teoretisk, n=10000
{
  n = 10000
  x = rnorm(n)
  x_df = data.frame(x)
  q = seq(-5, 5, length.out = n)
  norm_dens = sapply(
    q, function(q) {1/(sqrt(2*pi)) * exp(-0.5*q^2)}
  )
  q_df = data.frame(q = q, norm_dens = norm_dens)
}
plot(x, pch = 16, cex = 0.2)

{
  n = length(x)
  r <- r_hat(x, n)
  
  wiggle <- wiggle_gauss_vec(x, r)
  
  hn_hat_opt <- hn_hat(wiggle, n)
  cat(paste("h_opt = ", hn_hat_opt))
  
  h = hn_hat_opt
  f_hat <- kern_dens(x, h*sqrt(5)) ## Multiply by sqrt(5) to compare with density()
  f_hat_dens <- density(x, h, kernel = "epanechnikov", cut = 3*sqrt(5)) ## We multiplied h by sqrt(5) in kern_dens()
  
  png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/Assignment 1.1_group2/images/stdnorm.png', width=1800, height=1200, res=300)
    hist(x, prob=TRUE, main = paste("h_opt = ", round(hn_hat_opt, 4), "n=10000"))
    points(q, norm_dens, col='blue')
    lines(f_hat, col='yellow', type = "l", lwd = 4)
    lines(f_hat_dens, col = 'red', lwd = 2)
    legend(-5, 0.39, 
           legend = c("N(0,1)", "kern_dens", "density"), 
           col = c('blue', 
                   'yellow',
                   'red'), 
           #pch = c(17,19), 
           # bty = "n", 
           lty = c(1, 1, 1),
           lwd=2 ,
           pt.cex = 1, 
           cex = 1, 
           text.col = "black", 
           horiz = F , 
           inset = c(0.3, 0.4))
  dev.off()
}
{
png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/Assignment 1.1_group2/images/diff.png', width=1800, height=1200, res=300)
  plot(f_hat$x, f_hat$y - f_hat_dens$y, 
  type = "l", 
  lwd = 2,
  xlab = "x", 
  ylab = "Difference")
dev.off()
}

## Standard normal, teoretisk, n=1000
{
  n = 1000
  x = rnorm(n)
  x_df = data.frame(x)
  q = seq(-5, 5, length.out = n)
  norm_dens = sapply(
    q, function(q) {1/(sqrt(2*pi)) * exp(-0.5*q^2)}
  )
  q_df = data.frame(q = q, norm_dens = norm_dens)
}
plot(x, pch = 16, cex = 0.2)

{
  n = length(x)
  r <- r_hat(x, n)
  
  #wiggle <- wiggle_gauss_outer(x, r)
  wiggle <- wiggle_gauss_vec(x, r)
  
  hn_hat_opt <- hn_hat(wiggle, n)
  cat(paste("h_opt = ", hn_hat_opt))
  
  h = hn_hat_opt
  f_hat <- kern_dens(x, h*sqrt(5)) ## Multiply by sqrt(5) to compare with density()
  f_hat_dens <- density(x, h, kernel = "epanechnikov", cut = 3*sqrt(5)) ## We multiplied h by sqrt(5) in kern_dens()
  
  hist(x, prob=TRUE, xlim = c(min(x), max(x)), main = paste("h_opt = ", round(hn_hat_opt, 4)))
  points(q, norm_dens, col='blue')
  lines(f_hat, type = "l", lwd = 4)
  lines(f_hat_dens, col = "red", lwd = 2)
}


## Standard normal, mixed
{
  n = 10000
  x1 = rnorm(floor(0.5*n), -0.3, 0.1)
  x2 = rnorm(floor(0.5*n), 0.3, 0.1)
  x = c(x1, x2)
  x_df = data.frame(x)
  q = seq(-5, 5, length.out = n)
  norm_dens = sapply(
    q, function(q) {1/(sqrt(2*pi)) * exp(-0.5*q^2)}
  )
  q_df = data.frame(q = q, norm_dens = norm_dens)
}
plot(x, pch = 16, cex = 0.2)

{
  n = length(x)
  r <- r_hat(x, n)
  
  wiggle <- wiggle_gauss_vec(x, r)
  
  hn_hat_opt <- hn_hat(wiggle, n)
  cat(paste("h_opt = ", hn_hat_opt))
  
  #png('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/2020/Assignments/Assignment 1.1_group2/images/bimodal.png', width=1800, height=1200, res=300)
    f_hat <- kern_dens(x, h*sqrt(5)) ## Multiply by sqrt(5) to compare with density()
    f_hat_dens <- density(x, h, kernel = "epanechnikov", cut = 3*sqrt(5)) ## We multiplied h by sqrt(5) in kern_dens()
    hist(x, prob=TRUE, xlim = c(min(x), max(x)), main = paste("h_opt = ", round(hn_hat_opt, 4)))
    lines(f_hat, type = "l", lwd = 4)
    lines(f_hat_dens, col = "red", lwd = 2)
  #dev.off()
}


plot(f_hat$x, f_hat$y - f_hat_dens$y, 
     type = "l", 
     lwd = 2,
     xlab = "x", 
     ylab = "Difference"
)

## GED,
library(fGarch)
{
  n = 1000
  x = rged(floor(0.5*n), 0, 0.2, 0.5)
}
plot(x, pch = 16, cex = 0.2)

{
  n = length(x)
  r <- r_hat(x, n)
  
  wiggle <- wiggle_gauss_outer(x, r)
  
  hn_hat_opt <- hn_hat(wiggle, n)
  cat(paste("h_opt = ", hn_hat_opt))
  
  h = hn_hat_opt
  f_hat <- kern_dens(x, h*sqrt(5)) ## Multiply by sqrt(5) to compare with density()
  f_hat_dens <- density(x, h, kernel = "epanechnikov", cut = 3*sqrt(5)) ## We multiplied h by sqrt(5) in kern_dens()
  
  hist(x, prob=TRUE, xlim = c(min(x), max(x)), main = paste("h_opt = ", round(hn_hat_opt, 4)))
  lines(f_hat, type = "l", lwd = 4)
  lines(f_hat_dens, col = "red", lwd = 2)
}

## GED, bimodal
library(fGarch)
{
  n = 1000
  x1 = rged(floor(0.5*n), -0.8, 0.2, 5)
  x2 = rged(floor(0.5*n), 0.8, 0.2, 5)
  x = c(x1, x2)
}
plot(x, pch = 16, cex = 0.2)

{
  n = length(x)
  r <- r_hat(x, n)
  
  wiggle <- wiggle_gauss_outer(x, r)
  
  hn_hat_opt <- hn_hat(wiggle, n)
  cat(paste("h_opt = ", hn_hat_opt))
  
  h = hn_hat_opt
  f_hat <- kern_dens(x, h*sqrt(5)) ## Multiply by sqrt(5) to compare with density()
  f_hat_dens <- density(x, h, kernel = "epanechnikov", cut = 3*sqrt(5)) ## We multiplied h by sqrt(5) in kern_dens()
  
  hist(x, prob=TRUE, xlim = c(min(x), max(x)), main = paste("h_opt = ", round(hn_hat_opt, 4)))
  lines(f_hat, type = "l", lwd = 4)
  lines(f_hat_dens, col = "red", lwd = 2)
}


## KEDD's h
# {
#   n = length(x)
#   r <- r_hat(x, n)
#   
#   wiggle <- wiggle_gauss1(x, r)
#   
#   hn_hat_opt <- hn_hat(wiggle, n)
#   hn_hat_opt
#   
#   h = kedd::h.amise(x, kernel = "epanechnikov")$h
#   f_hat <- kern_dens(x, h*sqrt(5)) ## Multiply by sqrt(5) to compare with density()
#   f_hat_dens <- density(x, h, kernel = "epanechnikov", cut = 3*sqrt(5)) ## We multiplied h by sqrt(5) in kern_dens()
#   hist(x, prob=TRUE, xlim = c(min(x), max(x)))
#   lines(f_hat, type = "l", lwd = 4)
#   lines(f_hat_dens, col = "red", lwd = 2)
# }


## Profiling ----

## Standard normal, teoretisk
{
  n = 1000
  x = rnorm(n)
  r <- r_hat(x, n)
}

profvis(wiggle_gauss(x, r), interval = 0.01)
profvis(wiggle_gauss_vec(x, r), interval = 0.01)
profvis(wiggle_gauss_outer(x, r), interval = 0.01)


wiggle_gauss(x, r)
wiggle_gauss_vec(x, r)
wiggle_gauss_outer(x, r)




## Benchmarking ----

{
  n = 1000
  x = rnorm(n)
  r <- r_hat(x, n)
}


## Wiggle

wiggle_bench <- microbenchmark(wiggle_gauss(x, r), wiggle_gauss_vec(x, r), wiggle_gauss_outer(x, r))
levels(wiggle_bench$expr) <- c("nested", "sum", "outer")

wiggle_bench_plot <- autoplot(wiggle_bench) +
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  theme(legend.position = "none") 
wiggle_bench_plot
#labs(title = "Kernel smoother", subtitle = "")


