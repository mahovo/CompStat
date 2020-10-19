# Assignment 4-1 ----
# Optimization: Logistic regression smoothing ----

# *** ----
# SETUP ----

{
  rm(list=ls())
  if(dev.cur() != 1){dev.off()}
  
  # Work. dir. and libraries
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

  # Source functions
  #source("/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/EKSAMEN/3/r/3_functions.R")
}
{
  # save the workspace to the file .RData in the cwd 
  # save.image()
  
  # save the workspace to the file "workspace_HjOpg2.RData" in the cwd 
  # save.image("workspace_HjOpg2.RData")
  
  # Load a workspace into the current session.
  # If you don't specify the path, the cwd is assumed.
  # load("workspace_HjOpg2.RData") 
}


# *** ----
# BENCHMARKING ----
# Tests:


# 1) ====
# ==== Sparse matrices ====


# ---->> Benchmark sparse f(): Sparse matrix created outside benchmark ----
# ---->>> Horse data ----
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 200)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  
  f_sparse_bench_1 <- microbenchmark(
    f(beta, dat2$temp, inner_knots),
    f_sparse(beta, dmat)
  )
  levels(f_sparse_bench_1$expr) = c("H", "H_sparse")
  f_sparse_bench_1
}
num_zero_entries(splineDesign(knots(inner_knots), dat2$temp))
{
  # ---->> Benchmark sparse f(): Sparse matrix created inside benchmark ----
  # ---->>> Horse data ----
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  
  f_sparse_bench_2 <- microbenchmark(
    {for(i in 1:100){f(beta, dat2$temp, inner_knots)}},
    {dmat = design_matrix(beta, dat2$temp, inner_knots, knots);for(i in 1:100){f_sparse(beta, dmat)}}
  )
  levels(f_sparse_bench_2$expr) = c("H", "H_sparse")
  f_sparse_bench_2
}
num_zero_entries(splineDesign(knots(inner_knots), dat2$temp))
{
  # ---->> Plot sparse f() benchmark ----
  # ---->>> Horse data ----
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  f_sparse_bench_3 = NULL
  f_bench_median = NULL
  f_sparse_bench_median = NULL
  lambda = 0.1
  n = 6
  for(i in 1:n){
    f_sparse_bench_3 <- microbenchmark(
      {for(j in 1:(2^(i-1))){f(beta, dat2$temp, inner_knots)}},
      {dmat = design_matrix(beta, dat2$temp, inner_knots, knots);for(j in 1:(2^(i-1))){f_sparse(beta, dmat)}},
      unit = "us"
    )
    f_bench_median[i] = summary(f_sparse_bench_3)$median[1]
    f_sparse_bench_median[i] = summary(f_sparse_bench_3)$median[2]
  }
}
{
  f_sparse_bench_df = data.frame(
    num_iter = 2^(0:(n-1)),
    f_bench_median = f_bench_median,
    f_sparse_bench_median = f_sparse_bench_median
  )
}
{
  plot_f_sparse_bench <- ggplot(data = f_sparse_bench_df, aes(x = num_iter, y = f_bench_median)) +
    geom_point(aes(x = num_iter, y = f_bench_median), size = 3, col = "blue") + 
    geom_line(aes(x = num_iter, y = f_bench_median), col = "blue") +
    geom_point(aes(x = num_iter, y = f_sparse_bench_median), size = 3, col = "red") + 
    geom_line(aes(x = num_iter, y = f_sparse_bench_median), col = "red") +
    scale_x_continuous(trans='log2') +
    scale_y_continuous(trans='log2') +
    labs(x = "Number of iterations", y = "Median runtime (microseconds)") +
    labs(title = 'f() with sparse matrices', subtitle = 'blue = not sparse, red = sparse. 20 basis functions.')
  plot_f_sparse_bench
}
{
  # Write sparse f() benchmark to png
  png('plot_f_sparse_bench.png', width=1800, height=1200, res=300)
  plot_f_sparse_bench
  dev.off()
}



# ----> Benchmark sparse H() ----
# ---->>> Horse data ----
{
  # Test that H and H_sparse produce same output
  {
    inner_knots = dat2$temp
    beta = rnorm(length(inner_knots) + 2)
    lambda = 0.1
    dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
    test1 <- H(beta, dat2$temp, dat2$dead, lambda, inner_knots)
    test2 <- H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots)
    list("H_out" = test1, "H_sparse_out" = test2, "diff" = test1[1,1] - test2[1,1])
  }
}
num_zero_entries(splineDesign(knots(inner_knots), dat2$temp))
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  H_sparse_bench <- microbenchmark(
    H(beta, dat2$temp, dat2$dead, lambda, inner_knots),
    H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots)
  )
  levels(H_sparse_bench$expr) <- c("H", "H_sparse")
  H_sparse_bench
}

# 20 basis functions gives a design matrix with 82% zeros, so good
# balance between sparseness and time for testing.
num_zero_entries(splineDesign(knots(inner_knots), dat2$temp))

# OBS: Plot as function of n
# Result: Sparse is worse, and getting increasingly worse with n 

# Profile:
profvis(for(i in 1:20){H(beta, dat2$temp, dat2$dead, lambda, inner_knots)})
profvis(for(i in 1:20){H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots)})

# Identifies problem in line:
# obj_func = -sum(y * f_val - log(1 + exp(f_val)))

# Check that sparse and original matrices are identical:
sum(tmp_f_val-tmp_f_val_sparse)
# Ok!


# source('Tracer.R')
# H_tracer = tracer(c('f_val'), N = 10)
# H_out =  H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots, cb = H_tracer$trace)
# H_out
# summary(H_tracer)

# Pinpoint the sore spot in the code:
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  tmp_f_val = f(beta, dat2$temp, inner_knots)
  tmp_f_val_sparse = f_sparse(beta, dmat)
}

microbenchmark(
  -sum(dat2$dead * tmp_f_val - log(1 + exp(tmp_f_val))),
  -sum(dat2$dead * tmp_f_val_sparse - log(1 + exp(tmp_f_val_sparse)))
)

# Test sparse matrix in exp function
microbenchmark(
  exp(tmp_f_val),
  exp(tmp_f_val_sparse)
)
# Sparse is slower!

# Test sparse matrix in log function
tmp = 1 + exp(tmp_f_val)
tmp_sparse = 1 + exp(tmp_f_val_sparse)
microbenchmark(
  log(tmp),
  log(tmp_sparse)
)
# Sparse is slower!

# Conclusion:
# f() is faster, when using sparse matrix.
# But f_sparse() outputs a vector as a sparse matrix.
# Using a sparse matrix instead of a vector in H_sparse, is slower than not.
# Over all, using sparse matrix makes Newton (much) slower, because H_sparse is a bottle neck.
# A fix is to wrap last line of f_sparse() in as.vector()

{
  # ---->> Plot sparse H() benchmark ----
  # ---->>> Horse data ----
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  H_sparse_bench = NULL
  H_bench_median = NULL
  H_sparse_bench_median = NULL
  lambda = 0.1
  n = 6
  for(i in 1:n){
    H_sparse_bench <- microbenchmark(
      {for(j in 1:(2^(i-1))){H(beta, dat2$temp, dat2$dead, lambda, inner_knots)}},
      {dmat = design_matrix(beta, dat2$temp, inner_knots, knots);for(j in 1:(2^(i-1))){H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots)}},
      unit = "us"
    )
    H_bench_median[i] = summary(H_sparse_bench)$median[1]
    H_sparse_bench_median[i] = summary(H_sparse_bench)$median[2]
  }
}
{
  H_sparse_bench_df = data.frame(
    num_iter = 2^(0:(n-1)),
    H_bench_median = H_bench_median,
    H_sparse_bench_median = H_sparse_bench_median
  )
  plot_H_sparse_bench <- ggplot(data = H_sparse_bench_df, aes(x = num_iter, y = H_bench_median)) +
    geom_point(aes(x = num_iter, y = H_bench_median), size = 3, col = "blue") + 
    geom_line(aes(x = num_iter, y = H_bench_median), col = "blue") +
    geom_point(aes(x = num_iter, y = H_sparse_bench_median), size = 3, col = "red") + 
    geom_line(aes(x = num_iter, y = H_sparse_bench_median), col = "red") +
    scale_x_continuous(trans='log2') +
    scale_y_continuous(trans='log2') +
    labs(x = "Number of iterations", y = "Median runtime (microseconds)") +
    labs(title = 'H() with sparse matrices', subtitle = 'blue = not sparse, red = sparse. 20 basis functions.')
  plot_H_sparse_bench
}
# Observation: Creating sparse matrix is costly. If you are only going to run the function a few times,
# It is slower to use sparse matrices.

{
  # Write sparse H() benchmark to png
  png('plot_H_sparse_bench.png', width=1800, height=1200, res=300)
  plot_H_sparse_bench
  dev.off()
}


# ----> Benchmark sparse hessian_H() ----
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  hessian_H_sparse_bench <- microbenchmark(
    hessian_H(beta, dat2$temp, dat2$dead, lambda, inner_knots),
    {hes = hessian_H_make_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots);
    hessian_H_sparse(dmat, hes$f_val, hes$p, hes$W, hes$length)}
  )
  hessian_H_sparse_bench
}
{
  tmp_f_val = f_sparse(beta, dmat)
  tmp_p = logit_inv(tmp_f_val)
  tmp_W = diag(as.vector(tmp_p * (1 - tmp_p)))
  tmp_W_sparse = Matrix(diag(as.vector(tmp_p * (1 - tmp_p))), sparse = T)
}
microbenchmark(
  crossprod(dmat, tmp_W %*% dmat),
  crossprod(dmat, tmp_W_sparse %*% dmat)
)
microbenchmark(
  -((-crossprod(dmat, tmp_W %*% dmat) - lambda * pen_mat(inner_knots)) / length(dat2$temp)),
  -((-crossprod(dmat, tmp_W_sparse %*% dmat) - lambda * pen_mat(inner_knots)) / length(dat2$temp))
)
{
  # ---->> Plot sparse H() benchmark ----
  # ---->>> Horse data ----
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  hes_sparse_bench = NULL
  hes_bench_median = NULL
  hes_sparse_bench_median = NULL
  n = 6
  for(i in 1:n){
    hessian_H_sparse_bench <- microbenchmark(
      {for(j in 1:(2^(i-1))){hessian_H(beta, dat2$temp, dat2$dead, lambda, inner_knots)}},
      {dmat = design_matrix(beta, dat2$temp, inner_knots, knots);
      hes = hessian_H_make_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots);
      for(j in 1:(2^(i-1))){hessian_H_sparse(dmat, hes$f_val, hes$p, hes$W, hes$length)}},
      unit = "ms"
    )
    hes_bench_median[i] = summary(hessian_H_sparse_bench)$median[1]
    hes_sparse_bench_median[i] = summary(hessian_H_sparse_bench)$median[2]
  }
}
{
  hes_sparse_bench_df = data.frame(
    num_iter = 2^(0:(n-1)),
    hes_bench_median = hes_bench_median,
    hes_sparse_bench_median = hes_sparse_bench_median
  )
  plot_hes_sparse_bench <- ggplot(data = hes_sparse_bench_df, aes(x = num_iter, y = hes_bench_median)) +
    geom_point(aes(x = num_iter, y = hes_bench_median), size = 3, col = "blue") + 
    geom_line(aes(x = num_iter, y = hes_bench_median), col = "blue") +
    geom_point(aes(x = num_iter, y = hes_sparse_bench_median), size = 3, col = "red") + 
    geom_line(aes(x = num_iter, y = hes_sparse_bench_median), col = "red") +
    scale_x_continuous(trans='log2') +
    scale_y_continuous(trans='log2') +
    labs(x = "Number of iterations", y = "Median runtime (ms)") +
    labs(title = 'hessian_H() with sparse matrices', subtitle = 'blue = not sparse, red = sparse. 20 basis functions.')
  plot_hes_sparse_bench
}
{
  # Write sparse hessian_H() benchmark to png
  png('plot_hes_H_sparse_bench.png', width=1800, height=1200, res=300)
  plot_hes_sparse_bench
  dev.off()
}


# ----> Benchmark sparse pen_mat() ----
# Compare non-sparse with sparse matrix in crossprod
{
  tmp_knots = sort(knots(inner_knots))
  tmp_d = diff(inner_knots)  ## the vector of knot differences; b - a 
  tmp_g_ab = splineDesign(tmp_knots, inner_knots, derivs = 2)
  tmp_g_ab_sparse = Matrix(splineDesign(tmp_knots, inner_knots, derivs = 2), sparse = T)
  tmp_knots_mid = inner_knots[-length(inner_knots)] + tmp_d / 2
  tmp_g_ab_mid = splineDesign(tmp_knots, tmp_knots_mid, derivs = 2)
  tmp_g_ab_mid_sparse = Matrix(splineDesign(tmp_knots, tmp_knots_mid, derivs = 2), sparse = T)
  tmp_g_a = tmp_g_ab[-nrow(tmp_g_ab), ]
  tmp_g_a_sparse = tmp_g_ab_sparse[-nrow(tmp_g_ab_sparse), ]
  tmp_g_b = tmp_g_ab[-1, ]
  tmp_g_b_sparse = tmp_g_ab_sparse[-1, ]
  tmp_g_ab = splineDesign(tmp_knots, inner_knots, derivs = 2)
  tmp_g_ab_sparse = Matrix(splineDesign(tmp_knots, inner_knots, derivs = 2), sparse = T)
}

pen_bench = microbenchmark(
  pen_mat(inner_knots),
  {
    tmp_pen = pen_mat_make_sparse(inner_knots); 
    pen_mat_sparse(tmp_pen$d, tmp_pen$g_ab_mid,  tmp_pen$g_a,  tmp_pen$g_b)
  } 
)

levels(pen_bench$expr) = c("pen_mat", "pen_mat_sparse")
pen_bench
# Sparse is more than 10x slower

# Visualize the pen_mat_sparse matrix
plot_pen = image(pen_mat_sparse(tmp_pen$d, tmp_pen$g_ab_mid,  tmp_pen$g_a,  tmp_pen$g_b))
{
  # Write sparse hessian_H() benchmark to png
  png('plot_pen.png', width=1800, height=1200, res=300)
  plot_pen
  dev.off()
}

str(pen_mat_sparse(tmp_pen$d, tmp_pen$g_ab_mid,  tmp_pen$g_a,  tmp_pen$g_b))

# Number of zeros in g_a matrix
num_zero_entries(tmp_g_a)

pen_bench = microbenchmark(
  crossprod(tmp_g_a,  tmp_g_a),
  crossprod(tmp_g_a_sparse,  tmp_g_a_sparse)
)
levels(pen_bench$expr) = c("g_a", "g_a_sparse")
pen_bench
# Surprisingly, the sparse version is much slower, eventhough the sparse
# matrix is generated outside microbenchmark().
tmp_g_a_sparse

# Check that matrices are identical
sum(tmp_g_a - tmp_g_a_sparse)

# Conclusion: I find no benefit from sparse matrix in pen_mat().



# Test logit_inv with sparse           
microbenchmark(
  logit_inv( f(beta, dat2$temp, inner_knots)),
  logit_inv( f_sparse(beta, dmat))
)
# Sparse is clearly faster
is.vector(tmp_p * (1 - tmp_p))





# --->> Benchmark Newton_sparse ----
# ---->>> Horse data ----
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  # Condition number of hessian_H (but ok, because we never take inverse)
  # kappa(hessian_H(beta, dat2$temp, dat2$dead, lambda = 0.8, inner_knots)) # inf
  newton_sparse_bench <- microbenchmark(
    Newton(beta, dat2$temp, dat2$dead, lambda, inner_knots, H, grad_H, hessian_H),
    Newton_sparse(beta, dat2$temp, dat2$dead, lambda, inner_knots, H_sparse, grad_H_sparse, hessian_H_sparse, knots)
  )
  levels(newton_sparse_bench$expr) <- c("Newton", "Newton_sparse")
  newton_sparse_bench
}
# Plot sparse Newton() benchmark
{
  plot_Newton_sparse_bench <- autoplot(newton_sparse_bench) +
    geom_jitter(position = position_jitter(0.2, 0), 
                aes(color = expr), alpha = 0.4) + 
    aes(fill = I("gray")) + 
    theme(legend.position = "none")
  plot_Newton_sparse_bench
}
{
  # Write sparse Newton() benchmark to png
  png('Newton_sparse_bench.png', width=1800, height=1200, res=300)
  plot_Newton_sparse_bench
  #grid.arrange(plot1, plot2, ncol = 2)
  dev.off()
}

profvis(Newton(beta, dat2$temp, dat2$dead, lambda, inner_knots, H, grad_H, hessian_H))
profvis(Newton_sparse(beta, dat2$temp, dat2$dead, lambda, inner_knots, H_sparse, grad_H_sparse, hessian_H_sparse, knots))

# Conclusion: There are benefits to using sparse matrix for f(). 
# There can be benefits for H() and hessian_H(), if n is large enough. 
# There are no benefits for pen_mat, at least for this data set.
# The over all effect on Newton() is, that there are no benefits from 
# sparse matrices. In fact, it is much slower. I attribute this to the 
# fact, that Newton converges too in too few iterations for the marginal 
# benefit to outweigh the computational cost of generatingsparse matrices.
# (See test 3 below)

{
  
  # *** ----
  # DATA ----
  
  # Data for testing ====
  
  dat1 = read.csv('Horses.csv', header = TRUE)
  x = dat1[ , 1]
  y = dat1[ , 2]
  
  # Removing NAs
  na = is.na(x)
  x = x[!na]
  y = y[!na]
  dat2 = data.frame(temp = x, dead = y)
  
  breaks = c(0, seq(37, 39, 0.125), Inf)
  dat2 = mutate(dat2, cat.temp = cut(temp, breaks = breaks, dig.lab = 4))
  means = data.frame(summarize(group_by(dat2, cat.temp),
                               mean.temp = mean(temp), 
                               sample.p = mean(dead), n = length(temp)))
  
  # Inspecting data ====
  head(dat2)
  plot1 = ggplot(data = dat2) +
    geom_point(aes(x = temp, y = dead)) +
    geom_point(data = means, mapping = aes(x = mean.temp, y = sample.p), col = 'red')
  plot1
}

# *** ----
# BENCHMARKING ----
# Tests:


# 1) ====
# ==== Sparse matrices ====


# ---->> Benchmark sparse f(): Sparse matrix created outside benchmark ----
# ---->>> Horse data ----
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 200)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  
  f_sparse_bench_1 <- microbenchmark(
    f(beta, dat2$temp, inner_knots),
    f_sparse(beta, dmat)
  )
  levels(f_sparse_bench_1$expr) = c("H", "H_sparse")
  f_sparse_bench_1
}
num_zero_entries(splineDesign(knots(inner_knots), dat2$temp))
{
  # ---->> Benchmark sparse f(): Sparse matrix created inside benchmark ----
  # ---->>> Horse data ----
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  
  f_sparse_bench_2 <- microbenchmark(
    {for(i in 1:100){f(beta, dat2$temp, inner_knots)}},
    {dmat = design_matrix(beta, dat2$temp, inner_knots, knots);for(i in 1:100){f_sparse(beta, dmat)}}
  )
  levels(f_sparse_bench_2$expr) = c("H", "H_sparse")
  f_sparse_bench_2
}
num_zero_entries(splineDesign(knots(inner_knots), dat2$temp))
{
  # ---->> Plot sparse f() benchmark ----
  # ---->>> Horse data ----
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  f_sparse_bench_3 = NULL
  f_bench_median = NULL
  f_sparse_bench_median = NULL
  lambda = 0.1
  n = 6
  for(i in 1:n){
    f_sparse_bench_3 <- microbenchmark(
      {for(j in 1:(2^(i-1))){f(beta, dat2$temp, inner_knots)}},
      {dmat = design_matrix(beta, dat2$temp, inner_knots, knots);for(j in 1:(2^(i-1))){f_sparse(beta, dmat)}},
      unit = "us"
    )
    f_bench_median[i] = summary(f_sparse_bench_3)$median[1]
    f_sparse_bench_median[i] = summary(f_sparse_bench_3)$median[2]
  }
}
{
  f_sparse_bench_df = data.frame(
    num_iter = 2^(0:(n-1)),
    f_bench_median = f_bench_median,
    f_sparse_bench_median = f_sparse_bench_median
  )
}
{
  plot_f_sparse_bench <- ggplot(data = f_sparse_bench_df, aes(x = num_iter, y = f_bench_median)) +
    geom_point(aes(x = num_iter, y = f_bench_median), size = 3, col = "blue") + 
    geom_line(aes(x = num_iter, y = f_bench_median), col = "blue") +
    geom_point(aes(x = num_iter, y = f_sparse_bench_median), size = 3, col = "red") + 
    geom_line(aes(x = num_iter, y = f_sparse_bench_median), col = "red") +
    scale_x_continuous(trans='log2') +
    scale_y_continuous(trans='log2') +
    labs(x = "Number of iterations", y = "Median runtime (microseconds)") +
    labs(title = 'f() with sparse matrices', subtitle = 'blue = not sparse, red = sparse. 20 basis functions.')
  plot_f_sparse_bench
}
{
  # Write sparse f() benchmark to png
  png('plot_f_sparse_bench.png', width=1800, height=1200, res=300)
  plot_f_sparse_bench
  dev.off()
}



# ----> Benchmark sparse H() ----
# ---->>> Horse data ----
{
  # Test that H and H_sparse produce same output
  {
    inner_knots = dat2$temp
    beta = rnorm(length(inner_knots) + 2)
    lambda = 0.1
    dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
    test1 <- H(beta, dat2$temp, dat2$dead, lambda, inner_knots)
    test2 <- H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots)
    list("H_out" = test1, "H_sparse_out" = test2, "diff" = test1[1,1] - test2[1,1])
  }
}
num_zero_entries(splineDesign(knots(inner_knots), dat2$temp))
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  H_sparse_bench <- microbenchmark(
    H(beta, dat2$temp, dat2$dead, lambda, inner_knots),
    H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots)
  )
  levels(H_sparse_bench$expr) <- c("H", "H_sparse")
  H_sparse_bench
}

# 20 basis functions gives a design matrix with 82% zeros, so good
# balance between sparseness and time for testing.
num_zero_entries(splineDesign(knots(inner_knots), dat2$temp))

# OBS: Plot as function of n
# Result: Sparse is worse, and getting increasingly worse with n 

# Profile:
profvis(for(i in 1:20){H(beta, dat2$temp, dat2$dead, lambda, inner_knots)})
profvis(for(i in 1:20){H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots)})

# Identifies problem in line:
# obj_func = -sum(y * f_val - log(1 + exp(f_val)))

# Check that sparse and original matrices are identical:
sum(tmp_f_val-tmp_f_val_sparse)
# Ok!


# source('Tracer.R')
# H_tracer = tracer(c('f_val'), N = 10)
# H_out =  H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots, cb = H_tracer$trace)
# H_out
# summary(H_tracer)

# Pinpoint the sore spot in the code:
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  tmp_f_val = f(beta, dat2$temp, inner_knots)
  tmp_f_val_sparse = f_sparse(beta, dmat)
}

microbenchmark(
  -sum(dat2$dead * tmp_f_val - log(1 + exp(tmp_f_val))),
  -sum(dat2$dead * tmp_f_val_sparse - log(1 + exp(tmp_f_val_sparse)))
)

# Test sparse matrix in exp function
microbenchmark(
  exp(tmp_f_val),
  exp(tmp_f_val_sparse)
)
# Sparse is slower!

# Test sparse matrix in log function
tmp = 1 + exp(tmp_f_val)
tmp_sparse = 1 + exp(tmp_f_val_sparse)
microbenchmark(
  log(tmp),
  log(tmp_sparse)
)
# Sparse is slower!

# Conclusion:
# f() is faster, when using sparse matrix.
# But f_sparse() outputs a vector as a sparse matrix.
# Using a sparse matrix instead of a vector in H_sparse, is slower than not.
# Over all, using sparse matrix makes Newton (much) slower, because H_sparse is a bottle neck.
# A fix is to wrap last line of f_sparse() in as.vector()

{
  # ---->> Plot sparse H() benchmark ----
  # ---->>> Horse data ----
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  H_sparse_bench = NULL
  H_bench_median = NULL
  H_sparse_bench_median = NULL
  lambda = 0.1
  n = 6
  for(i in 1:n){
    H_sparse_bench <- microbenchmark(
      {for(j in 1:(2^(i-1))){H(beta, dat2$temp, dat2$dead, lambda, inner_knots)}},
      {dmat = design_matrix(beta, dat2$temp, inner_knots, knots);for(j in 1:(2^(i-1))){H_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots)}},
      unit = "us"
    )
    H_bench_median[i] = summary(H_sparse_bench)$median[1]
    H_sparse_bench_median[i] = summary(H_sparse_bench)$median[2]
  }
}
{
  H_sparse_bench_df = data.frame(
    num_iter = 2^(0:(n-1)),
    H_bench_median = H_bench_median,
    H_sparse_bench_median = H_sparse_bench_median
  )
  plot_H_sparse_bench <- ggplot(data = H_sparse_bench_df, aes(x = num_iter, y = H_bench_median)) +
    geom_point(aes(x = num_iter, y = H_bench_median), size = 3, col = "blue") + 
    geom_line(aes(x = num_iter, y = H_bench_median), col = "blue") +
    geom_point(aes(x = num_iter, y = H_sparse_bench_median), size = 3, col = "red") + 
    geom_line(aes(x = num_iter, y = H_sparse_bench_median), col = "red") +
    scale_x_continuous(trans='log2') +
    scale_y_continuous(trans='log2') +
    labs(x = "Number of iterations", y = "Median runtime (microseconds)") +
    labs(title = 'H() with sparse matrices', subtitle = 'blue = not sparse, red = sparse. 20 basis functions.')
  plot_H_sparse_bench
}
# Observation: Creating sparse matrix is costly. If you are only going to run the function a few times,
# It is slower to use sparse matrices.

{
  # Write sparse H() benchmark to png
  png('plot_H_sparse_bench.png', width=1800, height=1200, res=300)
  plot_H_sparse_bench
  dev.off()
}


# ----> Benchmark sparse hessian_H() ----
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  hessian_H_sparse_bench <- microbenchmark(
    hessian_H(beta, dat2$temp, dat2$dead, lambda, inner_knots),
    {hes = hessian_H_make_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots);
    hessian_H_sparse(dmat, hes$f_val, hes$p, hes$W, hes$length)}
  )
  hessian_H_sparse_bench
}
{
  tmp_f_val = f_sparse(beta, dmat)
  tmp_p = logit_inv(tmp_f_val)
  tmp_W = diag(as.vector(tmp_p * (1 - tmp_p)))
  tmp_W_sparse = Matrix(diag(as.vector(tmp_p * (1 - tmp_p))), sparse = T)
}
microbenchmark(
  crossprod(dmat, tmp_W %*% dmat),
  crossprod(dmat, tmp_W_sparse %*% dmat)
)
microbenchmark(
  -((-crossprod(dmat, tmp_W %*% dmat) - lambda * pen_mat(inner_knots)) / length(dat2$temp)),
  -((-crossprod(dmat, tmp_W_sparse %*% dmat) - lambda * pen_mat(inner_knots)) / length(dat2$temp))
)
{
  # ---->> Plot sparse H() benchmark ----
  # ---->>> Horse data ----
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  hes_sparse_bench = NULL
  hes_bench_median = NULL
  hes_sparse_bench_median = NULL
  n = 6
  for(i in 1:n){
    hessian_H_sparse_bench <- microbenchmark(
      {for(j in 1:(2^(i-1))){hessian_H(beta, dat2$temp, dat2$dead, lambda, inner_knots)}},
      {dmat = design_matrix(beta, dat2$temp, inner_knots, knots);
      hes = hessian_H_make_sparse(beta, dat2$temp, dat2$dead, dmat, lambda, inner_knots);
      for(j in 1:(2^(i-1))){hessian_H_sparse(dmat, hes$f_val, hes$p, hes$W, hes$length)}},
      unit = "ms"
    )
    hes_bench_median[i] = summary(hessian_H_sparse_bench)$median[1]
    hes_sparse_bench_median[i] = summary(hessian_H_sparse_bench)$median[2]
  }
}
{
  hes_sparse_bench_df = data.frame(
    num_iter = 2^(0:(n-1)),
    hes_bench_median = hes_bench_median,
    hes_sparse_bench_median = hes_sparse_bench_median
  )
  plot_hes_sparse_bench <- ggplot(data = hes_sparse_bench_df, aes(x = num_iter, y = hes_bench_median)) +
    geom_point(aes(x = num_iter, y = hes_bench_median), size = 3, col = "blue") + 
    geom_line(aes(x = num_iter, y = hes_bench_median), col = "blue") +
    geom_point(aes(x = num_iter, y = hes_sparse_bench_median), size = 3, col = "red") + 
    geom_line(aes(x = num_iter, y = hes_sparse_bench_median), col = "red") +
    scale_x_continuous(trans='log2') +
    scale_y_continuous(trans='log2') +
    labs(x = "Number of iterations", y = "Median runtime (ms)") +
    labs(title = 'hessian_H() with sparse matrices', subtitle = 'blue = not sparse, red = sparse. 20 basis functions.')
  plot_hes_sparse_bench
}
{
  # Write sparse hessian_H() benchmark to png
  png('plot_hes_H_sparse_bench.png', width=1800, height=1200, res=300)
  plot_hes_sparse_bench
  dev.off()
}


# ----> Benchmark sparse pen_mat() ----
# Compare non-sparse with sparse matrix in crossprod
{
  tmp_knots = sort(knots(inner_knots))
  tmp_d = diff(inner_knots)  ## the vector of knot differences; b - a 
  tmp_g_ab = splineDesign(tmp_knots, inner_knots, derivs = 2)
  tmp_g_ab_sparse = Matrix(splineDesign(tmp_knots, inner_knots, derivs = 2), sparse = T)
  tmp_knots_mid = inner_knots[-length(inner_knots)] + tmp_d / 2
  tmp_g_ab_mid = splineDesign(tmp_knots, tmp_knots_mid, derivs = 2)
  tmp_g_ab_mid_sparse = Matrix(splineDesign(tmp_knots, tmp_knots_mid, derivs = 2), sparse = T)
  tmp_g_a = tmp_g_ab[-nrow(tmp_g_ab), ]
  tmp_g_a_sparse = tmp_g_ab_sparse[-nrow(tmp_g_ab_sparse), ]
  tmp_g_b = tmp_g_ab[-1, ]
  tmp_g_b_sparse = tmp_g_ab_sparse[-1, ]
  tmp_g_ab = splineDesign(tmp_knots, inner_knots, derivs = 2)
  tmp_g_ab_sparse = Matrix(splineDesign(tmp_knots, inner_knots, derivs = 2), sparse = T)
}

pen_bench = microbenchmark(
  pen_mat(inner_knots),
  {
    tmp_pen = pen_mat_make_sparse(inner_knots); 
    pen_mat_sparse(tmp_pen$d, tmp_pen$g_ab_mid,  tmp_pen$g_a,  tmp_pen$g_b)
  } 
)

levels(pen_bench$expr) = c("pen_mat", "pen_mat_sparse")
pen_bench
# Sparse is more than 10x slower

# Visualize the pen_mat_sparse matrix
plot_pen = image(pen_mat_sparse(tmp_pen$d, tmp_pen$g_ab_mid,  tmp_pen$g_a,  tmp_pen$g_b))
{
  # Write sparse hessian_H() benchmark to png
  png('plot_pen.png', width=1800, height=1200, res=300)
  plot_pen
  dev.off()
}

str(pen_mat_sparse(tmp_pen$d, tmp_pen$g_ab_mid,  tmp_pen$g_a,  tmp_pen$g_b))

# Number of zeros in g_a matrix
num_zero_entries(tmp_g_a)

pen_bench = microbenchmark(
  crossprod(tmp_g_a,  tmp_g_a),
  crossprod(tmp_g_a_sparse,  tmp_g_a_sparse)
)
levels(pen_bench$expr) = c("g_a", "g_a_sparse")
pen_bench
# Surprisingly, the sparse version is much slower, eventhough the sparse
# matrix is generated outside microbenchmark().
tmp_g_a_sparse

# Check that matrices are identical
sum(tmp_g_a - tmp_g_a_sparse)

# Conclusion: I find no benefit from sparse matrix in pen_mat().



# Test logit_inv with sparse           
microbenchmark(
  logit_inv( f(beta, dat2$temp, inner_knots)),
  logit_inv( f_sparse(beta, dmat))
)
# Sparse is clearly faster
is.vector(tmp_p * (1 - tmp_p))





# --->> Benchmark Newton_sparse ----
# ---->>> Horse data ----
{
  rg = range(dat2$temp)
  inner_knots = seq(rg[1], rg[2], length.out = 20)
  beta = rnorm(length(inner_knots) + 2)
  lambda = 0.1
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  # Condition number of hessian_H (but ok, because we never take inverse)
  # kappa(hessian_H(beta, dat2$temp, dat2$dead, lambda = 0.8, inner_knots)) # inf
  newton_sparse_bench <- microbenchmark(
    Newton(beta, dat2$temp, dat2$dead, lambda, inner_knots, H, grad_H, hessian_H),
    Newton_sparse(beta, dat2$temp, dat2$dead, lambda, inner_knots, H_sparse, grad_H_sparse, hessian_H_sparse, knots)
  )
  levels(newton_sparse_bench$expr) <- c("Newton", "Newton_sparse")
  newton_sparse_bench
}
# Plot sparse Newton() benchmark
{
  plot_Newton_sparse_bench <- autoplot(newton_sparse_bench) +
    geom_jitter(position = position_jitter(0.2, 0), 
                aes(color = expr), alpha = 0.4) + 
    aes(fill = I("gray")) + 
    theme(legend.position = "none")
  plot_Newton_sparse_bench
}
{
  # Write sparse Newton() benchmark to png
  png('Newton_sparse_bench.png', width=1800, height=1200, res=300)
  plot_Newton_sparse_bench
  #grid.arrange(plot1, plot2, ncol = 2)
  dev.off()
}

profvis(Newton(beta, dat2$temp, dat2$dead, lambda, inner_knots, H, grad_H, hessian_H))
profvis(Newton_sparse(beta, dat2$temp, dat2$dead, lambda, inner_knots, H_sparse, grad_H_sparse, hessian_H_sparse, knots))

# Conclusion: There are benefits to using sparse matrix for f(). 
# There can be benefits for H() and hessian_H(), if n is large enough. 
# There are no benefits for pen_mat, at least for this data set.
# The over all effect on Newton() is, that there are no benefits from 
# sparse matrices. In fact, it is much slower. I attribute this to the 
# fact, that Newton converges too in too few iterations for the marginal 
# benefit to outweigh the computational cost of generatingsparse matrices.
# (See test 3 below)





# 2) ====
# ==== Logit ====
logit_inv_bench = microbenchmark(
  1 / (1 + exp(-x)),
  exp(x) / (1 + exp(x))
)
logit_inv_bench
# Ratio
summary(logit_inv_bench)$median[1]/summary(logit_inv_bench)$median[2]
# Conclusion: The first expression is almost twice as fast!
# However, not a bottleneck in algorithms where large matrices are bottlenecks.


# 3) ====
# ==== Number of iterations for convergence ====
rg = range(dat2$temp)
inner_knots = seq(rg[1], rg[2], length.out = 20)
beta = rnorm(length(inner_knots) + 2)
lambda = 0.1
tmp = Newton_sparse(beta, dat2$temp, dat2$dead, lambda, inner_knots, H_sparse, grad_H_sparse, hessian_H_sparse, knots)
tmp$n
# RESULT: n = 8


