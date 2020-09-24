
## Epanechnikov kernel
epan <- function(x){
  val <- 0.75*(1 - x^2)
  val * (abs(x) < 1)
}

## Kernel density
# kern_dens <- function(x, h, m = 512) {
#   rg <- range(x)
#   n <- length(x)
#   # xx is equivalent to grid points in 'density()'
#   xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m) ## cut=3 default i density()
#   y <- numeric(m) # The evaluations, initialized as a vector of zeros
#   # The actual computation is done using nested for-loops. The outer loop
#   # is over the grid points, and the inner loop is over the data points.
#   for (i in seq_along(xx))
#     for (j in seq_along(x))
#       y[i] <- y[i] + epan((xx[i] - x[j])/h)
#   y <- y / (n*h)
#   list(x = xx, y = y)
# }

# kern_dens <- function(x, h, m = 512, kernel = epan) {
#   rg <- range(x)
#   n <- length(x)
#   ## xx is equivalent to grid points in 'density()'
#   xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m) ## cut=3 default i density()
#   y <- numeric(m) # The evaluations, initialized as a vector of zeros
#   ## The actual computation is done using nested for-loops. The outer loop
#   ## is over the grid points, and the inner loop is over the data points.
#   for (i in seq_along(xx))
#     for (j in seq_along(x))
#       y[i] <- y[i] + kernel((xx[i] - x[j])/h)
#   y <- y / (n*h)
#   list(x = xx, y = y)
# }

kern_dens <- function(x, h, m = 512, kernel = epan) {
  rg <- range(x)
  n <- length(x)
  grid_points <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m) ## cut=3 default in density()
  y <- numeric(m)
  for (i in seq_along(grid_points)) {
    y[i] <- sum(kernel((grid_points[i] - x)/h))
  }
  y <- y / (n*h)
  list(x = grid_points, y = y)
}



## 2.3.3
r_hat <- function(x, n) {
  sigma_hat <- sd(x)
  IQR <- quantile(x, 0.75) - quantile(x, 0.25)
  IQR_theoretical <- qnorm(0.75) - qnorm(0.25)
  #sigma_tilde <- min(sigma_hat, IQR/1.34) # 1.34 by tradition, rather than 1.35
  sigma_tilde <- min(sigma_hat, IQR/IQR_theoretical)
  #0.9 * sigma_tilde * n^(-0.2) ## Silverman's rule
  (4/3)^(1/5) * sigma_tilde * n^(-0.2) ## (4/3)^(1/5) = 1.059224
}

## H double prime, gaussian kernel
# H_dp <- function(x) {
#   exp(-0.5 * x^2) * (x^2 - 1) / sqrt(2 * pi)
# }

## 2.3.3
## "wiggle" here is the 2-norm of f_tilde
## Version 1, Clive R. Loader: "Bandwidth Selection: Classical or Plug-In?"
# wiggle_gauss1 <- function(x, r) {
#   wiggle = 0
#   for(i in seq_along(x)){
#       for(j in seq_along(x)) {
#         z <- (x[i] - x[j]) / (sqrt(2) * r)
#         wiggle <- wiggle + exp(-0.5 * z^2) * (z^4 - 6 * z^2 + 3)/sqrt(2*pi)
#       }
#   }
#   wiggle <- wiggle / (n^2 * (sqrt(2)*r)^5)
#   wiggle
# }

# wiggle_gauss2 <- function(x, r) {
#   wiggle = 0
#   z = numeric(0L)
#   count = 1
#   for(i in seq_along(x)){
#     for(j in seq_along(x)) {
#       z[count] <- (x[i] - x[j]) / (sqrt(2) * r)
#       count = count + 1
#     }
#   }
#   w <- sapply(z , function(z){exp(-0.5 * z^2) * (z^4 - 6 * z^2 + 3)/sqrt(2*pi)})
#   #wiggle <- sum(exp(-0.5 * z^2) * (z^4 - 6 * z^2 + 3)/sqrt(2*pi))
#   wiggle <- sum(w) / (n^2 * (sqrt(2)*r)^5)
#   wiggle
# }

# wiggle_gauss3 <- function(x, r) {
#   wiggle = 0
#   for(i in seq_along(x)){
#     for(j in seq_along(x)) {
#       z <- (x[i] - x[j])/r
#       wiggle <- wiggle + exp(-0.25 * z^2) * (0.25*z^4 - 3 * z^2 + 3)/sqrt(2*pi)
#     }
#   }
#   wiggle <- wiggle / (n^2 * (sqrt(2)*r)^5)
#   wiggle
# }

## > Nested ----
wiggle_gauss <- function(x, r) {
  wiggle = 0
  c1 = -1/(4*r^2)
  c2 = 1/(sqrt(2)*r)
  for(i in seq_along(x)){
    for(j in seq_along(x)) {
      z <- (x[i] - x[j])/c2
      wiggle <- wiggle + exp(c1 * z^2) * ((c2*z)^4 - 6 * (c2*z)^2 + 3)
    }
  }
  wiggle <- wiggle / (8 * n^2 * r^5 * sqrt(pi))
  wiggle
}

## > sum ----
wiggle_gauss_vec <- function(x, r) {
  wiggle = 0
  c1 = -1/(4*r^2)
  c2 = 1/(sqrt(2)*r)
  for(i in seq_along(x)){
      z <- (x[i] - x)/c2
      wiggle <- wiggle + sum(exp(c1 * z^2) * ((c2*z)^4 - 6 * (c2*z)^2 + 3))
  }
  wiggle <- wiggle / (8 * n^2 * r^5 * sqrt(pi))
  wiggle
}

## > outer ----
wiggle_gauss_outer <- function(x, r) {
  c1 = -1/(4*r^2)
  c2 = 1/(sqrt(2)*r)
    wiggle <- outer(x, x, function(ww, w){
      z <- (ww - w)/c2
      exp(c1 * z^2) * ((c2*z)^4 - 6 * (c2*z)^2 + 3)
    })
  sum(wiggle) / (8 * n^2 * r^5 * sqrt(pi))
}


## 2.3.3
## (2.3)
# hn_hat <- function(wiggle, n) {
#   K_2norm <- 0.5 / sqrt(pi)
#   (K_2norm / wiggle)^(0.2) * n^(-0.2) ## sigma_K^4 = 1
# }

hn_hat <- function(wiggle, n) {
  K_2norm <-0.6 ## See exercise 2.2
  (25 * K_2norm / wiggle)^(0.2) * n^(-0.2) ## sigma_K^4 = 1
}








## Slides:
## http://web.math.ku.dk/~richard/courses/compstat2020/Nearest_neighbor_static.html#7
# run_mean <- function(y, k) {
#   n <- length(y)
#   m <- floor((k - 1) / 2)
#   k <- 2 * m + 1           # Ensures k to be odd and m = (k - 1) / 2
#   y <- y / k
#   s <- rep(NA, n)
#   s[m + 1] <- sum(y[1:k])
#   for(i in (m + 1):(n - m - 1)) 
#     s[i + 1] <- s[i] - y[i - m] + y[i + 1 + m]
#   s
# }

## Slides
## http://web.math.ku.dk/~richard/courses/compstat2020/Nearest_neighbor_static.html#13
# loocv <- function(k, y) {
#   f_hat <- run_mean(y, k)
#   mean(((y - f_hat) / (1 - 1/k))^2, na.rm = TRUE) 
# }
