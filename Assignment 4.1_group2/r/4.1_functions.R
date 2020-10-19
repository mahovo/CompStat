
## Functions for optimization task ====

## Handling boundary issues with splines
knots = function(inner_knots) {
  c(rep(range(inner_knots), 3), inner_knots)
}

## Inverse logit for getting probabilities
## logit_inv = function(x) {
##   exp(x) / (1 + exp(x))
## }
logit_inv = function(x) {
  1 / (1 + exp(-x))
}


## Smooth estimate of logit conditional probability
f = function(par, x, inner_knots) {
  if(length(par) != length(inner_knots) + 2) {
    stop(paste('\ndimension of parameter vector (', length(par), ') and the spline design matrix (', length(inner_knots) + 2), ') does not fit.', sep = '')
  }
  X = splineDesign(knots(inner_knots), x)
  X %*% par
}

### Test:
# {
#   inner_knots = dat2$temp
#   beta = rnorm(length(inner_knots) + 2)
#   f(beta, dat2$temp, dat2$temp)
#   
#   inner_knots = dat2$temp
#   X = splineDesign(knots(inner_knots), dat2$temp)
#   dim(X)[1] * dim(X)[2] - sum(X != 0)
# }

## Penalty matrix
## CSwR, 3.5.2
pen_mat = function(inner_knots) {
  knots = sort(knots(inner_knots))
  d = diff(inner_knots)  ## the vector of knot differences; b - a 
  g_ab = splineDesign(knots, inner_knots, derivs = 2) 
  knots_mid = inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid = splineDesign(knots, knots_mid, derivs = 2)
  g_a = g_ab[-nrow(g_ab), ]
  g_b = g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 
}
### Test:
# {
#   rg = range(dat2$temp)
#   inner_knots = seq(rg[1], rg[2], length.out = 5)
#   lambda = 0.2
#   beta = rnorm(length(inner_knots) + 2)
#   dim(pen_mat(inner_knots))
#   t(beta) %*% pen_mat(inner_knots) %*% beta
# }

## H: Penalized log-likelihood criterion to minimize
H = function(par, z, y, lambda, inner_knots) {
  f_val = f(par, z, inner_knots)
  obj_func = -sum(y * f_val - log(1 + exp(f_val)))
  penalty = lambda * crossprod(par, pen_mat(inner_knots) %*% par)
  (obj_func + penalty) / length(z)
}
### Test
# {
#   inner_knots = dat2$temp
#   beta = rep(0, length(inner_knots) + 2)
#   lambda = 0.2
#   H(beta, dat2$temp, dat2$dead, lambda, inner_knots)
# }

## Derivatives of H
grad_H = function(par, z, y, lambda, inner_knots) {
  phi = splineDesign(knots(inner_knots), z)
  f_val = f(par, z, inner_knots)
  p = logit_inv(f_val)
  -((crossprod(phi, y - p) - lambda * pen_mat(inner_knots) %*% par) / length(z))
}
hessian_H = function(par, z, y, lambda, inner_knots) {
  phi = splineDesign(knots(inner_knots), z)
  f_val = f(par, z, inner_knots)
  p = logit_inv(f_val)
  W = diag(as.vector(p * (1 - p)))
  -((-crossprod(phi, W %*% phi) - lambda * pen_mat(inner_knots)) / length(z))
}
### Test
# {
#   rg = range(dat2$temp)
#   inner_knots = c(rg[1], (rg[2] - rg[1]) / 2, rg[2])
#   lambda = 0.2
#   beta = rep(0, length(inner_knots) + 2)
#   grad_H(beta, dat2$temp, dat2$dead, lambda, inner_knots)
#   beta = beta - grad_H(beta, dat2$temp, dat2$dead, lambda, inner_knots) * 0.2
#   grad_H(beta, dat2$temp, dat2$dead, lambda, inner_knots)
#   dim(hessian_H(beta, dat2$temp, dat2$dead, lambda, inner_knots))
#   eigen(hessian_H(beta, dat2$temp, dat2$dead, lambda, inner_knots))$values
# }





## Optimization algorithms ====
##
## ----> Newton optimization algorithm  ----

Newton = function(
  par, x, y, lambda = 0.2, inner_knots, H, grad_H, Hessian_H, maxiter = 50, 
  d = 0.8, c = 0.1, gamma0 = 1, epsilon = 1e-6, stop = 'func', cb = NULL) {
  n = 1
  
  ## Containers for function and paramter values
  f_val = rep(NA, maxiter)        
  par_val = matrix(NA, nrow = length(par), ncol = maxiter)
  par_val[ , 1] = par
  
  ## Descent
  repeat {
    f_val[n] = H(par, x, y, lambda, inner_knots)
    grad = grad_H(par, x, y, lambda, inner_knots)
    
    ## Stopping criteria
    if (n > maxiter) {
      break
    } else if (stop == 'func') {
      if (n > 2 && f_val[n - 1] - f_val[n] < epsilon * (f_val[n] + epsilon)) break
    } else if (stop == 'grad') {
      if (sum(grad^2) <= epsilon || n > maxiter) break
    } else if (stop == 'par') {
      if (n > 2 && sum((par_val[ , n - 1] - par_val[ , n - 2])^2) <= 
          epsilon * ((sum(par_val[ , n - 1])^2) + epsilon)) break
    }
    
    ## Descent direction and step
    Hessian = Hessian_H(beta, x, y, lambda, inner_knots)
    rho = -drop(solve(Hessian, grad)) 
    gamma = gamma0
    par1 = par + gamma * rho
    h_prime = t(grad) %*% rho 
    
    ## Tracer
    if(!is.null(cb)) cb()
    
    ## Backtracking
    while(H(par1, x, y, lambda, inner_knots) > f_val[n] +  c * gamma * h_prime) { 
      gamma = d * gamma 
      par1 = par + gamma * rho
    }
    par_val[ , n] = par1
    par = par1 
    n = n + 1
  }
  list(par = par, f_val = f_val[!is.na(f_val)], par_val = par_val[!is.na(par_val)])
}

## Sparse matrices ====

num_zero_entries = function(mat){
  x = sum(mat == 0) # Number of entries equal to 
  y = dim(mat)[[1]]*dim(mat)[[2]] # Total number of entries in matrix
  z = x/y
  list("numZero" = x, "totEntries" = y, "ratio" = z)
}

## Generate design matrix
design_matrix = function(par, x, inner_knots, knots){
  if(length(par) != length(inner_knots) + 2) {
    stop('\ndimension of parameter vector and the spline design matrix does not fit.')
  }
  Matrix(splineDesign(knots(inner_knots), x), sparse = TRUE)
}

## dmat is the design matrix
f_sparse = function(par, dmat) {
  # crossprod(dmat, par)
  as.vector(dmat %*% par) # Important to make into vector. 
  # If output is sparse matrix, H_sparse() will be slower than H()!
}
H_sparse = function(par, x, y, dmat, lambda, inner_knots, cb = NULL) {
  f_val = f_sparse(par, dmat)
  obj_func = -sum(y * f_val - log(1 + exp(f_val)))
  penalty = lambda * crossprod(par, pen_mat(inner_knots) %*% par)
  (obj_func + penalty) / length(x)
}
pen_mat_make_sparse = function(inner_knots) {
  knots = sort(knots(inner_knots))
  d = diff(inner_knots)  ## the vector of knot differences; b - a 
  g_ab = Matrix(splineDesign(knots, inner_knots, derivs = 2), sparse = T)
  knots_mid = inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid = Matrix(splineDesign(knots, knots_mid, derivs = 2), sparse = T)
  g_a = g_ab[-nrow(g_ab), ]
  g_b = g_ab[-1, ]
  list("d" = d, "g_ab_mid" = g_ab_mid, "g_a" = g_a, "g_b" = g_b)
}
pen_mat_sparse = function(d, g_ab_mid, g_a, g_b) {
  (crossprod(d * g_a,  g_a) + 4 * crossprod(d * g_ab_mid, g_ab_mid) + 
     crossprod(d * g_b, g_b)) / 6
}
grad_H_sparse = function(par, z, y, phi, lambda, inner_knots) {
  # phi = Matrix(splineDesign(knots(inner_knots), z), sparse = TRUE)
  f_val = f_sparse(par, dmat)
  p = logit_inv(f_val)
  -((crossprod(phi, y - p) - lambda * pen_mat(inner_knots) %*% par) / length(z))
}
hessian_H_make_sparse = function(par, z, y, phi, lambda, inner_knots) {
  # phi = splineDesign(knots(inner_knots), z)
  # f_val = f(par, z, inner_knots)
  f_val = f_sparse(par, phi)
  p = logit_inv(f_val)
  W = Matrix(diag(as.vector(p * (1 - p))), sparse = T)
  # W = diag(as.vector(p * (1 - p)))
  length = length(z)
  list("f_val" = f_val, "p" = p, "W" = W, "length" = length)
}
hessian_H_sparse = function(phi, f_val, p, W, length){
  -((-crossprod(phi, W %*% phi) - lambda * pen_mat(inner_knots)) / length)
}

Newton_sparse = function(
  par, x, y, lambda = 0.2, inner_knots, H_sparse, grad_H_sparse, hessian_H_sparse, knots, maxiter = 50, 
  d = 0.8, c = 0.1, gamma0 = 1, epsilon = 1e-6, stop = 'func', cb = NULL) {
  n = 1
  
  # Containers for function and paramter values
  f_val = rep(NA, maxiter)        
  par_val = matrix(NA, nrow = length(par), ncol = maxiter)
  par_val[ , 1] = par
  
  # Sparse design matrix
  dmat = design_matrix(beta, dat2$temp, inner_knots, knots)
  
  # Descent
  repeat {
    f_val[n] = H_sparse(par, x, y, dmat, lambda, inner_knots)
    grad = grad_H_sparse(par, x, y, dmat, lambda, inner_knots)
    
    # Stopping criteria
    if (n > maxiter) {
      break
    } else if (stop == 'func') {
      if (n > 2 && f_val[n - 1] - f_val[n] < epsilon * (f_val[n] + epsilon)) break
    } else if (stop == 'grad') {
      if (sum(grad^2) <= epsilon || n > maxiter) break
    } else if (stop == 'par') {
      if (n > 2 && sum((par_val[ , n - 1] - par_val[ , n - 2])^2) <= 
          epsilon * ((sum(par_val[ , n - 1])^2) + epsilon)) break
    }
    
    ## Descent direction and step
    ## Hessian = hessian_H_sparse(beta, x, y, dmat, lambda, inner_knots)
    hes = hessian_H_make_sparse(beta, x, y, dmat, lambda, inner_knots)
    Hessian = hessian_H_sparse(dmat, hes$f_val, hes$p, hes$W, hes$length)
    rho = -drop(solve(Hessian, grad)) 
    gamma = gamma0
    par1 = par + gamma * rho
    h_prime = t(grad) %*% rho 
    
    ## Tracer
    if(!is.null(cb)) cb()
    
    ## Backtracking
    while((H_sparse(par1, x, y, dmat, lambda, inner_knots))[1, 1] > (f_val[n] +  c * gamma * h_prime)[1, 1]) { 
      gamma = d * gamma 
      par1 = par + gamma * rho
    }
    par_val[ , n] = par1
    par = par1 
    n = n + 1
  }
  list(par = par, f_val = f_val[!is.na(f_val)], par_val = par_val[!is.na(par_val)], n = n)
}

#####################################################################
## Fra assignment 1.1 (brugbart??) ====
#####################################################################
smoother_1 = function(y, lambda, Phi, Omega) {
  Phi %*% solve(
    crossprod(Phi) + lambda * Omega,
    t(Phi)
    %*% y)
}



## 2) smoother_dr (smoother with Demmler-Reinsch basis via SVD)
## To calculate the optimal lambda, we need to create a lot of splines.
## For this we want a more efficient implementation of the spline smoother.


#source('/Users/mhvpbp13/Library/Mobile Documents/com~apple~CloudDocs/CBS/cbs/Semester k1/CompStat/Hjemmeopgaver/Hjemmeopgave 1/Tracer.R')
#smoother_tracer = tracer(c('test'), N = 1)

# smoother_dr = function(y, lambda, U_tilde, Gamma, Omega) {
#   beta_hat = crossprod(U_tilde, y)
#   beta_hat_lambda = beta_hat * solve((1+(lambda * Gamma)))
#   U_tilde %*% beta_hat_lambda
# }
# smoother_dr(Y, lambda, DR_out$U_tilde, DR_out$Gamma, Omega, smoother_tracer$trace)


# smoother_dr = function(y, lambda, U_tilde, Gamma, Omega) {
#   beta_hat = crossprod(U_tilde, y)
#   beta_hat_lambda = beta_hat # Place holder
#   for(i in seq_along(beta_hat_lambda)){
#     beta_hat_lambda[i] = beta_hat[i]/(1+(lambda * Gamma[i]))
#   }
#   U_tilde %*% beta_hat_lambda
# }

smoother_dr = function(y, lambda, U_tilde, Gamma) {
  beta_hat = crossprod(U_tilde, y)
  beta_hat_lambda = beta_hat/(1+(lambda * Gamma))
  U_tilde %*% beta_hat_lambda
}

## -> Demmler-Reinsch basis ----
## Output:
## 1) U_tilde (the Demmler-Reinsch basis)
## 2) Gamma
## 3) Omega_tilde
## 4) W
DR = function(X, inner_knots, Omega){
  Phi <- splineDesign(c(rep(range(inner_knots), 3), inner_knots), X)
  Phi_svd = svd(Phi)
  Omega_tilde = t(crossprod(Phi_svd$v, Omega %*% Phi_svd$v)) / Phi_svd$d
  Omega_tilde = t(Omega_tilde) / Phi_svd$d
  Omega_tilde_svd = svd(Omega_tilde)  
  U_tilde = Phi_svd$u %*% Omega_tilde_svd$u
  W = Omega_tilde_svd$u
  # abs() only affects numeric noise which happens to be negative.
  # All eigenvalues will be non-negative.
  Gamma = abs(diag(crossprod(W, Omega_tilde %*% W)))
  list(U_tilde = U_tilde, Gamma = Gamma, Omega_tilde = Omega_tilde, W = W)
}


## Detect zero entries in matrix ----

## Output:
## 1) Number of zero entries
## 2) Total number of entries
## 3) Ratio of number of zero entries to total number of entries

num_zero_entries = function(mat){
  x = sum(mat == 0) # Number of entries equal to 
  y = dim(mat)[[1]]*dim(mat)[[2]] # Total number of entries in matrix
  z = x/y
  list(x, y, z)
}

## Lambda-tuner ----
## Lambda-tuner using LOOCV
## Compute B-splines with splineDesign(). Riemann-sum approximations, Simpson's rule with breakpoints in the knots.
## Alternative: bsplinepen()
## CSwR 3.5.2 (slide 23-25)

## GCV ----
## CSwR, 3.5.2 (slide 23)

gcv <- function(lambda, y) {
  S <- Phi %*% solve(crossprod(Phi) + lambda * Omega, t(Phi))
  df <- sum(diag(S))  ## The trace of the smoother matrix
  sum(((y - S %*% y) / (1 - df / length(y)))^2, na.rm = TRUE) 
}
