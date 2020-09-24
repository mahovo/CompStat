##########################################################
## Saif Babrak                                          ##
## Computational Statistics                             ##
## Functions used in assignment 1 - density estimation  ##
##########################################################

EpanechnikovKernel <- function(x) {
  ifelse(abs(x) < 1, 3/4 * (1 - x**2), 0)
}


KernelDensity <- function(x, h, m = 512) {
  ## This is an implementation of the function 'KernelDensity' that computes 
  ## evaluations of Epanechnikov kernel density estimates in a grid of points.
  
  n <- length(x)
  rg <- range(x)
  
  GridPoints <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m) 
  
  const = h * n
  for (i in seq_along(GridPoints))
    for(j in seq_along(x)) {
      y[i] <-  sum(EpanechnikovKernel((GridPoints[i] - x[j]) / h)) / const
    }
  list(GridPoints = GridPoints, y = y)
}


VectorizedKernelDensity <- function(x, h, m = 512) {
  # This is an implementation of the function 'VectorizedKernelDensity' that computes 
  # evaluations of Epanechnikov kernel density estimates in a grid of points.
  
  n <- length(x)
  rg <- range(x)
  
  GridPoints <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m) 
  
  const = h * n
  for (i in seq_along(GridPoints))
    y[i] <-  sum(EpanechnikovKernel((GridPoints[i] - x) / h)) / const
  list(GridPoints = GridPoints, y = y)
}


# Implementing AMISE for optimal bandwidth choice

r_hat <- function(x) {
  EmpericalSigma <- sd(x)
  
  EmpericalIQR <- quantile(x, 0.75) - quantile(x, 0.25)
  TheoreticalIQR <- qnorm(0.75) - qnorm(0.25)
  
  SigmaTilde <- min(EmpericalSigma, EmpericalIQR / TheoreticalIQR)
  
  r <- (4 / 3)**(1 / 5) * (SigmaTilde) * length(x)**(-1 / 5)
}


GaussKernel4thOrderDerivative <- function(x) {
  constant <- sqrt(2 * pi)
  temp <- exp((-x**2) / 2) * (x**4 - 6 * x**2 + 3)
  temp / constant
}


EstimatedDensityTilde <- function(x, r) {
  
  result <- 0
  n <- length(x)
  c1 <- n**2 * (sqrt(2) * r)**5
  c2 <- (sqrt(2) * r)
  for(i in seq_along(x)) {
    for(j in seq_along(x)) {
      result <- result + GaussKernel4thOrderDerivative((x[i] - x[j]) / c2)
    }
  }
  
  result / c1
}


EstimatedDensityTildeVectorized <- function(x, r) {
  
  result <- 0
  n <- length(x)
  c1 <- n**2 * (sqrt(2) * r)**5
  c2 <- (sqrt(2) * r)
  for(i in seq_along(x)) {
    result <- result + sum(GaussKernel4thOrderDerivative((x[i] - x) / c2))
  }
  result / c1
}


# Equation 2.3
Optimal_hHat <- function(EstimatedDensity, NumberOfObservations) {
  K_2Norm <- 0.6 # Calculated in class from Epanechnikov kernel
  SigmaK4 <- (1 / 5)**2 # Calculated in class from Epanechnikov kernel
  
  hn_opt <- (K_2Norm / (EstimatedDensity * SigmaK4))**(1 / 5) * NumberOfObservations**(-1 / 5)
  hn_opt
}



# Optimal bandwidth using CV

LeaveOneOutCrossValidation <- function(x, h, kernel = EpanechnikovKernel) {
  
  n_i <- length(x) - 1
  constant <- h * n_i
  EstimaedDensity <- c()
  
  for(i in seq_along(x)) {
    LeaveOutPoint <- x[i]
    RestOfPoints <- x[-i]
    tmp <- 0
    for(j in seq_along(RestOfPoints)) {
      tmp <- tmp + kernel((LeaveOutPoint - RestOfPoints[j]) / h)
    }
    EstimaedDensity[i] <- tmp / constant
  }
  EstimaedDensity
}


LeaveOneOutCrossValidationVectorized <- function(x, h, kernel = EpanechnikovKernel) {
  n_i <- length(x) - 1
  constant <- h * n_i
  EstimaedDensity <- c()
  
  for(i in seq_along(GridPoints)) {
    LeaveOutPoint <- x[i]
    RestOfPoints <- x[-i]
    EstimaedDensity[i] <- sum(kernel((LeaveOutPoint - RestOfPoints) / h)) / constant
  }
  
  EstimaedDensity
}


