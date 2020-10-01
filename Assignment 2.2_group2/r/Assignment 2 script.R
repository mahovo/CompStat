## Monte Carlo functions ====

x_mat = function(n, m, rfunc = runif(n*m, -1.9, 2.0)){
  matrix(rfunc, n, m)
}

defaultTimesMatrix <- function(x_mat, startCap)
{
  default <- rep(0, ncol(x_mat))
  
  for (i in 1:ncol(x_mat)) {
    path <- startCap + cumsum(x_mat[ , i])
    if (sum(as.numeric(path<0)) >0) default[i] <- 1
  }
  
  return(default)
}

defaultTimes <- function(x, startCap = 30, n = 100)
{
  default <- 0
  path <- startCap + cumsum(x)
  if (any(path <= 0)) default <- 1
  return(default)
}


gMat <- function(X_mat, theta = 0.4)
{
  phi <- (exp(-1.9*theta)*(exp(3.9*theta))-1)/(theta)
  n <- 100
  res <- rep(0, 100)
  
  for (i in 1:n) {
    res[i]<- (1 / phi) * exp(theta * sum(X_mat[,i])) 
  }
  return (res)
}

f <- function(x)
{
  return(30 + cumsum(x))
}

g <- function(x, theta = -0.2)
{
  n <- length(x)
  phi <- (exp(2*theta) - exp(-1.9*theta))/(theta)
  res<- (1 / phi^n) * exp(theta * sum(x)) 

  return (res)
}

 # Returns n g-distributed random variables by using inverse transform sampling
gInverseSampling <- function(theta = -0.2, n = 100)
{
  res <- replicate(n, (log(runif(1)*(exp(2*theta)-exp(-1.9*theta)) + exp(-1.9*theta))) / theta)
  #u <- runif(1)
  #res <- (log(u*(exp(2*theta)-exp(-1.9*theta)) + exp(-1.9*theta))) / theta
  return (res)
}

importanceSampler <- function(sampler, evalFunc, targetDens, sampledDens)
{
  sample <- sampler
  weight <- targetDens / sampledDens
  return (evalFunc*weight)
}



############################ Monte-Carlo Integration ########################
set.seed(42)
n <- 1e5
tmp <- replicate(n, defaultTimes(runif(100, -1.9, 2)))
mu_hat <- (cumsum(tmp) / (1:n))
mu_hat[n]

## confidence interval
upperConf <- mu_hat + 1.96 * sqrt(mu_hat * (1 - mu_hat) / n)
lowerConf <- mu_hat - 1.96 * sqrt(mu_hat * (1 - mu_hat) / n)


library(ggplot2)
qplot(1:n, mu_hat) + 
  geom_ribbon(
    mapping = aes(
      ymin = upperConf, 
      ymax = lowerConf
    ), fill = "gray") +
  coord_cartesian(ylim = c(0, 0.004)) +
  geom_line() + 
  geom_point() 



############################# Importance Sampling ############################


# Samples from the g-function must be done by inverse transformation
# find the inverse g and input U(0,1)
# next input these samples in the given g-density and the f-density which is 
# constant 1/(2 - (-1.9) ) = 0.2564103 and calc w* = f(x)/g(x)  

# Find a suitable value for theta ie. one that generates some ruins
n <- 1e5
hist(replicate(n, gInverseSampling(theta = -0.2)))
hist(replicate(n, gInverseSampling(theta = 0.2 )))

hist(replicate(n, g(gInverseSampling(), theta = -0.2)))
# Since theta = -0.2 give most values <0 ie. thus most ruins theta = -2 is used

#The importance sampling
fDens <- 0.2564103^100 

n <- 1e5
tmp <- replicate(n, defaultTimes(gInverseSampling())*(fDens / g(gInverseSampling(), theta = -0.2)))
mu_hat <- (cumsum(tmp) / (1:n))
mu_hat[n]

## confidence interval
upperConf <- mu_hat + 1.96 * sqrt(mu_hat * (1 - mu_hat) / n)
lowerConf <- mu_hat - 1.96 * sqrt(mu_hat * (1 - mu_hat) / n)


library(ggplot2)
qplot(1:n, mu_hat) + 
  geom_ribbon(
    mapping = aes(
      ymin = upperConf, 
      ymax = lowerConf
    ), fill = "gray") +
  coord_cartesian(ylim = c(0, 0.4)) +
  geom_line() + 
  geom_point()




### to do
# make rcpp implementation and compare runtime

# see that theta = 0 makes g = f, and make the gInverse defined for theta = 0
# compare runtime/accuracy tradeoff for using MC integration and IS
# find optimal theta by using different values and compare runtime





