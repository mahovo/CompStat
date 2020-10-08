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

defaultTimes <- function(x, startCap = 30)
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
  if (theta == 0) {res <- 0.2564103^100}
  else{
    n <- length(x)
    phi <- (exp(2*theta) - exp(-1.9*theta))/(theta)
    res<- (1 / phi^n) * exp(theta * sum(x)) 
    
    return (res)
  }
}

 # Returns n g-distributed random variables by using inverse transform sampling
gInverseSampling <- function(theta = -0.2, n = 100)
{
  if (theta == 0){res <- runif(n, -1.9, 2)}
  else {
    res <- replicate(n, (log(runif(1)*(exp(2*theta)-exp(-1.9*theta)) + exp(-1.9*theta))) / theta)
  }
  return (res)
}



############################ Monte-Carlo Integration ########################
set.seed(42)
n <- 1e5
tmp <- replicate(n, defaultTimes(runif(100, -1.9, 2)))
mu_hat <- (cumsum(tmp) / (1:n))
mu_hat[n]

## confidence interval
upperConf <- mu_hat + 1.96 * sqrt(mu_hat * (1 - mu_hat) / 1:n)
lowerConf <- mu_hat - 1.96 * sqrt(mu_hat * (1 - mu_hat) / 1:n)

# could have used, emperical var by mu_hat + 1.96 * sd(tmp) / sqrt(1:n)
# but this gives approx same result as using that mu_hat is average on 1 or 0 variables

library(ggplot2)
qplot(1:n, mu_hat) + 
  geom_ribbon(
    mapping = aes(
      ymin = upperConf, 
      ymax = lowerConf
    ), fill = "gray") +
  coord_cartesian(ylim = c(0.001, 0.0025)) +
  geom_line() + 
  geom_point() 


##############################################################################
############################# Importance Sampling ############################
###############################################################################

# Samples from the g-function must be done by inverse transformation
# find the inverse g and input U(0,1)
# next input these samples in the given g-density and the f-density which is 
# constant 1/(2 - (-1.9) ) = 0.2564103 and calc w* = f(x)/g(x)  

# Find a suitable value for theta ie. one that generates some ruins
n <- 1e3
hist(replicate(n, gInverseSampling(theta = -0.2)))
hist(replicate(n, gInverseSampling(theta = 0.2 )))

hist(replicate(n, g(gInverseSampling(), theta = -0.2)))
plot(replicate(n, defaultTimes(gInverseSampling(theta = - 0.2))))
plot(replicate(n, defaultTimes(gInverseSampling(theta = 0.2))))
# Since theta = -0.2 give most values <0 ie. thus most ruins theta = -2 is used



############## The importance sampling with theta = -0.2
fDens <- 0.2564103^100 
n <- 1e5
set.seed(42)
tmp <- vector(mode = "numeric", length = n) #Preallocation
for (i in 1:n) { #Preallocated for-loop should be as fast as replicate
  samples <- gInverseSampling()
  tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = -0.2))
}
mu_hat <- (cumsum(tmp) / (1:n))
mu_hat[n]

############### confidence interval
upperConf <- mu_hat + 1.96 * sd(tmp) / sqrt(1:n)
lowerConf <- mu_hat - 1.96 * sd(tmp) / sqrt(1:n)


library(ggplot2)
qplot(1:n, mu_hat) + 
  geom_ribbon(
    mapping = aes(
      ymin = upperConf, 
      ymax = lowerConf
    ), fill = "gray") +
  coord_cartesian(ylim = c(0.001, 0.0025)) +
  geom_line() + 
  geom_point()


################################ IS with theta = 0
fDens <- 0.2564103^100 

n <- 1e5
set.seed(42)
tmp <- vector(mode = "numeric", length = n) #Preallocation
for (i in 1:n) { #Preallocated for-loop should be as fast as replicate
  samples <- gInverseSampling(theta = 0)
  tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = 0))
}
mu_hat <- (cumsum(tmp) / (1:n))
mu_hat[n]

############### confidence interval
upperConf <- mu_hat + 1.96 * sd(tmp) / sqrt(1:n)
lowerConf <- mu_hat - 1.96 * sd(tmp) / sqrt(1:n)


library(ggplot2)
qplot(1:n, mu_hat) + 
  geom_ribbon(
    mapping = aes(
      ymin = upperConf, 
      ymax = lowerConf
    ), fill = "gray") +
  coord_cartesian(ylim = c(0.001, 0.0025)) +
  geom_line() + 
  geom_point()


################ comparing sd of tmp with different values of theta


thetaValues <- seq(-1, 0.1, length.out =  100)
fDens <- 0.2564103^100 
m <- 1e3
sdValues <-rep(NA, length(thetaValues))

set.seed(42)
for(j in seq_along(thetaValues)) {
  tmp <- rep(NA, m)
  
  for(i in 1:m) {
    samples <- gInverseSampling(theta = thetaValues[j])
    tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = thetaValues[j]))
  }
  print(j)
  sdValues[j] <- sd(tmp)/sqrt(m)
}

library(ggplot2)
thetaVarPlot <- qplot(thetaValues, sdValues) + 
  coord_cartesian(ylim = c(0,0.0006)) +
  geom_point() + geom_vline(xintercept = -0.2, color = 'red') #+ scale_y_continuous(trans = 'log10')

thetaVarPlot
# seems that theta = -0.2 and -0.25 gives lowest sd and therefore best estimation
# but sd is lower for theta <-0.75 and theta > 0, why are those not optimal? See below for the estimated means when using these values


################ comparing emperical sd as used in conf interval with different theta
fDens <- 0.2564103^100 
n <- 1e3


# theta = -0.6
set.seed(42)
tmp <- vector(mode = "numeric", length = n) 
for (i in 1:n) { 
  samples <- gInverseSampling(theta = -0.6)
  tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = -0.6))
}
sdThetaNeg0.6 <- sd(tmp) / sqrt(1:n)

# theta = -0.4
set.seed(42)
tmp <- vector(mode = "numeric", length = n) 
for (i in 1:n) { 
  samples <- gInverseSampling(theta = -0.4)
  tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = -0.4))
}
sdThetaNeg0.4 <- sd(tmp) / sqrt(1:n)

# theta = -0.2
set.seed(42)
tmp <- vector(mode = "numeric", length = n) 
for (i in 1:n) { 
  samples <- gInverseSampling(theta = -0.2)
  tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = -0.2))
}
sdThetaNeg0.2 <- sd(tmp) / sqrt(1:n)

# theta = 0
set.seed(42)
tmp <- vector(mode = "numeric", length = n) 
for (i in 1:n) { 
  samples <- gInverseSampling(theta = 0)
  tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = 0))
}
sdTheta0 <- sd(tmp) / sqrt(1:n)

# theta = 0.2
set.seed(42)
tmp <- vector(mode = "numeric", length = n) 
for (i in 1:n) { 
  samples <- gInverseSampling(theta = 0.2)
  tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = 0.2))
}
sdThetaPos0.2 <- sd(tmp) / sqrt(1:n)


df <- data.frame(n, sdThetaNeg0.6, sdThetaNeg0.4, sdThetaNeg0.2, sdTheta0, sdThetaPos0.2)
ggplot(df, aes(x=n, y = sdThetaNeg0.6)) + 
  geom_line()


library(ggplot2)
colors <- c("theta = -0.6" = "red", "theta = -0.4" = "aquamarine4", "theta = -0.2" = "black", "theta = 0" = "darkblue",
            "theta = 0.2" = 'chocolate3')
thetaSdPlot <- qplot(1:n, sdThetaNeg0.2 ) + 
  coord_cartesian() +
  geom_line() + 
  geom_point() +
  geom_line(aes(y = sdThetaNeg0.4), color = "aquamarine4") +
  geom_point(aes(y = sdThetaNeg0.4), color = "aquamarine4") +
  geom_line(aes(y = sdTheta0), color = 'darkblue') +
  geom_point(aes(y = sdTheta0), color = 'darkblue') +
  geom_line(aes(y = sdThetaPos0.2), color = 'chocolate3') +
  geom_point(aes(y = sdThetaPos0.2), color = 'chocolate3') +
  geom_line(aes(y = sdThetaNeg0.6), color = 'red') +
  geom_point(aes(y = sdThetaNeg0.6), color = 'red') +
  scale_y_continuous(trans = 'log10') + 
  scale_x_continuous(trans = 'log10')
  labs(x = "Paths",
  y = "sd",
  color = "Legend") +
  scale_color_manual(values = c("a", "b"))

#legend not working
thetaSdPlot


############################
##### Saifs version ########
############################
# Standard deviation for every theta

thet <- c(-0.4, -0.2, 0, 0.2, 0.4)
m <- 1e5
stds <- matrix(rep(0, length(thet)), nrow = m, ncol = length(thet) + 1)

for(j in seq_along(thet)) {
  tmp <- rep(NA, m)
  set.seed(42)
  for(i in 1:m) {
    samples <- gInverseSampling(theta = thet[j])
    tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = thet[j]))
  }
  print(j)
  stds[ ,j] <- sd(tmp) / sqrt(1:m)
}
stds[, length(thet) + 1] <- 1:m

library(reshape2)
df <- as.data.frame(stds)
colnames(df) <- c('theta = -0.4', 'theta = -0.2', 'theta = 0', 'theta = 0.2', 'theta = 0.4', 'Paths')

mdf <- melt(df, id.vars = "Paths")
ggplot(mdf,aes(x=Paths, y=value, colour=variable, group=variable)) +
  geom_line() + 
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10')


set.seed(42)
sampl <- gInverseSampling(theta = -0.2)
g(sampl, theta = -0.2)
defaultTimes(sampl)


######################## Compare mean values for using different values of theta #######
thetaValues <- seq(-1, 0.1, length.out =  100)
fDens <- 0.2564103^100 
m <- 1e3
muValues <-rep(NA, length(thetaValues))

set.seed(42)
for(j in seq_along(thetaValues)) {
  tmp <- rep(NA, m)
  
  for(i in 1:m) {
    samples <- gInverseSampling(theta = thetaValues[j])
    tmp[i] <- defaultTimes(samples)*(fDens / g(samples, theta = thetaValues[j]))
  }
  print(j)
  muValues[j] <- mean(tmp)
}

library(ggplot2)
thetaVarPlot <- qplot(thetaValues, muValues) + 
  coord_cartesian(ylim = c(0,0.003)) +
  geom_point() + geom_vline(xintercept = -0.2, color = 'red') +
  geom_hline(yintercept = 0.00169)#+ scale_y_continuous(trans = 'log10')

thetaVarPlot
#black line is result from MC integration


######################### Yet to do
# make rcpp implementation and compare runtime

# compare runtime/accuracy tradeoff for using MC integration and IS
# find optimal theta by using different values and compare runtime
# make IS class taking functions f, g and h as input and having 
# mean, and confInterval as methods
# Camparing runtime of generating samples with different choices of theta



