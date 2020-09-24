

# Remove all
rm(list=ls())


# Get libraries
library(ggplot2)
library(microbenchmark)
library(kedd)

# Get functions from source file
source("~/Desktop/CompStat/Assigntment 1.2/Functions assignment 1 - density estimation.R")

# Setwidth
setwd("~/Desktop/CompStat/Exercise A.5")

# Get data
infrared_data <- read.table("infrared.dat", header=TRUE)
logF12 <- log(infrared_data$F12)


# Making dataframe of the log F12 data using epanechinkov kernel
# and using R's density function
RDensFunc <- as.data.frame(density(logF12, bw = 0.5, 
                              kernel = "epanechnikov")[1:2])

# plot the density
p0 <- ggplot(data = infrared_data, 
             aes(x = log(F12))) + geom_histogram(aes(y = ..density..), 
                                                 alpha = 0.5)

# F12 data with smoothed density plot
p1 <- p0 + geom_line(data = RDensFunc, aes(x = x, y = y), color = 'blue')

# Visulized plot
p1

# Make plot of sorted F12 data
plot(sort(EpanechnikovKernel(logF12)))

# Comparing our density implementation with R's density function
EpanechnikovDensity <- as.data.frame(VectorizedKernelDensity(logF12, 0.5))
p1 + geom_line(data = EpanechnikovDensity, aes(x = GridPoints, y = y), color = 'red')

# Using AMISE to find optimal bandwidth
r_hats <- r_hat(logF12)
fTilde <- EstimatedDensityTilde(logF12, r_hats)
fTildeVec <- EstimatedDensityTildeVectorized(logF12, r_hats)
h_opt <- Optimal_hHat(fTilde, length(logF12))
h_opt_vec <- Optimal_hHat(fTildeVec, length(logF12))

h_opt
h_opt_vec


DensityWithOptimalh <- as.data.frame(VectorizedKernelDensity(logF12, h = h_opt))
p1 + geom_line(data = DensityWithOptimalh, aes(x = GridPoints, y = y), color = 'purple')




# Sequences of bandwidths
h <- seq(0.1, 5, length.out = 50)
CV <- c()
for(i in seq_along(h)){
  CV[i] <- sum(log(LeaveOneOutCrossValidationVectorized(x = logF12, h = h[i])))
}

# Finding opt h
h_opts <- h[which.max(CV)]
qplot(h, 
      CV, 
      main = 'Optimal bandwith',
      ylab = 'MLCV') + geom_line() + geom_vline(xintercept = h_opts,
                                                color = 'red')



