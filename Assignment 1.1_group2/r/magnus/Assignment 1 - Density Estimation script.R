setwd("C:/Users/PC/OneDrive - CBS - Copenhagen Business School/CBS/5. år/Comp Stat/Exercises/data")
setwd("C:/Users/Magnus/OneDrive - CBS - Copenhagen Business School/CBS/5. år/Comp Stat/Exercises/data")
infrared <- read.table("infrared.dat", header = TRUE)
F12 <- infrared$F12
logF12 <- log(F12)

source("C:/Users/PC/OneDrive - CBS - Copenhagen Business School/CBS/5. år/Comp Stat/Assignments/Assignment 1/Assignment 1 - Density Estimation functions.R")


#################################
## Calculations            ######
## for the epanechnikov
## kernel                 ######
#################################

## Calculate the variance 
#The variance is calculated by
# Var = int^[1]_[-1] x^2 K(x) dx 
# since the kernel is symmetric this can be rewritten as
# Var = 2*int^[1]_[0] x^2 K(x) dx 
# 2*int^{1}_{0} x^2*3/4(1-x^2) dx 
# 8/4*3/4 * [(x^3/3)-(x^5/5)]^[1]_[0]
# 24/16 * (1/3 - 1/5)
# = 1/5

# Calculate L2-norm(K)^2
# in general
# LP-norm(f(x)) = (Int^{inf}_{-inf} abs(f(x))^{P} dx )^{1/P}
# in this case we have 
# L2-norm(K) = (Int^{1}_{-1} abs(K)^2 dx)^{ (1/2) * 2}
# Int^{1}_{-1}abs(K)^2 dx
# Int^{1}_{-1} abs(3/4*(1-x^2))^2 dx
# since K is always positive and especially K^2 is
#then we can remove the absolute operator
# 3/4^2 * Int^{1}_{-1} (1-x^2))^2 dx
# 9/16 * Int^{1}_{-1} 1 + x^4 - 2x^2 dx
# 9/16  * (Int^{1}_{-1} 1 dx + Int^{1}_{-1} x^4 dx
#  - Int^{1}_{-1} 2x^2 dx )
# 9/16 * (2 + (1/5 - (-1)/5) - (2/3 - (-2)/3))
# 9/16 * (2 + 2/5 - 4/3)
# = 3/5 = 0.6


################################################
#### Calculate optimal bandwidth          ######
################################################

########## AMISE plug-in method

my_s_tilde <- sigma_tilde(sd(logF12), IQR(logF12))
r_hat <- silverman(my_s_tilde, length(logF12))

#calculate f'' using gauss pilot with r_hat
fmm <- fdblmark_pilot(logF12, r_hat)

#plug f'' into the oracle bandwidth formula and cumpute h
#see ###calculations### above for some of these values
K <- 0.6
sigma_K <- 1/5
bw <- oracle_bw(K, fmm, sigma_K, length(logF12))
bw


################################################
##### Make density from Epanecknikov kernel ####
################################################

dens <- as.data.frame(density(logF12, bw = bw,
                              kernel = "epanechnikov")[1:2])
#my_dens <- data.frame('x' = logF12,'y' =  ep_kernel_dens(logF12, 5))
my_dens <- as.data.frame(ep_kernel_dens(logF12, bw * sqrt(5)))

#plotting densities
library('ggplot2')
p0 <- ggplot(data = infrared, aes(x = logF12)) + 
  geom_histogram(aes( y = ..density..) , alpha = 0.4,
                 binwidth = 0.5)

p0
p1 <- p0 + geom_line(data = dens, aes(x = x, y = y),
                     color = "red", bw = bw)
p1
p2 <- p1 + geom_line(data = my_dens, aes(x = x, y = y),
                     color = "blue")
p2


#####################################################
########## Remarks for the comparison of density ####
########## and implemented epenechnikov kernel   ####
#####################################################

#### Implemented epKernel will not be similar to 
#R's version because R normalize all kernels to have
#variance 1, which can be shown is not the case with 
#the formula for the ep-kernel. This normalization
#is effectively a reparametrazasation of the bandwidth
#to get the same bw you'll have to multiply with 
# sqrt(5) - this number apperently comes from some 
#math.

#remember to multiply the bandwidth with sqrt(5) 
#when comparing performance with R's density, so that
#the compared functions actaully return the same reseult

#density is much faster because it uses binning.
#binning is first seeing how many times the same kernel is
#calculated and then only calculating this once.
#so it only calculates when there is a new set of points
#in the kernel. This increases performance in trade for a 
#small accuracy error.




###################
### Profiling #####
###################
library("profvis")
sim_data <- rnorm(1000)

#bw selector
profvis({
  my_s_tilde <- sigma_tilde(sd(sim_data), IQR(sim_data))
  r_hat <- silverman(my_s_tilde, length(sim_data))
  fmm <- fdblmark_pilot(sim_data, r_hat)
  K <- 0.6
  sigma_K <- 1/5
  bw <- oracle_bw(K, fmm, sigma_K, length(sim_data)) #bottleneck
})



#kern dens
profvis(ep_kernel_dens(sim_data, bw))  #bottleneck in ifelse and for-loop



###################
### Benchmark #####
###################

########### Test of similarity ##################
sim_data <- rnorm(2^11)

my_s_tilde <- sigma_tilde(sd(sim_data), IQR(sim_data))
r_hat <- silverman(my_s_tilde, length(sim_data))

fmm <- fdblmark_pilot(sim_data, r_hat)
K <- 0.6
sigma_K <- 1/5
binwidth <- oracle_bw(K, fmm, sigma_K, length(sim_data))

dens <- as.data.frame(density(sim_data, bw = binwidth,
                              kernel = "epanechnikov")[1:2])
my_dens <- as.data.frame(ep_kernel_dens(sim_data, binwidth * sqrt(5)))

# plot of densities
library('ggplot2')
sim_data <- as.data.frame(sim_data)
p0 <- ggplot(data = sim_data, aes(x = sim_data)) + 
  geom_histogram(aes( y = ..density..) , alpha = 0.4)

p0
p1 <- p0 + geom_line(data = dens, aes(x = x, y = y),
                     color = "red")
p1
p2 <- p1 + geom_line(data = my_dens, aes(x = x, y = y),
                     color = "blue")
p2 ###Seems to give exact same density


######## plot of differences
#dens_diff <- as.data.frame(my_dens)
#dens_diff$y <- my_dens$y - dens$y
#This does not work because the density is not 
#evaluated for the same x-values


##### benchmarking #######
sim_data <- rnorm(2^11)
# bw selector
K <- 0.6
sigma_K <- 1/5
bw_bench <- microbenchmark({
  my_s_tilde <- sigma_tilde(sd(sim_data), IQR(sim_data))
  r_hat <- silverman(my_s_tilde, length(sim_data))
  fmm <- fdblmark_pilot(sim_data, r_hat)
  bw <- oracle_bw(K, fmm, sigma_K, length(sim_data)) #bottleneck
})

summary(bw_bench)
# density estimation


library('microbenchmark')

#for n=2048
kern_bench <- microbenchmark(
  ep_kernel_dens(sim_data, binwidth * sqrt(5)),
  density(sim_data, bw = binwidth,
          kernel = "epanechnikov")
)

sum_kern_bench2048 <- summary(kern_bench)

# n = 1024
kern_bench <- microbenchmark(
  ep_kernel_dens(sim_data[1:1024], binwidth * sqrt(5)),
  density(sim_data[1:1024], bw = binwidth,
          kernel = "epanechnikov")
)

sum_kern_bench1024 <- summary(kern_bench)

# n = 512
kern_bench <- microbenchmark(
  ep_kernel_dens(sim_data[1:512], binwidth * sqrt(5)),
  density(sim_data[1:512], bw = binwidth,
          kernel = "epanechnikov")
)

sum_kern_bench512 <- summary(kern_bench)

# n = 128
kern_bench <- microbenchmark(
  ep_kernel_dens(sim_data[1:128], binwidth * sqrt(5)),
  density(sim_data[1:128], bw = binwidth,
          kernel = "epanechnikov")
)

sum_kern_bench128 <- summary(kern_bench)

# n = 256
kern_bench <- microbenchmark(
  ep_kernel_dens(sim_data[1:256], binwidth * sqrt(5)),
  density(sim_data[1:256], bw = binwidth,
          kernel = "epanechnikov")
)

sum_kern_bench256 <- summary(kern_bench)


# n = 64
kern_bench <- microbenchmark(
  ep_kernel_dens(sim_data[1:64], binwidth * sqrt(5)),
  density(sim_data[1:64], bw = binwidth,
          kernel = "epanechnikov")
)

sum_kern_bench64 <- summary(kern_bench)

# n = 32
kern_bench <- microbenchmark(
  ep_kernel_dens(sim_data[1:32], binwidth * sqrt(5)),
  density(sim_data[1:32], bw = binwidth,
          kernel = "epanechnikov")
)

sum_kern_bench32 <- summary(kern_bench)

my_dens_time <- c(sum_kern_bench32$median[1], sum_kern_bench64$median[1], sum_kern_bench128$median[1], 
  sum_kern_bench256$median[1], 
  sum_kern_bench512$median[1], sum_kern_bench1024$median[1], 
  sum_kern_bench2048$median[1])

r_dens_time <- c(sum_kern_bench32$median[2], sum_kern_bench64$median[2], sum_kern_bench128$median[2], 
            sum_kern_bench256$median[2], 
            sum_kern_bench512$median[2], sum_kern_bench1024$median[2], 
            sum_kern_bench2048$median[2])


## runtime plot for n= 32   64  128  256  512 1024 2048
n = 2^(5:11)
kern_benchmarks <- data.frame( my_dens_time, r_dens_time, n)
#r_dens_time <- data.frame(r_dens_time, n)

runtime_plot <- ggplot(data = kern_benchmarks, aes(x = n, y = my_dens_time)) + 
  geom_abline(intercept = 15, slope = 1, color = "gray", linetype = 2) + 
  geom_line(data = kern_benchmarks) + 
  stat_summary(data = kern_benchmarks) + scale_x_continuous(trans = "log2") + 
  scale_y_continuous("Time (ms)", trans = "log2", 
                     breaks = c(1e5, 1e6, 1e7, 1e8), 
                     labels = c("0.1", "1", "10", "100")) +
  scale_color_discrete("Function:") + 
  theme(legend.position="top")
  
runtime_plot
runtime_plot2 <- runtime_plot + geom_line(r_dens_time, aes(x = n, y = r_dens_time))


## autoplot
autoplot(kern_bench) + 
  geom_jitter(position = position_jitter(0.2, 0), aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  theme(legend.position = "none")



#systematic benchmarking
library('dplyr')
conf <- expand.grid(
  fun = c("ep_kernel_dens", "density"),
  n = 2^(5:11), 
  m = 2^c(5, 7, 9, 11)
)
set.seed(1234)
x <- rnorm(2^11)
calls <- paste0(conf[, 1], "(x[1:", conf[, 2], "], h = 0.2, m = ", conf[, 3], ")")
expr_list <- lapply(calls, function(x) parse(text = x)[[1]])
kern_benchmarks <- microbenchmark(list = expr_list, times = 40L)


class(kern_benchmarks) <- "data.frame"
kern_benchmarks <- dplyr::bind_cols(conf, expr = calls) %>% 
  dplyr::left_join(kern_benchmarks, .)
kern_benchmarks$m <- factor(
  kern_benchmarks$m, 
  levels = c(32, 128, 512, 2048),
  labels = c("m = 32", "m = 128", "m = 512", "m = 2048")
)
ggplot(kern_benchmarks, aes(x = n, y = time, color = fun)) + 
  geom_abline(intercept = 15, slope = 1, color = "gray", linetype = 2) +
  stat_summary(fun = "median", geom = "line") + 
  stat_summary(fun = "median", geom = "point") + 
  facet_wrap(~ m) + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous("Time (ms)", trans = "log2", 
                     breaks = c(1e5, 1e6, 1e7, 1e8), 
                     labels = c("0.1", "1", "10", "100")) +
  scale_color_discrete("Function:") + 
  theme(legend.position="top")



