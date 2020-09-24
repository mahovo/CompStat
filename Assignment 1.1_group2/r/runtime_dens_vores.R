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