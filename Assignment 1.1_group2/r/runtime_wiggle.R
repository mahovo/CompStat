
conf <- expand.grid(
  fun = c("wiggle_gauss", "wiggle_gauss_vec", "wiggle_gauss_outer"),
  n = 2^(5:11),
  r = r_hat(x, n)
)
set.seed(1234)
x <- rnorm(2^11)
calls <- paste0(conf[, 1], "(x[1:", conf[, 2], "], r = ", conf[, 3], ")")
expr_list <- lapply(calls, function(x) parse(text = x)[[1]])
kern_benchmarks <- microbenchmark(list = expr_list, times = 40L)


class(kern_benchmarks) <- "data.frame"
kern_benchmarks <- dplyr::bind_cols(conf, expr = calls) %>% 
  dplyr::left_join(kern_benchmarks, .)
# kern_benchmarks$m <- factor(
#   kern_benchmarks$m, 
#   levels = c(32, 128, 512, 2048),
#   labels = c("m = 32", "m = 128", "m = 512", "m = 2048")
# )
ggplot(kern_benchmarks, aes(x = n, y = time, color = fun)) + 
  geom_abline(intercept = 15, slope = 1, color = "gray", linetype = 2) +
  stat_summary(fun = "median", geom = "line") + 
  stat_summary(fun = "median", geom = "point") + 
  # facet_wrap(~ m) + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous("Time (ms)", trans = "log2", 
                     breaks = c(1e5, 1e6, 1e7, 1e8), 
                     labels = c("0.1", "1", "10", "100")) +
  scale_color_discrete("Function:") + 
  theme(legend.position="top")
