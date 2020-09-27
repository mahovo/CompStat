
## # Generate vector of simulated data for single sim path
#x <- 

## Sn
# Sn <- 


## g
#g <- 


## phi
# phi <- 


MCI_ruin <- function(n, m, a, b) {
  # n is number of days the company earns/loses money
  # m is the number withdrawings
  # a and b are the upper and lower bound where the distribution takes place
  numberOfDefaults <- 0
  probability_of_default <- 0
  
  for(i in 1:m) {
    startCapital <- 30
    earnings <- runif(n, a, b)
    for(j in seq_along(earnings)) {
      startCapital <- startCapital + earnings[j]
      if(startCapital <= 0) {
        numberOfDefaults <- numberOfDefaults + 1
        break
      }
    }
  }
  probabilityOfDefault <- numberOfDefaults / m
  probabilityOfDefault
}

