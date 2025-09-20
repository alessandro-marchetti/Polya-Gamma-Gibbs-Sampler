rm(list = ls())
library(bpr)

p = 5
set.seed(9)
sim_mix <- function(n, p, m, v)
{
  k <- sample(1:length(m), n, replace = TRUE, prob = p)
  return(rnorm(n, m[k], v[k]))
}
beta_p5 <- sim_mix(p, c(0.3,0.5,0.4,0.4), c(0.6,-0.1,2.5,-2.4), c(0.7,0.7,0.7,0.7))


seeds = 1:50

for(seed in seeds){
  n = 1000
  set.seed(seed)
  X <- matrix(rep(NA, n*p), n, p)
  X[,1] <- 1
  X[,2] <- runif(n, -2, 1)
  X[,3] <- rgamma(n, 2,2)
  X[,4] <- rnorm(n, 2,2)
  X[,5] = rexp(n, 2) + rpois(n, 1)
  
  for(j in c(2:5))
  {
    X[,j] = (X[,j] - mean(X[,j]))/sd(X[,j])
  }
  
  ix = ( exp(X%*%beta_p5) > 1 ) & ( exp(X%*%beta_p5) < 200 )
  if(sum(ix) < 200) stop()
  X = X[ix,]
  
  #-------------#   n=25   #-------------#
  n = 25
  ixn = sample.int(sum(ix),n)
  X_p5_n25 <- X[ixn,]
  
  set.seed(seed)
  y_p5_n25 <- rpois(n, exp(X_p5_n25%*%beta_p5))
  
  #-------------#   n=50   #-------------#
  n = 50
  ixn = sample.int(sum(ix),n)
  X_p5_n50 <- X[ixn,]
  
  set.seed(seed)
  y_p5_n50 <- rpois(n, exp(X_p5_n50%*%beta_p5))
  
  #-------------#   n=100   #-------------#
  n = 100
  ixn = sample.int(sum(ix),n)
  X_p5_n100 <- X[ixn,]
  
  set.seed(seed)
  y_p5_n100 <- rpois(n, exp(X_p5_n100%*%beta_p5))
  
  #-------------#   n=200   #-------------#
  n = 200
  ixn = sample.int(sum(ix),n)
  X_p5_n200 <- X[ixn,]
  
  set.seed(seed)
  y_p5_n200 <- rpois(n, exp(X_p5_n200%*%beta_p5))
  
}

b <- rep(0,p)
B <- diag(10,p)
R <- 3000

samples_bpr0_p5_n25 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p5_n25, X_p5_n25)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p5_n25[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p5_n25, file = "bpr0_p5_n25.RData")

samples_bpr0_p5_n50 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p5_n50, X_p5_n50)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p5_n50[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p5_n50, file = "bpr0_p5_n50.RData")


samples_bpr0_p5_n100 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p5_n100, X_p5_n100)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p5_n100[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p5_n100, file = "bpr0_p5_n100.RData")


samples_bpr0_p5_n200 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p5_n200, X_p5_n200)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p5_n200[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p5_n200, file = "bpr0_p5_n200.RData")

############ p = 10

rm(list = ls())
p = 10
set.seed(9)
sim_mix <- function(n, p, m, v)
{
  k <- sample(1:length(m), n, replace = TRUE, prob = p)
  return(rnorm(n, m[k], v[k]))
}
beta_p10 <- sim_mix(p, c(0.3,0.5,0.4,0.4), c(0.6,-0.1,2.5,-2.4), c(0.7,0.7,0.7,0.7))



seeds = 1:50

for(seed in seeds){
  n = 1000
  set.seed(seed)
  X <- matrix(rep(NA, n*p), n, p)
  for(j in 1:2)
  {
    X[,j] <- rnorm(n, 0, 1)
    X[,j+2] <- runif(n, -2, 1)
    X[,j+4] <- rgamma(n, 2,2)
    X[,j+6] <- rnorm(n, 2,2)
  }
  X[, 7] = rbinom(n, 1, 0.4)
  X[, 8] = rbinom(n, 1, 0.6)
  X[, 9] = rnorm(n, 5, 6)
  X[,10] = rexp(n, 2) + rpois(n, 1)
  X[,1] = 1
  
  for(j in c(2:6,9:10))
  {
    X[,j] = (X[,j] - mean(X[,j]))/sd(X[,j])
  }
  
  ix = ( exp(X%*%beta_p10) > 1 ) & ( exp(X%*%beta_p10) < 200 )
  if(sum(ix) < 200) stop()
  X = X[ix,]
  
  #-------------#   n=25   #-------------#
  n = 25
  ixn = sample.int(sum(ix),n)
  X_p10_n25 <- X[ixn,]
  
  set.seed(seed)
  y_p10_n25 <- rpois(n, exp(X_p10_n25%*%beta_p10))
  
  #-------------#   n=50   #-------------#
  n = 50
  ixn = sample.int(sum(ix),n)
  X_p10_n50 <- X[ixn,]
  
  set.seed(seed)
  y_p10_n50 <- rpois(n, exp(X_p10_n50%*%beta_p10))
  
  #-------------#   n=100   #-------------#
  n = 100
  ixn = sample.int(sum(ix),n)
  X_p10_n100 <- X[ixn,]
  
  set.seed(seed)
  y_p10_n100 <- rpois(n, exp(X_p10_n100%*%beta_p10))
  
  #-------------#   n=200   #-------------#
  n = 200
  ixn = sample.int(sum(ix),n)
  X_p10_n200 <- X[ixn,]
  
  set.seed(seed)
  y_p10_n200 <- rpois(n, exp(X_p10_n200%*%beta_p10))
}



b <- rep(0,p)
B <- diag(10,p)
R <- 3000

samples_bpr0_p10_n25 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p10_n25, X_p10_n25)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p10_n25[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p10_n25, file = "bpr0_p10_n25.RData")

samples_bpr0_p10_n50 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p10_n50, X_p10_n50)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p10_n50[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p10_n50, file = "bpr0_p10_n50.RData")


samples_bpr0_p10_n100 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p10_n100, X_p10_n100)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p10_n100[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p10_n100, file = "bpr0_p10_n100.RData")


samples_bpr0_p10_n200 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p10_n200, X_p10_n200)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p10_n200[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p10_n200, file = "bpr0_p10_n200.RData")



########## p = 20
rm(list = ls())


p = 20
set.seed(9)
sim_mix <- function(n, p, m, v)
{
  k <- sample(1:length(m), n, replace = TRUE, prob = p)
  return(rnorm(n, m[k], v[k]))
}
beta_p20 <- sim_mix(p, c(0.3,0.5,0.4,0.4), c(0.6,-0.1,2.5,-2.4), c(0.7,0.7,0.7,0.7))



seeds = 1:50

for(seed in seeds){
  n = 1000
  set.seed(seed)
  X <- matrix(rep(NA, n*p), n, p)
  for(j in 1:4)
  {
    X[,j] <- rnorm(n, 0, 1)
    X[,j+4] <- runif(n, -2, 1)
    X[,j+8] <- rgamma(n, 2,2)
    X[,j+12] <- rnorm(n, 2,2)
  }
  X[, 17] = rbinom(n, 1, 0.4)
  X[, 18] = rbinom(n, 1, 0.7)
  X[,19:20] = rexp(2*n, 2) + rpois(n*2, 1)
  X[,1] = 1
  head(X)
  
  for(j in c(2:16,19:20))
  {
    X[,j] = (X[,j] - mean(X[,j]))/sd(X[,j])
  }
  
  ix = ( exp(X%*%beta_p20) > 1 ) & ( exp(X%*%beta_p20) < 200 )
  if(sum(ix) < 200) stop()
  X = X[ix,]
  
  #-------------#   n=25   #-------------#
  n = 25
  ixn = sample.int(sum(ix),n)
  X_p20_n25 <- X[ixn,]
  
  set.seed(seed)
  y_p20_n25 <- rpois(n, exp(X_p20_n25%*%beta_p20))
  
  
  #-------------#   n=50   #-------------#
  n = 50
  ixn = sample.int(sum(ix),n)
  X_p20_n50 <- X[ixn,]
  
  set.seed(seed)
  y_p20_n50 <- rpois(n, exp(X_p20_n50%*%beta_p20))
  
  
  #-------------#   n=100   #-------------#
  n = 100
  ixn = sample.int(sum(ix),n)
  X_p20_n100 <- X[ixn,]
  
  set.seed(seed)
  y_p20_n100 <- rpois(n, exp(X_p20_n100%*%beta_p20))
  
  #-------------#   n=200   #-------------#
  n = 200
  ixn = sample.int(sum(ix),n)
  X_p20_n200 <- X[ixn,]
  
  set.seed(seed)
  y_p20_n200 <- rpois(n, exp(X_p20_n200%*%beta_p20))
  
  
}


b <- rep(0,p)
B <- diag(10,p)
R <- 3000


samples_bpr0_p20_n25 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p20_n25, X_p20_n25)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p20_n25[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p20_n25, file = "bpr0_p20_n25.RData")

samples_bpr0_p20_n50 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p20_n50, X_p20_n50)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p20_n50[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p20_n50, file = "bpr0_p20_n50.RData")


samples_bpr0_p20_n100 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p20_n100, X_p20_n100)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p20_n100[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p20_n100, file = "bpr0_p20_n100.RData")


samples_bpr0_p20_n200 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p20_n200, X_p20_n200)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5),
                    state = rep(0,p))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr0_p20_n200[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr0_p20_n200, file = "bpr0_p20_n200.RData")

