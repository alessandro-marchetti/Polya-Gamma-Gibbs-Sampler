rm(list=ls())

# The procedure is similar for all the three methods (Polya-Gamma Gibbs sampling, Hamiltonian Monte Carlo, Canale-D'Angelo's bpr):
# Prepare an empty matrix with p columns and rows the number of iterations of the chains times the number of loops of the method

library(BayesLogit)
library(rstan)
library(coda)
library(bpr)

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
max_r <- 100
trunc_lambda <- 1000
R <- 3000

############## n = 25
## Polya-Gamma Gibbs Sampler
samples_gs_p20_n25 <- matrix(NA, nrow = R * length(seeds), ncol = p)

for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  res <- Gibbs_fixed_r(R = R, burn_in = R, y = y_p20_n25, X = X_p20_n25,
                       B = B, b = b, max_r = max_r, trunc_lambda = trunc_lambda)
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_gs_p20_n25[idx, ] <- res$out
  if(i%%5 == 0){
    print(i) 
  }
}
save(samples_gs_p20_n25, file = "gs_p20_n25.RData")

## Hamiltonian Monte Carlo
data_p20_n25 <- list(
  N = 25,
  K = p,
  X = X_p20_n25,
  y = y_p20_n25
)

samples_stan_p20_n25 <- matrix(NA, nrow = 5000 * length(seeds), ncol = p)

for (i in seq_along(seeds)) {
  res <- stan(
    file = "C:/Users/almar/Downloads/poisson_regression.stan",
    data = data_p20_n25,
    iter = 2000,
    chains = 5,
    seed = seeds[i]
  )
  
  # Extract beta samples and store them
  idx <- ((i - 1) * 5000 + 1):(i * 5000)
  samples_stan_p20_n25[idx, ] <- extract(res)$beta
}
save(samples_stan_p20_n25, file = "stan_p20_n25.RData")

## CDMH - bpr
samples_bpr_p20_n25 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p20_n25, X_p20_n25)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr_p20_n25[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}

save(samples_bpr_p20_n25, file = "bpr_p20_n25.RData")


########## n = 50
## Polya-Gamma Gibbs Sampler
samples_gs_p20_n50 <- matrix(NA, nrow = R * length(seeds), ncol = p)

for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  res <- Gibbs_fixed_r(R = R, burn_in = R, y = y_p20_n50, X = X_p20_n50,
                       B = B, b = b, max_r = max_r, trunc_lambda = trunc_lambda)
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_gs_p20_n50[idx, ] <- res$out
  if(i%%5 == 0){
    print(i) 
  }
}
save(samples_gs_p20_n50, file = "gs_p20_n50.RData")

## Hamiltonian Monte Carlo
data_p20_n50 <- list(
  N = 50,
  K = p,
  X = X_p20_n50,
  y = y_p20_n50
)

samples_stan_p20_n50 <- matrix(NA, nrow = 5000 * length(seeds), ncol = p)

for (i in seq_along(seeds)) {
  res <- stan(
    file = "C:/Users/almar/Downloads/poisson_regression.stan",
    data = data_p20_n50,
    iter = 2000,
    chains = 5,
    seed = seeds[i]
  )
  
  # Extract beta samples and store them
  idx <- ((i - 1) * 5000 + 1):(i * 5000)
  samples_stan_p20_n50[idx, ] <- extract(res)$beta
}
save(samples_stan_p20_n50, file = "stan_p20_n50.RData")

## CDMH - bpr
samples_bpr_p20_n50 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p20_n50, X_p20_n50)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr_p20_n50[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}
save(samples_bpr_p20_n50, file = "bpr_p20_n50.RData")


######## n  = 100
## Polya-Gamma Gibbs Sampler
samples_gs_p20_n100 <- matrix(NA, nrow = R * length(seeds), ncol = p)

for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  res <- Gibbs_fixed_r(R = R, burn_in = R, y = y_p20_n100, X = X_p20_n100,
                       B = B, b = b, max_r = max_r, trunc_lambda = trunc_lambda)
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_gs_p20_n100[idx, ] <- res$out
  if(i%%5 == 0){
    print(i) 
  }
}
save(samples_gs_p20_n100, file = "gs_p20_n100.RData")

## Hamiltonian Monte Carlo
data_p20_n100 <- list(
  N = 100,
  K = p,
  X = X_p20_n100,
  y = y_p20_n100
)

samples_stan_p20_n100 <- matrix(NA, nrow = 5000 * length(seeds), ncol = p)

for (i in seq_along(seeds)) {
  res <- stan(
    file = "C:/Users/almar/Downloads/poisson_regression.stan",
    data = data_p20_n100,
    iter = 2000,
    chains = 5,
    seed = seeds[i]
  )
  
  # Extract beta samples and store them
  idx <- ((i - 1) * 5000 + 1):(i * 5000)
  samples_stan_p20_n100[idx, ] <- extract(res)$beta
}
save(samples_stan_p20_n100, file = "stan_p20_n100.RData")

## CDMH - bpr
samples_bpr_p20_n100 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p20_n100, X_p20_n100)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr_p20_n100[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}
save(samples_bpr_p20_n100, file = "bpr_p20_n100.RData")


########## n = 200
## Polya-Gamma Gibbs Sampler
samples_gs_p20_n200 <- matrix(NA, nrow = R * length(seeds), ncol = p)

for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  res <- Gibbs_fixed_r(R = R, burn_in = R, y = y_p20_n200, X = X_p20_n200,
                       B = B, b = b, max_r = max_r, trunc_lambda = trunc_lambda)
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_gs_p20_n200[idx, ] <- res$out
  if(i%%5 == 0){
    print(i) 
  }
}
save(samples_gs_p20_n200, file = "gs_p20_n200.RData")

## Hamiltonian Monte Carlo
data_p20_n200 <- list(
  N = 200,
  K = p,
  X = X_p20_n200,
  y = y_p20_n200
)

samples_stan_p20_n200 <- matrix(NA, nrow = 5000 * length(seeds), ncol = p)

for (i in seq_along(seeds)) {
  res <- stan(
    file = "C:/Users/almar/Downloads/poisson_regression.stan",
    data = data_p20_n200,
    iter = 2000,
    chains = 5,
    seed = seeds[i]
  )
  
  # Extract beta samples and store them
  idx <- ((i - 1) * 5000 + 1):(i * 5000)
  samples_stan_p20_n200[idx, ] <- extract(res)$beta
}
save(samples_stan_p20_n200, file = "stan_p20_n200.RData")

## CDMH - bpr
samples_bpr_p20_n200 <- matrix(NA, nrow = R * length(seeds), ncol = p)
data_train = data.frame("y" = y_p20_n200, X_p20_n200)

for (i in seq_along(seeds)) {
  res <- sample_bpr(y ~ . -1, data = data_train, iter = 2*R, burnin = R,
                    prior = list(type = "gaussian", b = b, B = B),
                    pars = list(max_dist = 0.5, max_dist_burnin = 0.5))
  
  # Extract beta samples and store them
  idx <- ((i - 1) * R + 1):(i * R)
  samples_bpr_p20_n200[idx, ] <- res$sim$beta[3001:6000,]
  if(i %% 5 == 0){
    print(i)
  }
}
save(samples_bpr_p20_n200, file = "bpr_p20_n200.RData")

