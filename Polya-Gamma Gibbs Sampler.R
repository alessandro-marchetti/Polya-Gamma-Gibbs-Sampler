### Polya-Gamma Gibbs Sampler

library(BayesLogit)

Gibbs_fixed_r <- function(R, burn_in, y, X, B, b, max_r, trunc_lambda) {
  p <- ncol(X)
  n <- nrow(X)
  out <- matrix(0, R, p) # Initialize an empty matrix to store the values
  lam <- matrix(0, R, n)
  
  P <- solve(B) # Prior precision matrix
  Pb <- P %*% b # Term appearing in the Gibbs sampling
  
  # Initialization
  beta <- rep(0, p)
  logr <- log(max_r)
  dif <- (y - max_r) / 2
  
  # Iterative procedure
  for (i in 1:(R + burn_in)) {
    
    # Sampling the Pólya-gamma latent variables
    eta <- c(X %*% beta)
    lambda_trunc <- pmin(exp(eta), trunc_lambda)
    
    omega <- rpg.devroye(num = n, h = y + max_r, z = eta - logr)
    kappa <- omega * logr + dif
    Xk <- crossprod(X, kappa)
    
    # Sampling beta
    eig <- eigen(crossprod(X * sqrt(omega)) + P, symmetric = TRUE)
    
    Sigma <- crossprod(t(eig$vectors) / sqrt(eig$values))
    mu <- Sigma %*% (Xk + Pb)
    
    A1 <- t(eig$vectors) / sqrt(eig$values)
    beta <- mu + c(matrix(rnorm(1 * p), 1, p) %*% A1)
    
    # Store the values after the burn-in period
    if (i > burn_in) {
      out[i - burn_in, ] <- beta
      lam[i - burn_in, ] <- lambda_trunc
    }
  }
  list(out = out, lam = lam)
}
