source("Combinatorics.R")
library(glmnet)
# Functions used while updating the estimates of SHIM 

# Criterion to stop the model
compute_rel_dif <- function(Q_old, Q_new)
{
  rel_dif <- abs(Q_old - Q_new) / abs(Q_old)
  return(rel_dif)
}


# Kappa functions
kappa0 <- function(x) {
  k0 <- rep(0, length(x))
  tt <- which((x <= 700) & (x != 0))
  k0[tt] <- log((exp(x[tt]) - 1) / x[tt])
  tt <- which(x > 700)
  k0[tt] <- x[tt] - log(x[tt])
  return(k0)
}


kappa1 <- function(x) {
  k1 <- rep(1 / 2, length(x))
  tt <- which(abs(x) <= 0.0001)
  k1[tt] <- 1 / 2 + x[tt] / 12 - x[tt] ^ 3 / 720 + x[tt] ^ 5 / 30240
  tt <- which(abs(x) > 0.0001)
  k1[tt] <- 1 - 1 / (x[tt]) - 1 / (1 - exp(x[tt]))
  return(k1)
}


kappa2 <- function(x) {
  k2 <- rep(1 / 12, length(x))
  tt <- which(abs(x) <= 0.015)
  k2[tt] <- 1 / 12 - x[tt] ^ 2 / 240 + x[tt] ^ 4 / 6048
  tt <- which(abs(x) > 0.015)
  k2[tt] <- 1 / (x[tt]) ^ 2 + 1 / (2 - 2 * cosh(x[tt]))
  return(k2)
}


#approximation of link function g
g.link <- function(x) {
  tt <- apply(
    as.matrix(x),
    1,
    FUN = function(v)
      min(max(v, 0.001), 0.999)
  )
  g <- 3.5 * tan(pi * (2 * tt - 1) / 2)
  return(g)
}


#Lasso for continuous Bernoulli distribution
irlasso.cb <- function(X,
                       Y,
                       lambda,
                       w.lambda = NULL,
                       beta0 = NULL,
                       centering = TRUE,
                       scaling = TRUE,
                       intercept = TRUE,
                       maxit = 10,
                       tol = 0.0545,
                       sd.tol = 1e-6,
                       verbose = FALSE,
                       C = 0) {
  if (verbose)
    print("Performing IRLASSO-NEW")
  
  
  # Get variables
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Center variables (column-wise)
  mu.X <- rep(0, ncol(X))
  mu.Y <- rep(0, ncol(Y))
  if (centering) {
    mu.X <- apply(X, 2, mean)
    X <- apply(X, 2, function(v) {
      v - mean(v)
    })
  }
  
  # Scale variables (column-wise)
  sd.X <- rep(1, ncol(X))
  sd.Y <- rep(1, ncol(Y))
  if (scaling) {
    sd.X <- apply(X, 2, function(v) {
      max(sd(v), sd.tol)
    })
    X <- apply(X, 2, function(v) {
      v / max(sd(v), sd.tol)
    })
  }
  
  # Intercept
  if (intercept) {
    X <- cbind(1, X)
  }
  
  # Get parameters
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # Initialization
  if (is.null(beta0)) {
    beta0 <- rep(0, p)
  }
  
  beta.old <- array(beta0, dim = c(p, 1, length(lambda)))
  Z.old <- array(0, dim = c(n, 1, length(lambda)))
  W.old <- array(0, dim = c(n, n, length(lambda)))
  
  beta <- array(beta0, dim = c(p, 1, length(lambda)))
  Z <- array(0, dim = c(n, 1, length(lambda)))
  W <- array(0, dim = c(n, n, length(lambda)))
  
  R <- list()
  
  nc <- 1:length(lambda)
  it.stop <- rep(0, length(lambda))
  
  # Run until convergence or stop
  counter <- 0
  
  repeat {
    if (verbose)
      print(paste("Performing IRLASSO iter", counter, sep = " "))
    
    # Keep old variables
    beta.old <- beta
    Z.old <- Z
    W.old <- W
    
    # Compute W and Z
    eta <- array(0, dim = c(n, 1, length(lambda)))
    k.p <- array(0, dim = c(n, 1, length(lambda)))
    k.pp <- array(0, dim = c(n, 1, length(lambda)))
    
    if (length(nc) == 0) {
      if (verbose)
        print("No lambda left...")
      break
    }
    
    if (verbose)
      print(paste("Lambda left", paste(nc, collapse = " ")))
    for (m in nc) {
      eta[, , m] <- X %*% beta[, , m] + C
      
      k.p[, , m] <- kappa1(eta[, , m])
      
      k.pp[, , m] <- kappa2(eta[, , m])
      k.pp[, , m][which(k.pp[, , m] < 0.005)] <- 0.005
      
      W[, , m] <- diag(as.vector(k.pp[, , m]))
      
      Z[, , m] <- eta[, , m] + diag(1 / as.vector(k.pp[, , m])) %*% (Y -
                                                                       as.vector(k.p[, , m]))
      
      # Weighted X.W and Z.W on largest component
      X.W <- as.matrix(sqrt(W[, , m]) %*% X)
      Z.W <- as.matrix(sqrt(W[, , m]) %*% Z[, , m])
      
      # Compute coefficients for Z.W ~ X.W
      # No center/scale X.W, Z.W
      # No intercept
      if (is.null(w.lambda))
        w.lambda <- rep(1, ncol(X.W))
      fit.lasso <- glmnet(
        x = X.W,
        y = Z.W,
        family = "gaussian",
        alpha = 1,
        lambda = lambda[m],
        standardize = FALSE,
        intercept = FALSE,
        penalty.factor = w.lambda
      )
      beta[, , m] <- as.numeric(fit.lasso$beta)
      
      # Compute model selection matrix
      if (intercept) {
        s.lasso <- which(as.numeric(fit.lasso$beta)[-1] != 0) # Not the intercept, at most 0 <= s.lasso <= p-1
        R[[m]] <- matrix(0, nrow = p, ncol = length(s.lasso))
        R[[m]][1, 1] <- 1                                    # Always take intercept
        if (length(s.lasso) > 0) {
          for (s in 1:length(s.lasso)) {
            i.lasso <- 1 + s.lasso[s]
            R[[m]][i.lasso, s] <- 1
          }
        }
      } else {
        s.lasso <- which(as.numeric(fit.lasso$beta) != 0)     # No intercept, at most 0 <= s.lasso <= p
        R[[m]] <- matrix(0, nrow = p, ncol = length(s.lasso))
        if (length(s.lasso) > 0) {
          for (s in 1:length(s.lasso)) {
            i.lasso <- s.lasso[s]
            R[[m]][i.lasso, s] <- 1
          }
        }
      }
      
    }
    
    epsilon <- sqrt(apply((beta - beta.old) ^ 2, 3, sum) / apply((beta.old) ^
                                                                   2, 3, sum))
    print(paste("Min Divergence", min(epsilon[nc]), sep = " "))
    
    log.like <- apply(beta, 3, function(v)
      sum((X %*% v) * Y - kappa0(X %*% v)))
    log.like.ratio <- log.like - apply(beta.old, 3, function(v)
      sum((X %*% v) * Y - kappa0(X %*% v)))
    print(paste("Min Loglike ratio", min(log.like.ratio[nc]), sep = " "))
    
    if (sum(is.nan(epsilon[nc])) > 0) {
      nan.stop <- which(is.nan(epsilon))
      if (verbose)
        print(paste("Divergence NaN comps", paste(nan.stop, collapse = " ")))
      for (m in nc) {
        beta[, , m] <- beta.old[, , m]
        Z[, , m] <- Z.old[, , m]
        W[, , m] <- W.old[, , m]
      }
      nc <- setdiff(nc, nan.stop)
    }
    
    if ((min(epsilon[nc]) < tol) | (min(log.like.ratio[nc]) < tol)) {
      nc.stop <- which((epsilon < tol) | (log.like.ratio < tol))
      it.stop[nc.stop] <- counter
      if (verbose)
        print(paste(
          "Divergence/Loglike stop comps",
          paste(nc.stop, collapse = " ")
        ))
      nc <- setdiff(nc, nc.stop)
    }
    
    if (counter == maxit) {
      if (verbose)
        print("Maximum iterarion, no convergence...")
      it.stop[which(it.stop == 0)] <- counter
      break
    }
    else {
      counter <- counter + 1
    }
    
  }
  
  # Compute coeffs for original variables
  beta.tilde <- beta * 0
  
  for (m in 1:length(lambda)) {
    if (intercept) {
      beta.tilde[, , m][1] <- sd.Y[1] * beta[, , m][1] + mu.Y[1] - sd.Y[1] * (mu.X /
                                                                                sd.X) %*% beta[, , m][2:p] # Intercept
      beta.tilde[, , m][2:p] <- sd.Y[1] * beta[, , m][2:p] / sd.X                                         # Core
      
    } else {
      for (mc in 1:m) {
        beta.tilde[, , m][1:p] <- sd.Y[1] * beta[, , m][1:p] / sd.X                                       # No intercept, Core
      }
    }
  }
  
  return(list(
    BETA = beta.tilde,
    beta = beta,
    Z = Z,
    W = W,
    R = R,
    it = it.stop
  ))
}


### CONTRIBUTION FUNCTIONS 

# In the followings, l1, l2, l3, l4 will be the number of levels of each factor

mains_contribution <- function(X,
                               beta_main,
                               l1 = 21,
                               l2 = 14,
                               l3 = 2,
                               l4 = 3)
{
  range_main <- unlist(get_ranges4(l1, l2, l3, l4)[1])
  mains_contrib <- X[, range_main] %*% beta_main
  return(mains_contrib)
}


two_ways_contribution <- function(X,
                                  gamma_vec,
                                  beta_vec_2way,
                                  l1 = 21,
                                  l2 = 14,
                                  l3 = 2,
                                  l4 = 3,
                                  already_multiplied = FALSE)
{
  if (already_multiplied == TRUE)
  {
    gamma_vec <- array(1, dim = length(gamma_vec))
  }
  range_2ways <- unlist(get_ranges4(l1, l2, l3, l4)[2])
  two_ways_contrib <- X[, range_2ways] %*% (beta_vec_2way * gamma_vec) ##last multiplication should be elementwise
  return(two_ways_contrib)
}


three_ways_contribution <- function(X,
                                    delta_vec,
                                    beta_vec_3way,
                                    l1 = 21,
                                    l2 = 14,
                                    l3 = 2,
                                    l4 = 3,
                                    already_multiplied = FALSE)
{
  if (already_multiplied == TRUE)
  {
    delta_vec <- array(1, dim = length(delta_vec))
  }
  range_3ways <- unlist(get_ranges4(l1, l2, l3, l4)[3])
  three_ways_contrib <- X[, range_3ways] %*% (beta_vec_3way * delta_vec) ##last multiplication should be elementwise
  return(three_ways_contrib)
}


four_ways_contribution <- function(X,
                                   tau_vec,
                                   beta_vec_4way,
                                   l1 = 21,
                                   l2 = 14,
                                   l3 = 2,
                                   l4 = 3,
                                   already_multiplied = FALSE)
  ##assumes gamma_vec is in same order as X and beta_vec_3way only prod of beta2way
{
  if (already_multiplied == TRUE)
  {
    tau_vec <- array(1, dim = length(tau_vec))
  }
  range_4ways <- unlist(get_ranges4(l1, l2, l3, l4)[4])
  four_ways_contrib <- (X[, range_4ways] %*%  (matrix(tau_vec, ncol = 1) *
                                                 matrix(beta_vec_4way, ncol = 1))) ##last multiplication should be elementwise
  return(four_ways_contrib)
}


# penalty for 1 vector
get_penalty <- function(vector,
                        weights,
                        lambda,
                        already_weighted = TRUE) {
  result = lambda * sum(abs(vector) * abs(weights))
  return(result)
}


# Q BERN FUNCTION -computes the penalized log likelihood for continuous Bernoulli
Q_bern <- function(X,
                   y,
                   beta,
                   gamma_vec,
                   delta_vec,
                   tau_vec,
                   lambda_beta,
                   lambda_gamma,
                   lambda_delta,
                   lambda_tau,
                   w_beta,
                   w_gamma = 1,
                   w_delta = 1,
                   w_tau = 1,
                   l1 = 21,
                   l2 = 14,
                   l3 = 2,
                   l4 = 3,
                   scaled = TRUE,
                   intercept = 0)
{
  if (length(beta) == l1 + l2 + l3 + l4)
    # Give beta for main and compute for the rest given gamma delta tau
  {
    beta_2way <- get_beta_vec_2way4(
      beta = beta,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      gamma = gamma_vec,
      only_beta = FALSE
    )
    beta_3way <- get_beta_vec_3way4(
      beta_2way = beta_2way,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      delta = delta_vec,
      only_beta = FALSE
    )
    beta_4way <- get_beta_vec_4way4(
      beta_3way = beta_3way,
      tau = tau_vec,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      only_beta = FALSE
    )
    beta <- c(beta, beta_2way, beta_3way, beta_4way)
  }
  
  penalty_beta <- get_penalty(vector = beta[unlist(get_ranges4(l1, l2, l3, l4)[1])],
                              weights = w_beta,
                              lambda = lambda_beta)
  penalty_gamma <- get_penalty(vector = gamma_vec,
                               weights = 1,
                               lambda = lambda_gamma)
  penalty_delta <- get_penalty(vector = delta_vec,
                               weights = 1,
                               lambda = lambda_delta)
  penalty_tau <- get_penalty(vector = tau_vec,
                             weights = 1,
                             lambda = lambda_tau)
  
  v = X %*% beta + intercept
  log.like <- sum(y * v - kappa0(v))
  if (scaled == TRUE)
    # be consistent with glmnet
  {
    log.like <- log.like / (2 * dim(X)[1])
  }
  loss <- -log.like + penalty_beta + penalty_gamma + penalty_delta + penalty_tau
  return(loss)
}


# Function to find the minimum for 1 dimensional update of the hierarchical constant for three-way interactions
# interval should depend on old estimator
# C is the offset term
minimizer_Q_bern_delta <- function(X,
                                   y,
                                   C,
                                   lambda,
                                   beta_old,
                                   weight = 1,
                                   scaled = TRUE)
  
{
  #The function that we want to minimize
  fct <- function(b)
  {
    penalty <- abs(b) * lambda * weight
    v = X * b + C
    log.like <- sum(y * v - kappa0(v))
    if (scaled == TRUE)
    {
      log.like <- log.like / (2 * length(X))
    }
    loss <- -log.like + penalty
    return(loss)
  }
  
  #interval for reasonable value of beta
  interval <- c(
    min(-beta_old / 2 - 5e-1, 5 * beta_old / 2 - 5e-1),
    max(-beta_old / 2 + 5e-1, 5 * beta_old / 2 + 5e-1)
  )
  
  result_optimize <- optimize(fct, interval = interval)
  minimum <- result_optimize$minimum
  
  # Check whether new result improves upon 0
  f_0 <- fct(0)
  if (f_0 <= fct(minimum) & f_0 <= fct(beta_old))
  {
    return(0)
  }
  
  
  # Keep 0s when improving is insignificant
  if ( f_0-fct(minimum) <=  10^(-12)*abs(fct(minimum))  &   f_0-fct(beta_old) <=  10^(-12)*abs(fct(minimum))  & beta_old==0)
  {return(0)}
  
  
  #Check whether new estimate improves upon the old one
  if (fct(beta_old) <= fct(minimum))
  {
    return(beta_old)
  }
  
  return(minimum)
}



# Similar function to find the minimum for 1 dimensional update of the hierarchical constant for two-way interactions
minimizer_Q_bern_gamma <- function(X,
                                   Z,
                                   y,
                                   C,
                                   lambda,
                                   beta_old,
                                   weight = 1,
                                   scaled = TRUE)
{
  fct <- function(b)
  {
    penalty <- abs(b) * lambda * weight
    v = X * b + Z * (b ^ 2) + C
    
    log.like <- sum(y * v - kappa0(v))
    if (scaled == TRUE)
    {
      log.like <- log.like / (2 * length(X))
    }
    loss <- -log.like + penalty
    return(loss)
  }
  
  interval <- c(
    min(-beta_old / 2 - 5e-1, 5 * beta_old / 2 - 5e-1),
    max(-beta_old / 2 + 5e-1, 5 * beta_old / 2 + 5e-1)
  )
  result_optimize <- optimize(fct, interval = interval)
  minimum <- result_optimize$minimum
  
  f_0 <- fct(0)
  if (f_0 <= fct(minimum) & f_0 <= fct(beta_old))
  {
    return(0)
  }
  
  # Keep 0s when improving is insignificant
  if ( f_0-fct(minimum) <=  10^(-9)*abs(fct(minimum))  &   f_0-fct(beta_old) <=  10^(-9)*abs(fct(minimum))  & beta_old==0)
  {return(0)}
  
  
  if (fct(beta_old) <= fct(minimum))
  {
    return(beta_old)
  }
  
  return(minimum)
}


# Function to find the minimum for 1 dimensional update of main effects
minimizer_Q_bern_beta <- function(X,
                                  Z,
                                  t,
                                  y,
                                  C,
                                  lambda,
                                  beta_old,
                                  weight = 1,
                                  scaled = TRUE)
{
  fct <- function(b)
  {
    penalty <- weight * abs(b) * lambda
    v = X * b + Z * (b ^ 2) + t * (b ^ 6) + C
    log.like <- sum(y * v - kappa0(v))
    if (scaled == TRUE)
    {
      log.like <- log.like / (2 * length(X))
    }
    loss <- -log.like + penalty
    return(loss)
  }
  
  interval <- c(min(-beta_old / 5 - 5e-1, 2 * beta_old - 5e-1),
                max(-beta_old / 5 + 5e-1, 2 * beta_old + 5e-1))
  result_optimize <- optimize(fct, interval = interval)
  minimum <- result_optimize$minimum
  
  f_0 <- fct(0)
  if (f_0 <= fct(minimum) & f_0 <= fct(beta_old))
  {
    return(0)
  }
  
  # Keep 0s when improving is insignificant
  if ( f_0-fct(minimum) <=  10^(-6)*abs(fct(minimum))  &   f_0-fct(beta_old) <=  10^(-6)*abs(fct(minimum))  & beta_old==0)
  {return(0)}
  
  if (fct(beta_old) <= fct(minimum))
  {
    return(beta_old)
  }
  
  return(minimum)
}
