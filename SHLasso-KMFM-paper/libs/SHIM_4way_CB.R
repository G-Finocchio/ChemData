source("Updates.R")

# SHIM CLASS for 4-way interactions and adaptive estimation of main effects
SHIM_4way <- function(X,
                      y,
                      beta_init,
                      gamma_init,
                      delta_init,
                      tau_init,
                      l1 = 21,
                      l2 = 14,
                      l3 = 2,
                      l4 = 3)
  
{
  self = list()
  
  self$beta_hat <- beta_init
  self$gamma_hat <- gamma_init
  self$delta_hat <- delta_init
  self$tau_hat <- tau_init
  
  self$means_X <- colMeans(X)
  self$stds_X <- apply(X, 2, sd)
  self$mean_y <- mean(y)
  self$scale = scale
  
  self$l1 = l1
  self$l2 = l2
  self$l3 = l3
  self$l4 = l4
  
  
  fit <- function(X,
                  y,
                  lambda_beta,
                  lambda_gamma,
                  lambda_delta,
                  lambda_tau,
                  w_beta = NULL,
                  w_gamma = 1,
                  w_delta = 1,
                  w_tau = 1,
                  tol = 5e-3,
                  max_iter = 10,
                  compute_Q = Q_bern,
                  intercept = 0,
                  use_intercept = TRUE) #Take care! Usually intercept has to be used!
  {
    # STEP 0 (INITIALIZE AND UPDATE THE INTERCEPT)
    beta_hat <- self$beta_hat
    gamma_hat <- self$gamma_hat
    delta_hat <- self$delta_hat
    tau_hat <- self$tau_hat
    Q_old <- 1e100
    beta_2way <- get_beta_vec_2way4(
      beta = beta_hat,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      gamma = gamma_hat,
      only_beta = FALSE
    )
    beta_3way <- get_beta_vec_3way4(
      beta_2way = beta_2way,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      delta = delta_hat,
      only_beta = FALSE
    )
    beta_4way <- get_beta_vec_4way4(
      beta_3way = beta_3way,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      tau = tau_hat,
      only_beta = FALSE
    )
    beta_all <- c(beta_hat, beta_2way, beta_3way, beta_4way)
    
    # Do the updates until convergence or maximum number of iterations is reached
    for (i in c(1:max_iter))
    { 
      cat(" Start iteration: ", i, " ; ")
      
      if (use_intercept == TRUE)
      {
        # Update the intercept
        intercept <- update_intercept(
          X = X,
          y = y,
          beta_all = beta_all,
          intercept_old = intercept
        ) 
      }
      
      # STEP 1 (UPDATE TAU HAT)
      tau_hat <- update_tau(
        X = X,
        y = y,
        beta_hat = beta_hat,
        gamma_hat = gamma_hat,
        delta_hat = delta_hat,
        tau_hat = tau_hat,
        lambda_tau = lambda_tau,
        l1 = self$l1,
        l2 = self$l2,
        l3 = self$l3,
        l4 = self$l4,
        intercept = intercept
      )
      
      # STEP 2 (UPDATE DELTA HAT)
      delta_hat <- update_delta(
        X = X,
        y = y,
        beta_hat = beta_hat,
        gamma_hat = gamma_hat,
        delta_hat = delta_hat,
        tau_hat = tau_hat,
        lambda_delta = lambda_delta,
        l1 = self$l1,
        l2 = self$l2,
        l3 = self$l3,
        l4 = self$l4,
        intercept = intercept
      )
      
      # STEP 3 (UPDATE GAMMA HAT)
      gamma_hat <- update_gamma(
        X = X,
        y = y,
        beta_hat = beta_hat,
        gamma_hat = gamma_hat,
        delta_hat = delta_hat,
        tau_hat = tau_hat,
        lambda_gamma = lambda_gamma,
        l1 = self$l1,
        l2 = self$l2,
        l3 = self$l3,
        l4 = self$l4,
        intercept = intercept
      )
      
      # STEP 4 (UPDATE BETA HAT)
      beta_hat <- update_beta(
        X = X,
        y = y,
        beta_hat = beta_hat,
        gamma_hat = gamma_hat,
        delta_hat = delta_hat,
        tau_hat = tau_hat,
        lambda_beta = lambda_beta,
        intercept = intercept,
        l1 = self$l1,
        l2 = self$l2,
        l3 = self$l3,
        l4 = self$l4,
        w = w_beta
      )
      
      
      # Keep track of the updates
      beta_2way <- get_beta_vec_2way4(beta = beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta = FALSE) 
      beta_3way <- get_beta_vec_3way4(beta_2way = beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta = FALSE) 
      beta_4way <- get_beta_vec_4way4(beta_3way = beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta = FALSE) 
      beta_all<-c(beta_hat, beta_2way, beta_3way, beta_4way)
      
      
      # STEP 5 (COMPUTE RELATIVE DIFFERENCE AND DECIDE IF THE ALGORITHM STOPS)
      Q_new <- compute_Q(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        lambda_delta = lambda_delta,
        lambda_tau = lambda_tau,
        w_beta = w_beta,
        w_gamma = w_gamma,
        w_delta = w_delta,
        w_tau = w_tau,
        l1 = self$l1,
        l2 = self$l2,
        l3 = self$l3,
        l4 = self$l4,
        intercept = intercept
      )
      if (Q_new == Q_old)
        # If equal then stop. This check ensures that the algorithm stops when Q_new=Q_old=0.
      {
        self$beta_hat <- beta_hat
        self$gamma_hat <- gamma_hat
        self$delta_hat <- delta_hat
        self$tau_hat <- tau_hat
        beta_2way <- get_beta_vec_2way4(
          beta = beta_hat,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          gamma = gamma_hat,
          only_beta = FALSE
        )
        beta_3way <- get_beta_vec_3way4(
          beta_2way = beta_2way,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          delta = delta_hat,
          only_beta = FALSE
        )
        beta_4way <- get_beta_vec_4way4(
          beta_3way = beta_3way,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          tau = tau_hat,
          only_beta = FALSE
        )
        beta_all <- c(beta_hat, beta_2way, beta_3way, beta_4way)
        self$intercept <- intercept
        print("The model has converged.")
        return (
          list(
            "beta_hat" = self$beta_hat,
            "gamma_hat" = self$gamma_hat ,
            "delta_hat" = self$delta_hat,
            "tau_hat" = tau_hat,
            "beta_all" = beta_all,
            "intercept" = self$intercept
          )
        )
      }
      # If small enough then stop
      rel_dif <- compute_rel_dif(Q_old = Q_old, Q_new = Q_new)
      cat("  Relative difference is ", rel_dif)
      if (abs(rel_dif) <= tol) {
        self$beta_hat <- beta_hat
        self$gamma_hat <- gamma_hat
        self$delta_hat <- delta_hat
        self$tau_hat <- tau_hat
        beta_2way <- get_beta_vec_2way4(
          beta = beta_hat,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          gamma = gamma_hat,
          only_beta = FALSE
        )
        beta_3way <- get_beta_vec_3way4(
          beta_2way = beta_2way,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          delta = delta_hat,
          only_beta = FALSE
        )
        beta_4way <- get_beta_vec_4way4(
          beta_3way = beta_3way,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          tau = tau_hat,
          only_beta = FALSE
        )
        beta_all <- c(beta_hat, beta_2way, beta_3way, beta_4way)
        self$intercept <- intercept
        print("The model has converged.")
        return (
          list(
            "beta_hat" = self$beta_hat,
            "gamma_hat" = self$gamma_hat ,
            "delta_hat" = self$delta_hat,
            "tau_hat" = tau_hat,
            "beta_all" = beta_all,
            "intercept" = self$intercept
          )
        )
      }
      Q_old <- Q_new #Update Q_old
    }
    cat(
      "It has not converged. The max number of iterations can be increased."
    )
    self$beta_hat <- beta_hat
    self$gamma_hat <- gamma_hat
    self$delta_hat <- delta_hat
    self$tau_hat <- tau_hat
    beta_2way <- get_beta_vec_2way4(
      beta = beta_hat,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      gamma = gamma_hat,
      only_beta = FALSE
    )
    beta_3way <- get_beta_vec_3way4(
      beta_2way = beta_2way,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      delta = delta_hat,
      only_beta = FALSE
    )
    beta_4way <- get_beta_vec_4way4(
      beta_3way = beta_3way,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      tau = tau_hat,
      only_beta = FALSE
    )
    beta_all <- c(beta_hat, beta_2way, beta_3way, beta_4way)
    self$intercept <- intercept
    return (
      list(
        "beta_hat" = self$beta_hat,
        "gamma_hat" = self$gamma_hat ,
        "delta_hat" = self$delta_hat,
        "tau_hat" = tau_hat,
        "beta_all" = beta_all,
        "intercept" = self$intercept
      )
    )
  }
  
  predict <- function(self, X_new)
  {
    beta_all <- self$beta_all
    v <-  X_new %*% beta_all + self$intercept
    y_pred <- kappa1(v)
    return(y_pred)
  }
  
  R2_score <- function(self, X_new, y_true, verbose = TRUE)
  {
    y_pred <- predict(self, X_new, scale = scale)
    if (verbose == TRUE)
    {
      cat("R2 score is", r2(y_true, y_pred), "\n")
      plot(array(y_pred), array(y_true),
           xlab = "Predicted Values",
           ylab = "True Values",
           main = "True vs Predicted Values")
      
    }
    
    return(r2(y_true, y_pred))
  }
  
  return(list(
    fit = fit,
    predict = predict,
    R2_score = R2_score,
    self = self
  ))
  
}
