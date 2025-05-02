source("Functions_for_updates.R")

################## Update the vector of hierarchical constants for the four-way interactions ##################

update_tau <- function(X,
                       y,
                       beta_hat,
                       gamma_hat,
                       delta_hat,
                       tau_hat,
                       lambda_tau,
                       l1 = 21,
                       l2 = 14,
                       l3 = 2,
                       l4 = 3,
                       intercept = 0)
{
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
    only_beta = TRUE
  )
  X_4way <- X[, unlist(get_ranges4(l1, l2, l3, l4)[4])]
  if (var(beta_4way) == 0)
    # Lasso does not work if predictor has variance 0
  {
    print("Four ways are 0.")
    return(beta_4way * 0)
  }
  X_tilde <- matrix(rep(beta_4way, each = nrow(X_4way)), nrow = nrow(X_4way)) *
    X_4way
  X_c <- X[, c(unlist(get_ranges4(l1, l2, l3, l4)[1]),
               unlist(get_ranges4(l1, l2, l3, l4)[2]),
               unlist(get_ranges4(l1, l2, l3, l4)[3]))] # Xc
  beta_c <- c(beta_hat, beta_2way, beta_3way)
  C <- X_c %*% beta_c + intercept # add intercept to C
  tau_hat_old <- tau_hat
  
  Q_old <- Q_bern(
    X = X,
    y = y,
    beta = beta_hat,
    gamma_vec = gamma_hat,
    delta_vec = delta_hat,
    tau_vec = tau_hat,
    lambda_beta = 0,
    lambda_gamma = 0,
    lambda_delta = 0,
    lambda_tau = lambda_tau,
    w_beta = 1,
    w_gamma = 1,
    w_delta = 1,
    w_tau = 1,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4,
    intercept = intercept
  )
  
  lasso_rez <- irlasso.cb(
    X = X_tilde,
    Y = y,
    lambda = lambda_tau,
    w.lambda = NULL,
    centering = FALSE,
    scaling = FALSE,
    intercept = F,
    maxit = 10,
    tol = 0.0545,
    sd.tol = 1e-6,
    verbose = F,
    C = C
  )
  tau_hat <- array(lasso_rez$BETA, dim = length(lasso_rez$BETA))
  
  Q_new <- Q_bern(
    X = X,
    y = y,
    beta = beta_hat,
    gamma_vec = gamma_hat,
    delta_vec = delta_hat,
    tau_vec = tau_hat,
    lambda_beta = 0,
    lambda_gamma = 0,
    lambda_delta = 0,
    lambda_tau = lambda_tau,
    w_beta = 1,
    w_gamma = 1,
    w_delta = 1,
    w_tau = 1,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4,
    intercept = intercept
  )
  
  if (Q_new - Q_old > abs(Q_old) * 1e-5) {
    tau_hat <- tau_hat_old
    print(
      "There might be numerical instability in update tau which was taken care of by using old tau. "
    )
  }
  if (Q_new - Q_old >= 0) {
    tau_hat <- tau_hat_old
  }
  print("Updated tau")
  
  return(tau_hat)
  
}


####################################### Update the intercept #########################################
update_intercept <- function(X, y, beta_all, intercept_old)
  #function to find the minimum for intercept
{
  fctint <- function(intercept)
  {
    v = as.matrix(X %*% beta_all + intercept)
    log.like <- sum(y * v - kappa0(v))
    loss <- -log.like #scale does not matter
    return(loss)
  }
  interval <- c(
    min(-intercept_old / 2 - 7, 5 * intercept_old / 2 - 7),
    max(-intercept_old / 2 + 7, 5 * intercept_old / 2 + 7)
  )
  result_optimize <- optimize(fctint, interval = interval)
  minimum <- result_optimize$minimum
  f_0 <- fctint(0)
  if (f_0 <= fctint(minimum) & f_0 <= fctint(intercept_old))
  {
    return(0)
  }
  
  if (fctint(intercept_old) <= fctint(minimum))
  {
    return(intercept_old)
  }
  print("Updated intercept")
  
  return(minimum)
}



################## Update the vector of hierarchical constants for the two-way interactions ##################
update_gamma <- function(X,
                         y,
                         beta_hat,
                         gamma_hat,
                         delta_hat,
                         tau_hat,
                         lambda_gamma,
                         l1 = 21,
                         l2 = 14,
                         l3 = 2,
                         l4 = 3,
                         w = 1,
                         intercept = 0)
{
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  X_main <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[1])]
  X_2way <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[2])]
  X_3way <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[3])]
  X_4way <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[4])]
  
  if (w == 1)
  {
    w = array(1, dim = length(gamma_hat))
  }
  
  # We will do the updates separately for all possible combinations of different factors i.e., f1f2, f1f3, ..., f1f4
  
  # i,j
  for (i in range1) {
    for (j in range2) {
      if (beta_hat[i] == 0 || beta_hat[j] == 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(c(i, j), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
      } else {
      discard_from_c_2way <- c(get_position_vec_from_theta_matrix4(
        c(i, j),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      ))
      
      discard_from_c_3way <- c()
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      # Contribution of 2-ways without gamma
      two_ways <- X_2way[, get_position_vec_from_theta_matrix4(
        c(i, j),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[i] * beta_hat[j]
      three_ways = 0
      
      # Contribution of 3-ways without gamma
      for (k in c(range3, range4))
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[j] * beta_hat[k]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      
      
      discard_from_c_4way <- c()
      four_ways = 0
      # Contribution of 4-ways without gamma
      for (k in range3)
      {
        for (l in range4)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, k),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
      
      X_tilde <- two_ways +  three_ways
      Z_tilde <- four_ways
      
      
      C <- intercept + X_main %*% beta_hat + X_2way[, -discard_from_c_2way] %*%
        beta_2way[-discard_from_c_2way] +
        X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
        X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] # C is the offset term
      
      
      
      Q_old <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      
      
      gamma_old_value <- gamma_hat[get_position_vec_from_theta_matrix4(
        c(i, j),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )]
      # Use the minimizer
      gamma_hat[get_position_vec_from_theta_matrix4(
        c(i, j),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )] <- minimizer_Q_bern_gamma(
        X = X_tilde,
        Z = Z_tilde,
        y = y,
        C = C,
        lambda = lambda_gamma,
        beta_old = gamma_old_value,
        weight = 1,
        scaled = TRUE
      ) # intercept is in C
      
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      Q_new <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      if (Q_new - Q_old > abs(Q_old) * 1e-2) {
        print("There might be numerical instability in update gamma, which was taken care of.")
      }
      # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
      if (Q_new - Q_old >= 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(i, j),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- gamma_old_value
      }
    }
  }
  }
    
  # i,k
  for (i in range1) {
    for (k in range3) {
      
      if (beta_hat[i] == 0 || beta_hat[k] == 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(c(i, k), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
      } else {
      discard_from_c_2way <- c(get_position_vec_from_theta_matrix4(
        c(i, k),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      ))
      
      discard_from_c_3way <- c()
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      # Contribution of 2-ways without gamma
      two_ways <- X_2way[, get_position_vec_from_theta_matrix4(
        c(i, k),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[i] * beta_hat[k]
      # Contribution of 3-ways without gamma
      three_ways = 0
      for (j in range2)
        # Compute 3 ways ijk
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[j] * beta_hat[k]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      for (l in range4)
        # Compute 3 ways ikl
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[k] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      
      discard_from_c_4way <- c()
      # Contribution of 4-ways without gamma
      four_ways = 0
      for (j in range2)
      {
        for (l in range4)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, j),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
      
      X_tilde <- two_ways +  three_ways
      Z_tilde <- four_ways
      
      C <- intercept + X_main %*% beta_hat + X_2way[, -discard_from_c_2way] %*%
        beta_2way[-discard_from_c_2way] +
        X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
        X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] # C is the offest term
      
      Q_old <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      # Use the minimizer
      gamma_old_value <- gamma_hat[get_position_vec_from_theta_matrix4(
        c(i, k),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )]
      
      gamma_hat[get_position_vec_from_theta_matrix4(
        c(i, k),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )] <- minimizer_Q_bern_gamma(
        X = X_tilde,
        Z = Z_tilde,
        y = y,
        C = C,
        lambda = lambda_gamma,
        beta_old = gamma_old_value,
        weight = 1,
        scaled = TRUE
      ) ## Intercept is in C
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      Q_new <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      if (Q_new - Q_old > abs(Q_old) * 1e-2) {
        print("There might be numerical instability in update gamma which was taken care of.")
      }
      # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
      if (Q_new - Q_old >= 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(i, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- gamma_old_value
      }
    }
  }
  }
  
  # i,l
  for (i in range1) {
    for (l in range4) {
      if (beta_hat[i] == 0 || beta_hat[l] == 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(c(i, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
      } else {
      discard_from_c_2way <- c(get_position_vec_from_theta_matrix4(
        c(i, l),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      ))
      
      discard_from_c_3way <- c()
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      # Contribution of 2-ways without gamma
      two_ways <- X_2way[, get_position_vec_from_theta_matrix4(
        c(i, l),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[i] * beta_hat[l]
      # Contribution of 3-ways without gamma
      three_ways = 0
      
      for (j in range2)
        # Compute 3 ways ijk
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[j] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      for (k in range3)
        # Compute 3 ways ikl
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[k] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      
      discard_from_c_4way <- c()
      # Contribution of 4-ways without gamma
      four_ways = 0
      for (j in range2)
      {
        for (k in range3)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, j),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
      
      X_tilde <- two_ways +  three_ways
      Z_tilde <- four_ways
      
      C <- intercept + X_main %*% beta_hat + X_2way[, -discard_from_c_2way] %*%
        beta_2way[-discard_from_c_2way] +
        X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
        X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] # C is the offset term
      
      Q_old <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      gamma_old_value <- gamma_hat[get_position_vec_from_theta_matrix4(
        c(i, l),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )]
      #use the minimizer
      gamma_hat[get_position_vec_from_theta_matrix4(
        c(i, l),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )] <- minimizer_Q_bern_gamma(
        X = X_tilde,
        Z = Z_tilde,
        y = y,
        C = C,
        lambda = lambda_gamma,
        beta_old = gamma_old_value,
        weight = 1,
        scaled = TRUE
      )
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      Q_new <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      if (Q_new - Q_old > abs(Q_old) * 1e-2) {
        print("There might be numerical instability in update gamma.")
      }
      # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
      if (Q_new - Q_old >= 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(i, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- gamma_old_value
      }
    }
    }
  }
  
  # j,k
  for (j in range2) {
    for (k in range3) {
      if (beta_hat[j] == 0 || beta_hat[k] == 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(c(j, k), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
      } else {
      discard_from_c_2way <- c(get_position_vec_from_theta_matrix4(
        c(j, k),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      ))
      
      discard_from_c_3way <- c()
      
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
      ) # This is  WITH tau
      beta_all <- c(beta_hat, beta_2way, beta_3way, beta_4way)
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      # Contribution of 2-ways without gamma
      two_ways <- X_2way[, get_position_vec_from_theta_matrix4(
        c(j, k),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[j] * beta_hat[k]
      # Contribution of 3-ways without gamma
      three_ways = 0
      for (i in range1)
        # Compute 3 ways ijk
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[j] * beta_hat[k]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      for (l in range4)
        # Compute 3 ways ikl
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(j, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[j] * beta_hat[k] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      
      discard_from_c_4way <- c()
      # Contribution of 4-ways without gamma
      four_ways = 0
      for (i in range1)
      {
        for (l in range4)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, j),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
      
      X_tilde <- two_ways +  three_ways
      Z_tilde <- four_ways
      
      C <- intercept + X_main %*% beta_hat + X_2way[, -discard_from_c_2way] %*%
        beta_2way[-discard_from_c_2way] +
        X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
        X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] # C is the offset term
      
      Q_old <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      # Use the minimizer
      gamma_old_value <- gamma_hat[get_position_vec_from_theta_matrix4(
        c(j, k),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )]
      gamma_hat[get_position_vec_from_theta_matrix4(
        c(j, k),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )] <- minimizer_Q_bern_gamma(
        X = X_tilde,
        Z = Z_tilde,
        y = y,
        C = C,
        lambda = lambda_gamma,
        beta_old = gamma_old_value,
        weight = 1,
        scaled = TRUE
      )
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      Q_new <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      if (Q_new - Q_old > abs(Q_old) * 1e-2) {
        print("There might be numerical instability in update gamma which was taken care of.")
      }
      # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
      if (Q_new - Q_old >= 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- gamma_old_value
      }
    }
  }
  }
  # j,l
  for (j in range2) {
    for (l in range4) {
      if (beta_hat[j] == 0 || beta_hat[l] == 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(c(j, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
      } else {
      discard_from_c_2way <- c(get_position_vec_from_theta_matrix4(
        c(j, l),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      ))
      
      discard_from_c_3way <- c()
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      # Contribution of 2-ways without gamma
      two_ways <- X_2way[, get_position_vec_from_theta_matrix4(
        c(j, l),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[j] * beta_hat[l]
      # Contribution of 3-ways without gamma
      three_ways = 0
      
      for (i in range1)
        # Compute 3 ways ijl
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[j] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      for (k in range3)
        # Compute 3 ways jkl
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(j, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[j] * beta_hat[k] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      
      discard_from_c_4way <- c()
      # Contribution of 4-ways without gamma
      four_ways = 0
      for (i in range1)
        #
      {
        for (k in range3)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, j),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
      
      X_tilde <- two_ways +  three_ways
      Z_tilde <- four_ways
      
      C <- intercept + X_main %*% beta_hat + X_2way[, -discard_from_c_2way] %*%
        beta_2way[-discard_from_c_2way] +
        X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
        X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] # C is the offset term
      Q_old <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      # Use the minimizer
      gamma_old_value <- gamma_hat[get_position_vec_from_theta_matrix4(
        c(j, l),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )]
      gamma_hat[get_position_vec_from_theta_matrix4(
        c(j, l),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )] <- minimizer_Q_bern_gamma(
        X = X_tilde,
        Z = Z_tilde,
        y = y,
        C = C,
        lambda = lambda_gamma,
        beta_old = gamma_old_value,
        weight = 1,
        scaled = TRUE
      ) # Intercept is in C
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      Q_new <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      if (Q_new - Q_old > abs(Q_old) * 1e-2) {
        print("There might be numerical instability in update gamma.")
      }
      # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
      if (Q_new - Q_old >= 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(j, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- gamma_old_value
      }
    }
  }
  }
  # k,l
  for (k in range3) {
    for (l in range4) {
      if (beta_hat[k] == 0 || beta_hat[l] == 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(c(k, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
      } else {
      discard_from_c_2way <- c(get_position_vec_from_theta_matrix4(
        c(k, l),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      ))
      
      discard_from_c_3way <- c()
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      # Contribution of 2-ways without gamma
      two_ways <- X_2way[, get_position_vec_from_theta_matrix4(
        c(k, l),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[k] * beta_hat[l]
      
      # Contribution of 3-ways without gamma
      three_ways = 0
      for (i in range1)
        #compute 3 ways ijl
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[k] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      for (j in range2)
        #compute 3 ways jkl
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(j, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[j] * beta_hat[k] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
      
      discard_from_c_4way <- c()
      # Contribution of 4-ways without gamma
      four_ways = 0
      for (i in range1)
        #compute 4 ways contrib
      {
        for (j in range2)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, j),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
      
      X_tilde <- two_ways +  three_ways
      Z_tilde <- four_ways
      
      C <- intercept + X_main %*% beta_hat + X_2way[, -discard_from_c_2way] %*%
        beta_2way[-discard_from_c_2way] +
        X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
        X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] #C is the offset term
      
      Q_old <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      gamma_old_value <- gamma_hat[get_position_vec_from_theta_matrix4(
        c(k, l),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )]
      # Use the minimizer
      gamma_hat[get_position_vec_from_theta_matrix4(
        c(k, l),
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4
      )] <- minimizer_Q_bern_gamma(
        X = X_tilde,
        Z = Z_tilde,
        y = y,
        C = C,
        lambda = lambda_gamma,
        beta_old = gamma_old_value,
        weight = 1,
        scaled = TRUE
      ) # Intercept is in C
      
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
      beta_all <- matrix(beta_all, ncol = 1)
      eta <- X %*% beta_all + intercept
      
      Q_new <- Q_bern(
        X = X,
        y = y,
        beta = beta_hat,
        gamma_vec = gamma_hat,
        delta_vec = delta_hat,
        tau_vec = tau_hat,
        lambda_beta = 0,
        lambda_gamma = lambda_gamma,
        lambda_delta = 0,
        lambda_tau = 0,
        w_beta = 1,
        w_gamma = 1,
        w_delta = 1,
        w_tau = 1,
        l1 = l1,
        l2 = l2,
        l3 = l3,
        l4 = l4,
        intercept = intercept
      )
      
      if (Q_new - Q_old > abs(Q_old) * 1e-2) {
        print("There might be numerical instability in update gamma which was taken care of.")
      }
      # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
      if (Q_new - Q_old >= 0) {
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- gamma_old_value
      }
    }
  }
  }
  print("Updated gamma")
  return (gamma_hat)
}


################## Update the vector of hierarchical constants for the three-way interactions ##################
update_delta <- function(X,
                         y,
                         beta_hat,
                         gamma_hat,
                         delta_hat,
                         tau_hat,
                         lambda_delta,
                         l1 = 21,
                         l2 = 14,
                         l3 = 2,
                         l4 = 3,
                         w = 1,
                         intercept = 0)
{
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  X_main <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[1])]
  X_2way <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[2])]
  X_3way <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[3])]
  X_4way <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[4])]
  
  
  
  if (w == 1)
  {
    w = array(1, dim = length(delta_hat))
  }
  
  # We will do the updates separately for all possible combinations of different factors i.e., f1f2f3, f1f2f4, ..., f2f3f4
  
  # i,j,k
  for (i in range1) {
    for (j in range2) {
      for (k in c(range3)) {
        if (
          gamma_hat[get_position_vec_from_theta_matrix4(c(i, j), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0 ||
          gamma_hat[get_position_vec_from_theta_matrix4(c(i, k), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0 ||
          gamma_hat[get_position_vec_from_theta_matrix4(c(j, k), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0
        ) {
          delta_hat[psi_table_position_to_vector_index4(c(i, j, k), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
        } else {
        discard_from_c_3way <- c()
        
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
        beta_all <- matrix(beta_all, ncol = 1)
        eta <- X %*% beta_all + intercept
        
        # Contribution of 3-ways without delta
        three_ways = 0
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[j] * beta_hat[k]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        ) # Actually just one
        assert(
          length(discard_from_c_3way == 1),
          "should discard only one three way in update delta"
        )
        
        discard_from_c_4way <- c()
        # Contribution of 4-ways without delta
        four_ways = 0
        for (l in range4)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, k),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, j),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
        
        X_tilde <- three_ways + four_ways
        C <- X_main %*% beta_hat + X_2way %*% beta_2way++X_3way[, -discard_from_c_3way] %*%
          beta_3way[-discard_from_c_3way] +
          X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] +
          intercept # C is the offset term
        
        Q_old <- Q_bern(
          X = X,
          y = y,
          beta = beta_hat,
          gamma_vec = gamma_hat,
          delta_vec = delta_hat,
          tau_vec = tau_hat,
          lambda_beta = 0,
          lambda_gamma = 0,
          lambda_delta = lambda_delta,
          lambda_tau = 0,
          w_beta = 1,
          w_gamma = 1,
          w_delta = 1,
          w_tau = 1,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          intercept = intercept
        )
        # Use the minimizer
        delta_old_value <- delta_hat[psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )]
        delta_hat[psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- minimizer_Q_bern_delta(
          X = X_tilde,
          y = y,
          C = C,
          lambda = lambda_delta,
          beta_old = delta_old_value,
          weight = 1,
          scaled = TRUE
        ) # Intercept is in C
        
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
        beta_all <- matrix(beta_all, ncol = 1)
        eta <- X %*% beta_all + intercept
        
        Q_new <- Q_bern(
          X = X,
          y = y,
          beta = beta_hat,
          gamma_vec = gamma_hat,
          delta_vec = delta_hat,
          tau_vec = tau_hat,
          lambda_beta = 0,
          lambda_gamma = 0,
          lambda_delta = lambda_delta,
          lambda_tau = 0,
          w_beta = 1,
          w_gamma = 1,
          w_delta = 1,
          w_tau = 1,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          intercept = intercept
        )
        
        if (Q_new - Q_old > abs(Q_old) * 1e-2) {
          print("There might be numerical instability in update delta which was taken care of.")
        }
        # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
        if (Q_new - Q_old >= 0) {
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] <- delta_old_value
       
        }
        
      }
    }
  }
  }
  #i,j,l
  for (i in range1) {
    for (j in range2) {
      for (l in range4) {
        if (
          gamma_hat[get_position_vec_from_theta_matrix4(c(i,j), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0 ||
          gamma_hat[get_position_vec_from_theta_matrix4(c(i,l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0 ||
          gamma_hat[get_position_vec_from_theta_matrix4(c(j,l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0
        ) {
          delta_hat[psi_table_position_to_vector_index4(c(i, j, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
        } else {
        discard_from_c_3way <- c()
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
        beta_all <- matrix(beta_all, ncol = 1)
        eta <- X %*% beta_all + intercept
        
        three_ways = 0
        # Contribution of 3-ways without delta
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[j] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        ) # Actually just one
        assert(
          length(discard_from_c_3way == 1),
          "should discard only one three way in update delta"
        )
        discard_from_c_4way <- c()
        # Contribution of 4-ways without delta
        four_ways = 0
        for (k in range3)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, k),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, j),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
        
        
        X_tilde <- three_ways + four_ways
        C <- X_main %*% beta_hat + X_2way %*% beta_2way++X_3way[, -discard_from_c_3way] %*%
          beta_3way[-discard_from_c_3way] +
          X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] +
          intercept # C is the offset term
        Q_old <- Q_bern(
          X = X,
          y = y,
          beta = beta_hat,
          gamma_vec = gamma_hat,
          delta_vec = delta_hat,
          tau_vec = tau_hat,
          lambda_beta = 0,
          lambda_gamma = 0,
          lambda_delta = lambda_delta,
          lambda_tau = 0,
          w_beta = 1,
          w_gamma = 1,
          w_delta = 1,
          w_tau = 1,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          intercept = intercept
        )
        
        # Use the minimizer
        delta_old_value <- delta_hat[psi_table_position_to_vector_index4(
          c(i, j, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )]
        delta_hat[psi_table_position_to_vector_index4(
          c(i, j, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- minimizer_Q_bern_delta(
          X = X_tilde,
          y = y,
          C = C,
          lambda = lambda_delta,
          beta_old = delta_old_value,
          weight = 1,
          scaled = TRUE
        ) # Intercept is in C
        
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
        beta_all <- matrix(beta_all, ncol = 1)
        eta <- X %*% beta_all + intercept
        
        Q_new <- Q_bern(
          X = X,
          y = y,
          beta = beta_hat,
          gamma_vec = gamma_hat,
          delta_vec = delta_hat,
          tau_vec = tau_hat,
          lambda_beta = 0,
          lambda_gamma = 0,
          lambda_delta = lambda_delta,
          lambda_tau = 0,
          w_beta = 1,
          w_gamma = 1,
          w_delta = 1,
          w_tau = 1,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          intercept = intercept
        )
        
        if (Q_new - Q_old > abs(Q_old) * 1e-2) {
          print("There might be numerical instability in update delta which was taken care of.")
        }
        # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
        if (Q_new - Q_old >= 0) {
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] <- delta_old_value
        }
        
      }
    }
  }
  }
  # i,k,l
  for (i in range1) {
    for (k in range3) {
      for (l in range4) {
        if (
          gamma_hat[get_position_vec_from_theta_matrix4(c(i, k), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0 ||
          gamma_hat[get_position_vec_from_theta_matrix4(c(i, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0 ||
          gamma_hat[get_position_vec_from_theta_matrix4(c(k, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0
        ) {
          delta_hat[psi_table_position_to_vector_index4(c(i, k, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
        } else {
        discard_from_c_3way <- c()
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
        beta_all <- matrix(beta_all, ncol = 1)
        eta <- X %*% beta_all + intercept
        
        
        # Contribution of 3-ways without delta
        three_ways = 0
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[k] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
        assert(
          length(discard_from_c_3way == 1),
          "should discard only one three way in update delta"
        )
        discard_from_c_4way <- c()
        # Contribution of 4-ways without delta
        four_ways = 0
        for (j in range2)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, k),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, j),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
          
        }
        
        X_tilde <- three_ways + four_ways
        C <- X_main %*% beta_hat + X_2way %*% beta_2way++X_3way[, -discard_from_c_3way] %*%
          beta_3way[-discard_from_c_3way] +
          X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] +
          intercept # C is the offset term
        
        Q_old <- Q_bern(
          X = X,
          y = y,
          beta = beta_hat,
          gamma_vec = gamma_hat,
          delta_vec = delta_hat,
          tau_vec = tau_hat,
          lambda_beta = 0,
          lambda_gamma = 0,
          lambda_delta = lambda_delta,
          lambda_tau = 0,
          w_beta = 1,
          w_gamma = 1,
          w_delta = 1,
          w_tau = 1,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          intercept = intercept
        )
        
        #Use the minimizer
        delta_old_value <- delta_hat[psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )]
        delta_hat[psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- minimizer_Q_bern_delta(
          X = X_tilde,
          y = y,
          C = C,
          lambda = lambda_delta,
          beta_old = delta_old_value,
          weight = 1,
          scaled = TRUE
        ) # Intercept is in C
        
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
        beta_all <- matrix(beta_all, ncol = 1)
        eta <- X %*% beta_all + intercept
        
        Q_new <- Q_bern(
          X = X,
          y = y,
          beta = beta_hat,
          gamma_vec = gamma_hat,
          delta_vec = delta_hat,
          tau_vec = tau_hat,
          lambda_beta = 0,
          lambda_gamma = 0,
          lambda_delta = lambda_delta,
          lambda_tau = 0,
          w_beta = 1,
          w_gamma = 1,
          w_delta = 1,
          w_tau = 1,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          intercept = intercept
        )
        
        if (Q_new - Q_old > abs(Q_old) * 1e-2) {
          print("There might be numerical instability in update delta which was taken care of.")
        }
        # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
        if (Q_new - Q_old >= 0) {
          delta_hat[psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] <- delta_old_value
          
          
        }
        
      }
    }
  }
  }
  #j,k,l
  for (j in range2) {
    for (k in range3) {
      for (l in range4) {
        if (
          gamma_hat[get_position_vec_from_theta_matrix4(c(j, k), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0 ||
          gamma_hat[get_position_vec_from_theta_matrix4(c(j, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0 ||
          gamma_hat[get_position_vec_from_theta_matrix4(c(k, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] == 0
        ) {
          delta_hat[psi_table_position_to_vector_index4(c(j, k, l), l1 = l1, l2 = l2, l3 = l3, l4 = l4)] <- 0
        } else {
        discard_from_c_3way <- c()
        beta_2way <- get_beta_vec_2way4(
          beta = beta_hat,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          gamma = gamma_hat,
          only_beta = FALSE
        ) ###This is with delta
        beta_3way <- get_beta_vec_3way4(
          beta_2way = beta_2way,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          delta = delta_hat,
          only_beta = FALSE
        ) #This is with gamma WITH delta
        beta_4way <- get_beta_vec_4way4(
          beta_3way = beta_3way,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          tau = tau_hat,
          only_beta = FALSE
        ) #This is  WITH tau
        beta_all <- c(beta_hat, beta_2way, beta_3way, beta_4way)
        beta_all <- matrix(beta_all, ncol = 1)
        eta <- X %*% beta_all + intercept
        # Contribution of 3-ways without delta
        three_ways = 0
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(j, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[j] * beta_hat[k] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        ) # Actually just one
        assert(
          length(discard_from_c_3way == 1),
          "should discard only one three way in update delta"
        )
        discard_from_c_4way <- c()
        # Contribution of 4-ways without delta
        four_ways = 0
        for (i in range1)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i] * beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, k),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, j),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
          
        }
        
        X_tilde <- three_ways + four_ways
        
        C <- X_main %*% beta_hat + X_2way %*% beta_2way++X_3way[, -discard_from_c_3way] %*%
          beta_3way[-discard_from_c_3way] +
          X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] +
          intercept # C is the offset term
        
        Q_old <- Q_bern(
          X = X,
          y = y,
          beta = beta_hat,
          gamma_vec = gamma_hat,
          delta_vec = delta_hat,
          tau_vec = tau_hat,
          lambda_beta = 0,
          lambda_gamma = 0,
          lambda_delta = lambda_delta,
          lambda_tau = 0,
          w_beta = 1,
          w_gamma = 1,
          w_delta = 1,
          w_tau = 1,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          intercept = intercept
        )
        
        delta_old_value <- delta_hat[psi_table_position_to_vector_index4(
          c(j, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )]
        # Use the minimizer
        delta_hat[psi_table_position_to_vector_index4(
          c(j, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- minimizer_Q_bern_delta(
          X = X_tilde,
          y = y,
          C = C,
          lambda = lambda_delta,
          beta_old = delta_old_value,
          weight = 1,
          scaled = TRUE
        ) # Intercept is in C
        
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
        beta_all <- matrix(beta_all, ncol = 1)
        eta <- X %*% beta_all + intercept
        
        Q_new <- Q_bern(
          X = X,
          y = y,
          beta = beta_hat,
          gamma_vec = gamma_hat,
          delta_vec = delta_hat,
          tau_vec = tau_hat,
          lambda_beta = 0,
          lambda_gamma = 0,
          lambda_delta = lambda_delta,
          lambda_tau = 0,
          w_beta = 1,
          w_gamma = 1,
          w_delta = 1,
          w_tau = 1,
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4,
          intercept = intercept
        )
        
        if (Q_new - Q_old > abs(Q_old) * 1e-2) {
          print("There might be numerical instability in update delta which was taken care of.")
        }
        # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
        if (Q_new - Q_old >= 0) {
          delta_hat[psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] <- delta_old_value
          
        }
        
      }
    }
  }
  }
  print("Updated delta")
  return(delta_hat)
  
}


################## Update the vector of main effects ##################
update_beta <- function(X,
                        y,
                        beta_hat,
                        gamma_hat,
                        delta_hat,
                        tau_hat,
                        lambda_beta,
                        l1 = 21,
                        l2 = 14,
                        l3 = 2,
                        l4 = 3,
                        w = 1,
                        intercept = 0)
{
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  X_main <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[1])]
  X_2way <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[2])]
  X_3way <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[3])]
  X_4way <- X[, unlist(get_ranges4(
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )[4])]
  
  if (length(w) == 1)
  {
    w = array(1, dim = length(beta_hat))
  }
  
  beta_2way <- get_beta_vec_2way4(
    beta = beta_hat,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4,
    gamma = gamma_hat,
    only_beta = FALSE
  ) ###This is with delta
  beta_3way <- get_beta_vec_3way4(
    beta_2way = beta_2way,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4,
    delta = delta_hat,
    only_beta = FALSE
  ) #This is with gamma WITH delta
  beta_4way <- get_beta_vec_4way4(
    beta_3way = beta_3way,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4,
    tau = tau_hat,
    only_beta = FALSE
  ) #This is  WITH tau
  beta_all <- c(beta_hat, beta_2way, beta_3way, beta_4way)
  beta_all <- matrix(beta_all, ncol = 1)
  eta <- X %*% beta_all + intercept
  
  # We will do the updates separately for all 4 factors in the order f1, f2, f3, f4
  
  #i
  for (i in range1) {
    discard_from_c_main <- c(i)
    mains <- matrix(X_main[, i], nrow = length(eta))
    
    
    discard_from_c_2way <- c()
    # Contribution of 2-ways without beta
    two_ways <- 0
    for (j in c(range2, range3, range4)) {
      discard_from_c_2way <- c(
        discard_from_c_2way,
        get_position_vec_from_theta_matrix4(
          c(i, j),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )
      )
      two_ways <- two_ways + X_2way[, get_position_vec_from_theta_matrix4(
        c(i, j),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[j] *
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(i, j),
          l1 = l1,
          l2 = l2 ,
          l3 = l3,
          l4 = l4
        )]
    }
    
    discard_from_c_3way <- c()
    # Contribution of 3-ways without beta
    three_ways = 0
    
    for (j in range2) {
      for (k in c(range3, range4))
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[j] * beta_hat[k]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
    }
    for (k in range3) {
      for (l in  range4)
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[l] * beta_hat[k]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
    }
    
    discard_from_c_4way <- c()
    # Contribution of 4-ways without beta
    four_ways = 0
    for (j in range2) {
      for (k in range3)
      {
        for (l in range4)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[j]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, k),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, j),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
    }
    
    X_tilde <- mains + two_ways
    Z_tilde <- three_ways
    T_tilde <- four_ways
    C <- intercept + X_main[, -i] %*% beta_hat[-i] + X_2way[, -discard_from_c_2way] %*%
      beta_2way[-discard_from_c_2way] +
      X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
      X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] # C is the offset term
    Q_old <- Q_bern(
      X = X,
      y = y,
      beta = beta_hat,
      gamma_vec = gamma_hat,
      delta_vec = delta_hat,
      tau_vec = tau_hat,
      lambda_beta = lambda_beta,
      lambda_gamma = 0,
      lambda_delta = 0,
      lambda_tau = 0,
      w_beta = w[i],
      w_gamma = 1,
      w_delta = 1,
      w_tau = 1,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      intercept = intercept
    )
    
    beta_old_value <- beta_hat[i]
    # Use the minimizer
    beta_hat[i] <- minimizer_Q_bern_beta(
      X = X_tilde,
      Z = Z_tilde,
      t = T_tilde,
      y = y,
      C = C,
      lambda = lambda_beta,
      beta_old = beta_old_value,
      weight = w[i],
      scaled = TRUE
    ) # Intercept is in C
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
    beta_all <- matrix(beta_all, ncol = 1)
    eta <- X %*% beta_all + intercept
    
    Q_new <- Q_bern(
      X = X,
      y = y,
      beta = beta_hat,
      gamma_vec = gamma_hat,
      delta_vec = delta_hat,
      tau_vec = tau_hat,
      lambda_beta = lambda_beta,
      lambda_gamma = 0,
      lambda_delta = 0,
      lambda_tau = 0,
      w_beta = w[i],
      w_gamma = 1,
      w_delta = 1,
      w_tau = 1,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      intercept = intercept
    )
    
    if (Q_new - Q_old > abs(Q_old) * 1e-2) {
      print("There might be numerical instability in update beta which was taken care of.")
    }
    if (Q_new - Q_old >= 0) {
      beta_hat[i] <- beta_old_value
    }
  }
  
  #j
  for (j in range2) {
    discard_from_c_main <- c(j)
    mains <- matrix(X_main[, j], nrow = length(eta))
    
    
    discard_from_c_2way <- c()
    # Contribution of 2-ways without beta
    two_ways <- 0
    for (ik in c(range1, range3, range4)) {
      small_idx = min(ik, j)
      big_idx = max(ik, j)
      
      discard_from_c_2way <- c(
        discard_from_c_2way,
        get_position_vec_from_theta_matrix4(
          c(small_idx, big_idx),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )
      )
      two_ways <- two_ways + X_2way[, get_position_vec_from_theta_matrix4(
        c(small_idx, big_idx),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[ik] *
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(small_idx, big_idx),
          l1 = l1,
          l2 = l2 ,
          l3 = l3,
          l4 = l4
        )]
    }
    
    discard_from_c_3way <- c()
    # Contribution of 3-ways without beta
    three_ways = 0
    for (i in range1) {
      for (k in c(range3, range4))
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[k]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
    }
    for (k in range3) {
      for (l in  range4)
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(j, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[l] * beta_hat[k]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
    }
    
    discard_from_c_4way <- c()
    # Contribution of 2-ways without beta
    four_ways = 0
    for (i in range1) {
      for (k in range3)
      {
        for (l in range4)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, k),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, j),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[k] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
    }
    
    X_tilde <- mains + two_ways
    Z_tilde <- three_ways
    T_tilde <- four_ways
    C <- intercept + X_main[, -j] %*% beta_hat[-j] + X_2way[, -discard_from_c_2way] %*%
      beta_2way[-discard_from_c_2way] +
      X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
      X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] # C is the offset term
    Q_old <- Q_bern(
      X = X,
      y = y,
      beta = beta_hat,
      gamma_vec = gamma_hat,
      delta_vec = delta_hat,
      tau_vec = tau_hat,
      lambda_beta = lambda_beta,
      lambda_gamma = 0,
      lambda_delta = 0,
      lambda_tau = 0,
      w_beta = w[j],
      w_gamma = 1,
      w_delta = 1,
      w_tau = 1,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      intercept = intercept
    )
    
    beta_old_value <- beta_hat[j]
    # Use the minimizer
    beta_hat[j] <- minimizer_Q_bern_beta(
      X = X_tilde,
      Z = Z_tilde,
      t = T_tilde,
      y = y,
      C = C,
      lambda = lambda_beta,
      beta_old = beta_old_value,
      weight = w[j],
      scaled = TRUE
    ) # Intercept is in C
    
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
    beta_all <- matrix(beta_all, ncol = 1)
    eta <- X %*% beta_all + intercept
    
    Q_new <- Q_bern(
      X = X,
      y = y,
      beta = beta_hat,
      gamma_vec = gamma_hat,
      delta_vec = delta_hat,
      tau_vec = tau_hat,
      lambda_beta = lambda_beta,
      lambda_gamma = 0,
      lambda_delta = 0,
      lambda_tau = 0,
      w_beta = w[j],
      w_gamma = 1,
      w_delta = 1,
      w_tau = 1,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      intercept = intercept
    )
    
    if (Q_new - Q_old > abs(Q_old) * 1e-2) {
      print("There might be numerical instability in update beta which was taken care of.")
    }
    # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
    if (Q_new - Q_old >= 0) {
      beta_hat[j] <- beta_old_value
    }
  }
  
  #k
  for (k in range3) {
    discard_from_c_main <- c(k)
    mains <- matrix(X_main[, k], nrow = length(eta))
    
    discard_from_c_2way <- c()
    # Contribution of 2-ways without beta
    two_ways <- 0
    for (il in c(range1, range2, range4)) {
      small_idx = min(il, k)
      big_idx = max(il, k)
      
      discard_from_c_2way <- c(
        discard_from_c_2way,
        get_position_vec_from_theta_matrix4(
          c(small_idx, big_idx),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )
      )
      two_ways <- two_ways + X_2way[, get_position_vec_from_theta_matrix4(
        c(small_idx, big_idx),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[il] *
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(small_idx, big_idx),
          l1 = l1,
          l2 = l2 ,
          l3 = l3,
          l4 = l4
        )]
    }
    
    discard_from_c_3way <- c()
    # Contribution of 3-ways without beta
    three_ways = 0
    for (i in c(range1, range2)) {
      for (l in  range4)
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[l]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
    }
    for (i in range1) {
      for (j in  range2)
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[j]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
    }
    
    discard_from_c_4way <- c()
    # Contribution of 4-ways without beta
    four_ways = 0
    for (i in range1) {
      for (j in range2)
      {
        for (l in range4)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, k),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, j),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[j] * beta_hat[l]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
    }
    
    X_tilde <- mains + two_ways
    Z_tilde <- three_ways
    T_tilde <- four_ways
    C <- intercept + X_main[, -k] %*% beta_hat[-k] + X_2way[, -discard_from_c_2way] %*%
      beta_2way[-discard_from_c_2way] +
      X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
      X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] # C is the offset term
    
    Q_old <- Q_bern(
      X = X,
      y = y,
      beta = beta_hat,
      gamma_vec = gamma_hat,
      delta_vec = delta_hat,
      tau_vec = tau_hat,
      lambda_beta = lambda_beta,
      lambda_gamma = 0,
      lambda_delta = 0,
      lambda_tau = 0,
      w_beta = w[k],
      w_gamma = 1,
      w_delta = 1,
      w_tau = 1,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      intercept = intercept
    )
    
    beta_old_value <- beta_hat[k]
    # Use the minimizer
    beta_hat[k] <- minimizer_Q_bern_beta(
      X = X_tilde,
      Z = Z_tilde,
      t = T_tilde,
      y = y,
      C = C,
      lambda = lambda_beta,
      beta_old = beta_old_value,
      weight = w[k],
      scaled = TRUE
    ) # Intercept is in C
    
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
    
    beta_all <- matrix(beta_all, ncol = 1)
    eta <- X %*% beta_all + intercept
    
    Q_new <- Q_bern(
      X = X,
      y = y,
      beta = beta_hat,
      gamma_vec = gamma_hat,
      delta_vec = delta_hat,
      tau_vec = tau_hat,
      lambda_beta = lambda_beta,
      lambda_gamma = 0,
      lambda_delta = 0,
      lambda_tau = 0,
      w_beta = w[k],
      w_gamma = 1,
      w_delta = 1,
      w_tau = 1,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      intercept = intercept
    )
    
    if (Q_new - Q_old > abs(Q_old) * 1e-2) {
      print("There might be numerical instability in update beta which was taken care of.")
    }
    # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
    if (Q_new - Q_old >= 0) {
      beta_hat[k] <- beta_old_value
    }
  }
  
  #l
  for (l in range4) {
    discard_from_c_main <- c(l)
    mains <- matrix(X_main[, l], nrow = length(eta))
    
    discard_from_c_2way <- c()
    # Contribution of 2-ways without beta
    two_ways <- 0
    for (i in c(range1, range2, range3)) {
      discard_from_c_2way <- c(
        discard_from_c_2way,
        get_position_vec_from_theta_matrix4(
          c(i, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )
      )
      two_ways <- two_ways + X_2way[, get_position_vec_from_theta_matrix4(
        c(i, l),
        l1 = l1,
        l2 = l2 ,
        l3 = l3,
        l4 = l4
      )] * beta_hat[i] *
        gamma_hat[get_position_vec_from_theta_matrix4(
          c(i, l),
          l1 = l1,
          l2 = l2 ,
          l3 = l3,
          l4 = l4
        )]
    }
    
    discard_from_c_3way <- c()
    # Contribution of 3-ways without beta
    three_ways = 0
    
    for (i in c(range1, range2)) {
      for (k in  range3)
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[k]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, k),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(k, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
    }
    for (i in range1) {
      for (j in  range2)
      {
        three_ways <- three_ways + X_3way[, psi_table_position_to_vector_index4(
          c(i, j, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] * (beta_hat[i] * beta_hat[j]) ^ 2 *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, j),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(i, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          gamma_hat[get_position_vec_from_theta_matrix4(
            c(j, l),
            l1 = l1,
            l2 = l2 ,
            l3 = l3,
            l4 = l4
          )] *
          delta_hat[psi_table_position_to_vector_index4(
            c(i, j, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        discard_from_c_3way <- c(
          discard_from_c_3way,
          psi_table_position_to_vector_index4(
            c(i, j, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )
        )
      }
    }
    
    discard_from_c_4way <- c()
    # Contribution of 4-ways without beta
    four_ways = 0
    for (i in range1) {
      for (j in range2)
      {
        for (k in range3)
        {
          four_ways <- four_ways + X_4way[, phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] * (beta_hat[i]) ^ 6 *
            (gamma_hat[get_position_vec_from_theta_matrix4(
              c(i, k),
              l1 = l1,
              l2 = l2 ,
              l3 = l3,
              l4 = l4
            )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, k),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(j, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(k, l),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )] *
              gamma_hat[get_position_vec_from_theta_matrix4(
                c(i, j),
                l1 = l1,
                l2 = l2 ,
                l3 = l3,
                l4 = l4
              )]) ^ 2 *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, k),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, j, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(i, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            delta_hat[psi_table_position_to_vector_index4(
              c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            tau_hat[phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
            (beta_hat[j] * beta_hat[k]) ^ 6
          discard_from_c_4way <- c(
            discard_from_c_4way,
            phi_table_position_to_vector_index4(
              c(i, j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )
          )
        }
      }
    }
    
    X_tilde <- mains + two_ways
    Z_tilde <- three_ways
    T_tilde <- four_ways
    
    C <- intercept + X_main[, -l] %*% beta_hat[-l] + X_2way[, -discard_from_c_2way] %*%
      beta_2way[-discard_from_c_2way] +
      X_3way[, -discard_from_c_3way] %*% beta_3way[-discard_from_c_3way] +
      X_4way[, -discard_from_c_4way] %*% beta_4way[-discard_from_c_4way] # C is the offset term
    
    Q_old <- Q_bern(
      X = X,
      y = y,
      beta = beta_hat,
      gamma_vec = gamma_hat,
      delta_vec = delta_hat,
      tau_vec = tau_hat,
      lambda_beta = lambda_beta,
      lambda_gamma = 0,
      lambda_delta = 0,
      lambda_tau = 0,
      w_beta = w[l],
      w_gamma = 1,
      w_delta = 1,
      w_tau = 1,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      intercept = intercept
    )
    
    beta_old_value <- beta_hat[l]
    # Use the minimizer
    beta_hat[l] <- minimizer_Q_bern_beta(
      X = X_tilde,
      Z = Z_tilde,
      t = T_tilde,
      y = y,
      C = C,
      lambda = lambda_beta,
      beta_old = beta_old_value,
      weight = w[l],
      scaled = TRUE
    ) # Intercept is in C
    
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
    beta_all <- matrix(beta_all, ncol = 1)
    eta <- X %*% beta_all + intercept
    
    Q_new <- Q_bern(
      X = X,
      y = y,
      beta = beta_hat,
      gamma_vec = gamma_hat,
      delta_vec = delta_hat,
      tau_vec = tau_hat,
      lambda_beta = lambda_beta,
      lambda_gamma = 0,
      lambda_delta = 0,
      lambda_tau = 0,
      w_beta = w[l],
      w_gamma = 1,
      w_delta = 1,
      w_tau = 1,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4,
      intercept = intercept
    )
    
    if (Q_new - Q_old > abs(Q_old) * 1e-2) {
      print("There might be numerical instability in update beta which was taken care of.")
    }
    # Keep the new estimator if only it is better than the old one. This ensures numerical stability.
    if (Q_new - Q_old >= 0) {
      beta_hat[l] <- beta_old_value
    }
  }
  
  print("Updated beta")
  
  return(beta_hat)
}
