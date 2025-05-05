source("Helpers.R")

######################## GENERAL FUNCTIONS ############################

assert <- function(condition, message) {
  if (!condition)
    stop(message)
}


# Function to get the range of x (position of levels of factor corresponding to x)
get_range4 <- function(x,
                       l1 = 21,
                       l2 = 14,
                       l3 = 2,
                       l4 = 3)
  #returns the indices from the range of x
{
  assert(x <= l1 + l2 + l3 + l4, "x should be in correct range")
  if (x <= l1)
  {
    return(c(1:l1))
  }
  if (x <= l1 + l2)
  {
    return(c((l1 + 1):(l1 + l2)))
  }
  if (x <= l1 + l2 + l3)
  {
    return(c((l1 + l2 + 1):(l1 + l2 + l3)))
  }
  return(c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4)))
}


#Function to get ranges of main effects , two-way interactions, three ways and four ways
get_ranges4 <- function(l1 = 21,
                        l2 = 14,
                        l3 = 2,
                        l4 = 3)
{
  l_main <- l1 + l2 + l3 + l4
  l_theta <- l1 * (l2 + l3 + l4) + l2 * (l3 + l4) + l3 * l4
  l_psi <- l1 * l2 * (l3 + l4) + l3 * l4 * (l1 + l2)
  l_phi <- l1 * l2 * l3 * l4
  
  range_main <- c(1:l_main)
  range_theta <- c((l_main + 1):(l_main + l_theta))
  range_psi <- c((l_main + l_theta + 1):(l_main + l_theta + l_psi))
  range_phi <- c((l_main + l_theta + l_psi + 1):(l_main + l_theta + l_psi + l_phi))
  
  return(list(range_main, range_theta, range_psi, range_phi))
}


#################### FUNCTIONS 2-way combinatorics ################################

# Function that transform matrix theta_hat in vec_theta
get_theta_vec_2way4 <- function(Theta_hat,
                                l1 = 21,
                                l2 = 14,
                                l3 = 2,
                                l4 = 3)
  
{
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  counter <- 0
  vec_theta <- array(0, l1 * (l2 + l3 + l4) + l2 * (l3 + l4) + l3 * l4)
  
  # Case 1 :a with b or c or d}
  for (i in range1)
  {
    for (j in c(range2, range3, range4))
    {
      counter <- counter + 1
      vec_theta[counter] <- (Theta_hat[i, j] + Theta_hat[j, i]) / 2
    }
  }
  
  # Case 2: b with c or d
  for (i in range2)
  {
    for (j in c(range3, range4))
    {
      counter <- counter + 1
      vec_theta[counter] <- (Theta_hat[i, j] + Theta_hat[j, i]) / 2
    }
  }
  
  # Case 3: c with d
  for (i in range3)
  {
    for (j in  range4)
    {
      counter <- counter + 1
      vec_theta[counter] <- (Theta_hat[i, j] + Theta_hat[j, i]) / 2
    }
  }
  
  return(vec_theta)
}


# Get theta in matrix notation from vector notation
get_theta_from_theta_vec_2way4 <- function(vec_theta,
                                           l1 = 21,
                                           l2 = 14,
                                           l3 = 2,
                                           l4 = 3)
  #get theta matrix from theta vec
{
  counter <- 1
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  Theta_hat <- matrix(0, nrow = l1 + l2 + l3 + l4, ncol = l1 + l2 + l3 + l4)
  
  # Case 1 :a with b or c or d
  for (i in range1)
  {
    for (j in c(range2, range3, range4))
    {
      Theta_hat[i, j] <- vec_theta[counter]
      Theta_hat[j, i] <- vec_theta[counter]
      counter <- counter + 1
    }
  }
  
  # Case 2: b with c or d
  for (i in range2)
  {
    for (j in c(range3, range4))
    {
      Theta_hat[i, j] <- vec_theta[counter]
      Theta_hat[j, i] <- vec_theta[counter]
      counter <- counter + 1
    }
  }
  
  # Case 3: c with d
  for (i in range3)
  {
    for (j in  range4)
    {
      Theta_hat[i, j] <- vec_theta[counter]
      Theta_hat[j, i] <- vec_theta[counter]
      counter <- counter + 1
    }
  }
  
  assert(counter == l1 * l2 + l2 * l3 + l3 * l1 + l4 * (l1 + l2 + l3) + 1,
         'smth wrong with counter')
  return(Theta_hat)
}



# Position in vector form to position in matrix form 2way
get_position_vec_from_theta_matrix4 <- function(position_tuple,
                                                l1 = 21,
                                                l2 = 14,
                                                l3 = 2,
                                                l4 = 3)
  # Takes into account / works only for possible combinations!
{
  x <- position_tuple[1]
  y <- position_tuple[2]
  
  range_x <- get_range4(
    x,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )
  range_y <- get_range4(
    y,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )
  
  
  assert(x <= l1 + l2 + l3 + l4, "x should be <=l1+l2")
  assert(x < y, "x<y")
  assert(y > l1, 'y should be >l1')
  
  
  if (all(range_x == c(1:l1)) == TRUE)
    #ab or ac or ad
  {
    position_vector <- (x - 1) * (l2 + l3 + l4) + (y - l1)
  }
  
  
  if (all (range_x == c((l1 + 1):(l1 + l2))) == TRUE)
    #bc or bd
    
  {
    position_vector <- l1 * (l2 + l3 + l4) + (x - l1 - 1) * (l3 + l4) + y - (l1 + l2)
  }
  
  
  if (all (range_x == c((l1 + l2 + 1):(l1 + l2 + l3))) == TRUE)
    #cd
  {
    position_vector <- l1 * (l2 + l3 + l4) + l2 * (l3 + l4) + (x - l1 - l2 - 1) * l4 + y - (l1 + l2 + l3)
  }
  return(position_vector)
}


# Get pairwise products of main effects from different factor. It has the same order as the columns of X
get_beta_vec_2way4 <- function(beta, l1, l2, l3, l4, gamma, only_beta = FALSE)
{
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  beta_vec2way <- array(0, dim = l1 * (l2 + l3 + l4) + l2 * l3 + l2 * l4 + l3 * l4)
  counter <- 1
  
  for (i in range1) {
    # ab ac ad
    for (j in  c(range2, range3, range4)) {
      beta_vec2way[counter] <- beta[i] * beta[j]
      counter <- counter + 1
    }
  }
  
  for (i in range2) {
    # bc bd
    for (j in  c(range3, range4)) {
      beta_vec2way[counter] <- beta[i] * beta[j]
      counter <- counter + 1
    }
  }
  
  for (i in range3) {
    # cd
    for (j in   range4) {
      beta_vec2way[counter] <- beta[i] * beta[j]
      counter <- counter + 1
    }
  }
  
  assert(counter == l1 * (l2 + l3 + l4) + l2 * l3 + l2 * l4 + l3 * l4 + 1)
  if (only_beta == FALSE)
  {
    beta_vec2way <- beta_vec2way * gamma
  }
  return(beta_vec2way)
}


# Gets ls of positions in vector form  from ls positions in matrix
get_positions_2way4 <- function(ls_positions,
                                l1 = 21,
                                l2 = 14,
                                l3 = 2,
                                l4 = 3) {
  all_positions <- c()
  for (tuple in ls_positions)
  {
    pos <- get_position_vec_from_theta_matrix4(
      position_tuple = tuple,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4
    )
    all_positions <- c(all_positions, pos)
  }
  return(all_positions)
}


#####################  FUNCTIONS 3-way combinatorics  ####################################

# Value of psi_vec at tablepsi[i,j,k]
psi_value_from_table_position <- function (table, i, j, k)
{
  return((table[i, j, k] + table[i, k, j] + table [j, i, k] + table[j, k, i] + table[k, i, j] + table[k, j, i]) / 6)
}

# Position in 3 dim table to vector index: works only for possible combinations
psi_table_position_to_vector_index4 <- function(position_tuple,
                                                l1 = 21,
                                                l2 = 14,
                                                l3 = 2,
                                                l4 = 3)
  # takes into account / works only for possible combinations
{
  x <- position_tuple[1]
  y <- position_tuple[2]
  z <- position_tuple[3]
  
  range_x <- get_range4(
    x,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )
  range_y <- get_range4(
    y,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )
  range_z <- get_range4(
    z,
    l1 = l1,
    l2 = l2,
    l3 = l3,
    l4 = l4
  )
  
  assert(x <= l1 + l2 + l3 + l4, "x should be <=l1+l2")
  assert(x < y, "x<y")
  assert(y > l1, 'y should be >l1')
  assert(y < z)
  
  if (all(range_x == c(1:l1)) == TRUE)
    # ab c/d, acd
  {
    if (all(range_y == c((l1 + 1):(l1 + l2))) == TRUE) {
      # ab c,d
      position_vector <- (x - 1) * (l2 * (l3 + l4) + l3 * l4) + (y - l1 -
                                                                   1) * (l3 + l4) + (z - l1 - l2)
    }
    if (all(range_y == c((l1 + l2 + 1):(l1 + l2 + l3))) == TRUE) {
      # acd
      position_vector <-    (x - 1) * (l2 * (l3 + l4) + l3 * l4) + l2 *
        (l3 + l4) + (y - l1 - l2 - 1) * l4 + (z - l1 - l2 - l3)
    }
  }
  
  if (all (range_x == c((l1 + 1):(l1 + l2))) == TRUE)
    # bcd
    
  {
    position_vector <- l1 * (l2 * (l3 + l4) + l3 * l4) + (x - l1 - 1) * (l3 *l4) +
      (y - (l1 + l2 + 1)) * l4 +  z - (l1 + l2 + l3)
  }
  return(position_vector)
}


# Get psi vec from 3dim table
get_psi_vec4 <- function(psi,
                         l1 = 21,
                         l2 = 14,
                         l3 = 2,
                         l4 = 3)
{
  assert(all(dim(psi) == l1 + l2 + l3 + l4), "Dimensions are not ok")
  
  psi_vec <- array(0, dim = c(l1 * l2 * l3 + l1 * l2 * l4 + l1 * l3 * l4 +
                                l2 * l3 * l4))
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  
  for (i in range1) {
    # abc
    for (j in range2) {
      for (k in range3) {
        psi_vec[psi_table_position_to_vector_index4(
          c(i, j, k),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- psi_value_from_table_position(psi, i, j, k)
      }
    }
  }
  
  for (i in range1) {
    # abd
    for (j in range2) {
      for (l in range4) {
        psi_vec[psi_table_position_to_vector_index4(
          c(i, j, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- psi_value_from_table_position(psi, i, j, l)
      }
    }
  }
  
  for (i in range1) {
    # abd
    for (k in range3) {
      for (l in range4) {
        psi_vec[psi_table_position_to_vector_index4(
          c(i, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- psi_value_from_table_position(psi, i, k, l)
      }
    }
  }
  
  for (j in range2) {
    # abd
    for (k in range3) {
      for (l in range4) {
        psi_vec[psi_table_position_to_vector_index4(
          c(j, k, l),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] <- psi_value_from_table_position(psi, j, k, l)
      }
    }
  }
  
  return(psi_vec)
}


# Get psi 3dim table from psi_vec
get_psi_from_psi_vec4 <- function(psi_vec,
                                  l1 = 21,
                                  l2 = 14,
                                  l3 = 2,
                                  l4 = 3)
{
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  counter <- 1
  psi <- array(0, dim = (c(l1 + l2 + l3 + l4, l1 + l2 + l3 + l4, l1 + l2 +
                             l3 + l4)))
  
  for (i in range1) {
    for (j in range2) {
      for (k in c(range3, range4)) {
        psi[i, j, k] <- psi_vec[counter]
        psi[i, k, j] <- psi_vec[counter]
        psi[k, i, j] <- psi_vec[counter]
        psi[k, j, i] <- psi_vec[counter]
        psi[j, i, k] <- psi_vec[counter]
        psi[j, k, i] <- psi_vec[counter]
        counter <- counter + 1
      }
    }
    
    for (j in range3) {
      for (k in range4) {
        psi[i, j, k] <- psi_vec[counter]
        psi[i, k, j] <- psi_vec[counter]
        psi[k, i, j] <- psi_vec[counter]
        psi[k, j, i] <- psi_vec[counter]
        psi[j, i, k] <- psi_vec[counter]
        psi[j, k, i] <- psi_vec[counter]
        counter <- counter + 1
      }
    }
    
  }
  
  for (i in range2) {
    for (j in range3) {
      for (k in range4) {
        psi[i, j, k] <- psi_vec[counter]
        psi[i, k, j] <- psi_vec[counter]
        psi[k, i, j] <- psi_vec[counter]
        psi[k, j, i] <- psi_vec[counter]
        psi[j, i, k] <- psi_vec[counter]
        psi[j, k, i] <- psi_vec[counter]
        counter <- counter + 1
      }
    }
  }
  
  assert(counter == 1 + l1 * l2 * (l3 + l4) + l3 * l4 * (l1 + l2))
  return(psi)
}


# Get products of two-way interactions in order pf columns of X for three-ways
get_beta_vec_3way4 <- function(beta_2way,
                               delta,
                               l1 = 21,
                               l2 = 14,
                               l3 = 2,
                               l4 = 3,
                               only_beta = FALSE)
{
  beta_vec3way <- array(0, dim = l1 * l2 * l3)
  counter <- 1
  range1 <- c(1:l1)
  range2 <- c ((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c ((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  
  #Iterate over possible positions
  for (i in range1) {
    # ab c/d
    for (j in range2) {
      for (k in c(range3, range4)) {
        beta_vec3way[counter] <- beta_2way[get_position_vec_from_theta_matrix4(
          position_tuple = c(i, j),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] *
          beta_2way[get_position_vec_from_theta_matrix4(
            position_tuple = c(i, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] *
          beta_2way[get_position_vec_from_theta_matrix4(
            position_tuple = c(j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        counter <- counter + 1
      }
    }
    
    for (j in range3) {
      for (k in  range4) {
        beta_vec3way[counter] <- beta_2way[get_position_vec_from_theta_matrix4(
          position_tuple = c(i, j),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] *
          beta_2way[get_position_vec_from_theta_matrix4(
            position_tuple = c(i, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] *
          beta_2way[get_position_vec_from_theta_matrix4(
            position_tuple = c(j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        counter <- counter + 1
      }
    }
  }
  
  for (i in range2) {
    #bcd
    for (j in range3) {
      for (k in  range4) {
        beta_vec3way[counter] <- beta_2way[get_position_vec_from_theta_matrix4(
          position_tuple = c(i, j),
          l1 = l1,
          l2 = l2,
          l3 = l3,
          l4 = l4
        )] *
          beta_2way[get_position_vec_from_theta_matrix4(
            position_tuple = c(i, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] *
          beta_2way[get_position_vec_from_theta_matrix4(
            position_tuple = c(j, k),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )]
        counter <- counter + 1
      }
    }
  }
  
  if (only_beta == FALSE)
  {
    beta_vec3way <- beta_vec3way * delta
  }
  assert(counter == l1 * l2 * (l3 + l4) + l3 * l4 * (l1 + l2) + 1)
  return(beta_vec3way)
}


# Get positions for list of tuples
get_positions_3way4 <- function(ls_positions,
                                l1 = 21,
                                l2 = 14,
                                l3 = 2,
                                l4 = 3)
{
  all_positions <- c()
  for (tuple in ls_positions)
  {
    pos <- psi_table_position_to_vector_index4(
      position_tuple = tuple,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4
    )
    all_positions <- c(all_positions, pos)
  }
  return(all_positions)
}


#####################  FUNCTIONS 4-way combinatorics  ####################################

phi_value_from_table_position <- function (table, i, j, k, l)
{
  return(
    (
      table[i, j, k, l] + table[i, j, l, k] + table [i, k, j, l] + table[i, k, l, j] + table[i, l, j, k] + table[i, l, k, j]
      + table[j, i, k, l] + table[j, i, l, k] + table [j, k, i, l] +
        table[j, k, l, i] + table[j, l, i, k] + table[j, l, k, i]
      + table[k, i, j, l] + table[k, i, l, j] + table [k, j, i, l] +
        table[k, j, l, i] + table[k, l, i, j] + table[k, l, j, i]
      + table[l, i, j, k] + table[l, i, k, j] + table [l, j, i, k] +
        table[l, j, k, i] + table[l, k, i, j] + table[l, k, j, i]
    ) / 24
  )
}


# Give position is phi_vec from position in 4d table for phi
phi_table_position_to_vector_index4 <- function(position_tuple,
                                                l1 = 21,
                                                l2 = 14,
                                                l3 = 2,
                                                l4 = 3)
  # Takes into account / works only for possible combinations!!!!
{
  x <- position_tuple[1]
  y <- position_tuple[2]
  z <- position_tuple[3]
  t <- position_tuple[4]
  
  assert(x <= l1, "x should be <=l1")
  assert(y <= l1 + l2, "y should be <=l1+l2")
  assert(z <= l1 + l2 + l3, "z should be <= than l1+l2+l3")
  assert(t <= l1 + l2 + l3 + l4, "t should be <= l1+l2+l3+l4")
  
  assert(x >= 1, "x should be >=1")
  assert(y > l1, "y should be >l1")
  assert(z > l1 + l2, 'z should be >l1+l2')
  assert(t > l1 + l2 + l3)
  
  position_phi <- (x - 1) * l2 * l3 * l4 + (y - l1 - 1) * l3 * l4 + (z -
                                                                       l1 - l2 - 1) * l4 + (t - l1 - l2 - l3)
  return(position_phi)
}


# phi vec from phi table
get_phi_vec4 <- function(phi,
                         l1 = 21,
                         l2 = 14,
                         l3 = 2,
                         l4 = 3)
{
  assert(all(dim(phi) == l1 + l2 + l3 + l4), "Dimensions are not ok")
  phi_vec <- array(0, dim = l1 * l2 * l3 * l4)
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  
  for (i in range1) {
    for (j in range2) {
      for (k in range3) {
        for (l in range4)
          phi_vec[phi_table_position_to_vector_index4(
            c(i, j, k, l),
            l1 = l1,
            l2 = l2,
            l3 = l3,
            l4 = l4
          )] <- phi_value_from_table_position(phi, i, j, k, l)
      }
    }
  }
  
  return(phi_vec)
}


# Get phi table from phi_vec
get_phi_from_phi_vec4 <- function(phi_vec,
                                  l1 = 21,
                                  l2 = 14,
                                  l3 = 2,
                                  l4 = 3)
{
  range1 <- c(1:l1)
  range2 <- c((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  counter <- 1
  phi <- array(0, dim = (c(
    l1 + l2 + l3 + l4, l1 + l2 + l3 + l4, l1 + l2 + l3 + l4, l1 + l2 + l3 +
      l4
  )))
  
  for (i in range1) {
    for (j in range2) {
      for (k in range3) {
        for (l in range4) {
          phi[i, j, k, l] <- phi_vec[counter]
          phi[i, j, l, k] <- phi_vec[counter]
          phi[i, k, j, l] <- phi_vec[counter]
          phi[i, k, l, j] <- phi_vec[counter]
          phi[i, l, j, k] <- phi_vec[counter]
          phi[i, l, k, j] <- phi_vec[counter]
          
          phi[j, i, k, l] <- phi_vec[counter]
          phi[j, i, l, k] <- phi_vec[counter]
          phi[j, k, i, l] <- phi_vec[counter]
          phi[j, k, l, i] <- phi_vec[counter]
          phi[j, l, i, k] <- phi_vec[counter]
          phi[j, l, k, i] <- phi_vec[counter]
          
          phi[k, i, j, l] <- phi_vec[counter]
          phi[k, i, l, j] <- phi_vec[counter]
          phi[k, j, i, l] <- phi_vec[counter]
          phi[k, j, l, i] <- phi_vec[counter]
          phi[k, l, i, j] <- phi_vec[counter]
          phi[k, l, j, i] <- phi_vec[counter]
          
          phi[l, i, j, k] <- phi_vec[counter]
          phi[l, i, k, j] <- phi_vec[counter]
          phi[l, j, i, k] <- phi_vec[counter]
          phi[l, j, k, i] <- phi_vec[counter]
          phi[l, k, i, j] <- phi_vec[counter]
          phi[l, k, j, i] <- phi_vec[counter]
          
          counter <- counter + 1
        }
      }
    }
  }
  assert(counter == 1 + l1 * l2 * l3 * l4, "counter error in phi")
  return(phi)
}


# Get products of three-way interactions in order pf columns of X for four-ways
get_beta_vec_4way4 <- function(beta_3way,
                               tau,
                               l1 = 21,
                               l2 = 14,
                               l3 = 2,
                               l4 = 3,
                               only_beta = FALSE)
{
  assert(length(beta_3way) == l1 * l2 * (l3 + l4) + l3 * l4 * (l1 + l2))
  beta_vec4way <- array(0, dim = l1 * l2 * l3 * l4)
  counter <- 1
  range1 <- c(1:l1)
  range2 <- c ((l1 + 1):(l1 + l2))
  range3 <- c((l1 + l2 + 1):(l1 + l2 + l3))
  range4 <- c ((l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4))
  
  # Iterate over possible positions
  for (i in range1) {
    #ab c/d
    for (j in range2) {
      for (k in range3) {
        {
          for (l in range4) {
            beta_vec4way[counter] <- beta_3way[psi_table_position_to_vector_index4(
              position_tuple = c(j, k, l),
              l1 = l1,
              l2 = l2,
              l3 = l3,
              l4 = l4
            )] *
              beta_3way[psi_table_position_to_vector_index4(
                position_tuple = c(i, j, k),
                l1 = l1,
                l2 = l2,
                l3 = l3,
                l4 = l4
              )] *
              beta_3way[psi_table_position_to_vector_index4(
                position_tuple = c(i, j, l),
                l1 = l1,
                l2 = l2,
                l3 = l3,
                l4 = l4
              )] *
              beta_3way[psi_table_position_to_vector_index4(
                position_tuple = c(i, k, l),
                l1 = l1,
                l2 = l2,
                l3 = l3,
                l4 = l4
              )]
            
            counter <- counter + 1
          }
        }
      }
    }
  }
  
  if (only_beta == FALSE)
  {
    beta_vec4way <- beta_vec4way * tau
  }
  assert(counter == l1 * l2 * l3 * l4 + 1)
  return(beta_vec4way)
}


# Get positions for list of tuples
get_positions_4way4 <- function(ls_positions,
                                l1 = 21,
                                l2 = 14,
                                l3 = 2,
                                l4 = 3)
{
  all_positions <- c()
  for (tuple in ls_positions)
  {
    pos <- phi_table_position_to_vector_index4(
      position_tuple = tuple,
      l1 = l1,
      l2 = l2,
      l3 = l3,
      l4 = l4
    )
    all_positions <- c(all_positions, pos)
  }
  return(all_positions)
}
