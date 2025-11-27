import numpy as np


def get_range3(x, l1=36, l2=3, l3=4):
    # Function to get the range of x (position of levels of factor corresponding to x)
    # returns the indices from the range of x
    
    assert x < l1 + l2 + l3 , "x should be in correct range"
    
    if x < l1:
        return list(range(0, l1))
    if x < l1 + l2:
        return list(range(l1, l1 + l2))
    if x < l1 + l2 + l3:
        return list(range(l1 + l2, l1 + l2 + l3))
    
    return list(range(l1 + l2 + l3, l1 + l2 + l3 ))


def get_ranges3(l1=36, l2=3, l3=4):
    # Function to get ranges of main effects, two-way interactions, three-way, and four-way
    
    l_main = l1 + l2 + l3
    l_theta = l1 * (l2 + l3 ) + l2 * l3
    l_psi = l1 * l2 * l3
    
    range_main = list(range(0, l_main))
    range_theta = list(range(l_main, l_main + l_theta))
    range_psi = list(range(l_main + l_theta, l_main + l_theta + l_psi))
    
    return [range_main, range_theta, range_psi]



#################### FUNCTIONS 2-way combinatorics ################################

# Function that transforms matrix Theta_hat into vec_theta
def get_theta_vec_2way3(Theta_hat, l1=36, l2=3, l3=4):
    """
    Transform a Theta_hat matrix into vector form for 2-way interactions
    """
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    range3 = list(range(l1 + l2, l1 + l2 + l3))
    counter = 0
    vec_theta = np.zeros(l1 * (l2 + l3 ) + l2 * l3)
    
    # Case 1: a with b, c, or d
    for i in range1:
        for j in range2 + range3 :
            vec_theta[counter] = (Theta_hat[i, j] + Theta_hat[j, i]) / 2
            counter += 1
    
    # Case 2: b with c or d
    for i in range2:
        for j in range3 :
            vec_theta[counter] = (Theta_hat[i, j] + Theta_hat[j, i]) / 2
            counter += 1
    



# Get theta matrix from vector notation
def get_theta_from_theta_vec_2way3(vec_theta, l1=36, l2=3, l3=4):
    """
    Get Theta_hat matrix from vec_theta for 2-way interactions
    """
    counter = 0
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    range3 = list(range(l1 + l2, l1 + l2 + l3))

    Theta_hat = np.zeros((l1 + l2 + l3 , l1 + l2 + l3 ))
    
    # Case 1: a with b, c
    for i in range1:
        for j in range2 + range3 :
            Theta_hat[i, j] = vec_theta[counter]
            Theta_hat[j, i] = vec_theta[counter]
            counter += 1
    
    # Case 2: b with c 
    for i in range2:
        for j in range3:
            Theta_hat[i, j] = vec_theta[counter]
            Theta_hat[j, i] = vec_theta[counter]
            counter += 1
    
    
    expected_len = l1*(l2 + l3) + l2*l3  
    assert counter == expected_len, "Something wrong with counter"
    return Theta_hat


# Position in vector form from matrix indices for 2-way
def matrix_position_to_vector_index_2way3(position_tuple, l1=36, l2=3, l3=4):
    """
    Returns position in vector form from a pair of matrix indices
    """
    x, y = position_tuple
    # 0-based indices
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    range3 = list(range(l1 + l2, l1 + l2 + l3))
    
    assert x < y, "x should be less than y"
    
    if x in range1:
        # ab, ac, ad
        return (x) * (l2 + l3 ) + (y - l1)
    elif x in range2:
        # bc, bd
        return l1 * (l2 + l3 ) + (x - l1) * l3  + (y - (l1 + l2))

    else:
        raise ValueError("Invalid x for 2-way interaction")


# Get pairwise products of main effects from different factors
def get_beta_vec_2way3(beta, l1, l2, l3,  gamma, only_beta=False):
    """
    Returns vector of pairwise products in the same order as columns of X
    """
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    range3 = list(range(l1 + l2, l1 + l2 + l3))
    
    beta_vec2way = np.zeros(l1*(l2 + l3) + l2*l3 )
    counter = 0
    
    for i in range1:
        for j in range2 + range3 :
            beta_vec2way[counter] = beta[i] * beta[j]
            counter += 1
    
    for i in range2:
        for j in range3 :
            beta_vec2way[counter] = beta[i] * beta[j]
            counter += 1
    

    
    expected_len = l1*(l2 + l3 ) + l2*l3 
    assert counter == expected_len, "Counter mismatch in get_beta_vec_2way4"
    
    if not only_beta:
        beta_vec2way *= gamma
    
    return beta_vec2way


# Get list of positions in vector form from list of positions in matrix
def get_positions_2way3(ls_positions, l1=36, l2=3, l3=4):
    all_positions = []
    for pos_tuple in ls_positions:
        pos_vec = matrix_position_to_vector_index_2way3(
            position_tuple=pos_tuple,
            l1=l1, l2=l2, l3=l3
        )
        all_positions.append(pos_vec)
    return all_positions




# ------------------- FUNCTIONS 3-way combinatorics -------------------

# Value of psi_vec at table[i,j,k] (0-based indexing)
def psi_value_from_table_position(table, i, j, k):
    return (table[i, j, k] + table[i, k, j] + table[j, i, k] +
            table[j, k, i] + table[k, i, j] + table[k, j, i]) / 6



# Map 3D table position to vector index (0-based)
def table_position_to_vector_index3(position_tuple, l1=36, l2=3, l3=4):
    x, y, z = position_tuple
    assert x < l1 + l2 + l3 , "x should be in correct range"
    assert x<y, "x should be smaller than y"
    assert y>=l1, "y should be >=l1"
    assert y<z, "y should be smaller than z"
    # ranges for reference
    range_x = get_range3(x, l1, l2, l3)
    range_y = get_range3(y, l1, l2, l3)
    range_z = get_range3(z, l1, l2, l3)


    return x*l2*l3 + (y-l1)*(l3) + z-(l1+l2)
    

# Get psi vector from 3D table (0-based)
def get_psi_vec3(psi, l1=36, l2=3, l3=4):
    vec_len = l1*l2*l3 
    psi_vec = np.zeros(vec_len, dtype=float)
    assert psi.shape == (l1+l2+l3, l1+l2+l3, l1+l2+l3), "Dimensions are not ok or psi is not 3D"

    range1 = np.arange(0, l1)
    range2 = np.arange(l1, l1+l2)
    range3 = np.arange(l1+l2, l1+l2+l3)

    # Fill psi_vec with values from 3D table
    for i in range1:
        for j in range2:
            for k in range3:
                idx = table_position_to_vector_index3((i, j, k), l1, l2, l3)
                psi_vec[idx] = psi_value_from_table_position(psi, i, j, k)

    
    return psi_vec



def get_psi_from_psi_vec3(psi_vec, l1=36, l2=3, l3=4):
    range1 = np.arange(l1)
    range2 = np.arange(l1, l1 + l2)
    range3 = np.arange(l1 + l2, l1 + l2 + l3)

    psi = np.zeros((l1 + l2 + l3 , l1 + l2 + l3 , l1 + l2 + l3 ))
    counter = 0

    for i in range1:
        for j in range2:
            for k in range3:
                for a, b, c in [(i,j,k),(i,k,j),(k,i,j),(k,j,i),(j,i,k),(j,k,i)]:
                    psi[a, b, c] = psi_vec[counter]
                counter += 1


    assert counter == l1 * l2 * l3 , "Counter mismatch"
    return psi


def get_beta_vec_3way3(beta_2way, delta, l1=36, l2=3, l3=4, only_beta=False):
    beta_vec3way = np.zeros(l1 * l2 * l3 )
    counter = 0
    range1 = np.arange(l1)
    range2 = np.arange(l1, l1 + l2)
    range3 = np.arange(l1 + l2, l1 + l2 + l3)


    for i in range1:
        for j in range2:
            for k in range3:
                beta_vec3way[counter] = (
                    beta_2way[matrix_position_to_vector_index_2way3((i,j), l1,l2,l3)] *
                    beta_2way[matrix_position_to_vector_index_2way3((i,k), l1,l2,l3)] *
                    beta_2way[matrix_position_to_vector_index_2way3((j,k), l1,l2,l3)]
                )
                counter += 1


    if not only_beta:
        beta_vec3way *= delta

    assert counter == l1 * l2 * l3  , "Counter mismatch"
    return beta_vec3way


def get_positions_3way3(ls_positions, l1=36, l2=3, l3=4):
    all_positions = []
    for position_tuple in ls_positions:
        pos = table_position_to_vector_index3(position_tuple, l1,l2,l3)
        all_positions.append(pos)
    return all_positions





