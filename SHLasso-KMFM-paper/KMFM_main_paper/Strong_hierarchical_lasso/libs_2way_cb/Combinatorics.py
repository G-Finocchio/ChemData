import numpy as np


def get_range2(x, l1, l2):
    # Function to get the range of x (position of levels of factor corresponding to x)
    # returns the indices from the range of x
    
    assert x < l1 + l2 , "x should be in correct range"
    
    if x < l1:
        return list(range(0, l1))
    if x < l1 + l2:
        return list(range(l1, l1 + l2))



def get_ranges2(l1, l2):
    # Function to get ranges of main effects, two-way interactions, three-way, and four-way
    
    l_main = l1 + l2 
    l_theta = l1 * l2 
    range_main = list(range(0, l_main))
    range_theta = list(range(l_main, l_main + l_theta))
    
    return [range_main, range_theta]



#################### FUNCTIONS 2-way combinatorics ################################

# Function that transforms matrix Theta_hat into vec_theta
def get_theta_vec_2way2(Theta_hat, l1, l2):
    """
    Transform a Theta_hat matrix into vector form for 2-way interactions
    """
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    counter = 0
    vec_theta = np.zeros(l1 * l2)
    
    # Case 1: a with b, c, or d
    for i in range1:
        for j in range2  :
            vec_theta[counter] = (Theta_hat[i, j] + Theta_hat[j, i]) / 2
            counter += 1
    

    



# Get theta matrix from vector notation
def get_theta_from_theta_vec_2way2(vec_theta, l1, l2):
    """
    Get Theta_hat matrix from vec_theta for 2-way interactions
    """
    counter = 0
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))

    Theta_hat = np.zeros((l1 + l2  , l1 + l2 ))
    
    # Case 1: a with b, c
    for i in range1:
        for j in range2  :
            Theta_hat[i, j] = vec_theta[counter]
            Theta_hat[j, i] = vec_theta[counter]
            counter += 1
    
    expected_len = l1*l2   
    assert counter == expected_len, "Something wrong with counter"
    return Theta_hat


# Position in vector form from matrix indices for 2-way
def matrix_position_to_vector_index_2way2(position_tuple, l1, l2):
    """
    Returns position in vector form from a pair of matrix indices
    """
    x, y = position_tuple
    # 0-based indices
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    
    assert x < y, "x should be less than y"
    assert x in range1, "x should be in range 1"
    assert y in range2, "y should be in range 2"

    return x * l2   + (y - l1)



# Get pairwise products of main effects from different factors
def get_beta_vec_2way2(beta, l1, l2,  gamma, only_beta=False):
    """
    Returns vector of pairwise products in the same order as columns of X
    """
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    
    beta_vec2way = np.zeros(l1*l2)
    counter = 0
    
    for i in range1:
        for j in range2 :
            beta_vec2way[counter] = beta[i] * beta[j]
            counter += 1

    

    
    expected_len = l1*l2  
    assert counter == expected_len, "Counter mismatch in get_beta_vec_2way4"
    
    if not only_beta:
        beta_vec2way *= gamma
    
    return beta_vec2way


# Get list of positions in vector form from list of positions in matrix
def get_positions_2way2(ls_positions, l1, l2):
    all_positions = []
    for pos_tuple in ls_positions:
        pos_vec = matrix_position_to_vector_index_2way2(
            position_tuple=pos_tuple,
            l1=l1, l2=l2
        )
        all_positions.append(pos_vec)
    return all_positions









