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


import numpy as np


def get_range4(x, l1=21, l2=14, l3=2, l4=3):
    # Function to get the range of x (position of levels of factor corresponding to x)
    # returns the indices from the range of x
    
    assert x < l1 + l2 + l3 + l4, "x should be in correct range"
    
    if x < l1:
        return list(range(0, l1))
    if x < l1 + l2:
        return list(range(l1, l1 + l2))
    if x < l1 + l2 + l3:
        return list(range(l1 + l2, l1 + l2 + l3))
    
    return list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))


def get_ranges4(l1=21, l2=14, l3=2, l4=3):
    # Function to get ranges of main effects, two-way interactions, three-way, and four-way
    
    l_main = l1 + l2 + l3 + l4
    l_theta = l1 * (l2 + l3 + l4) + l2 * (l3 + l4) + l3 * l4
    l_psi = l1 * l2 * (l3 + l4) + l3 * l4 * (l1 + l2)
    l_phi = l1 * l2 * l3 * l4
    
    range_main = list(range(0, l_main))
    range_theta = list(range(l_main, l_main + l_theta))
    range_psi = list(range(l_main + l_theta, l_main + l_theta + l_psi))
    range_phi = list(range(l_main + l_theta + l_psi, l_main + l_theta + l_psi + l_phi))
    
    return [range_main, range_theta, range_psi, range_phi]





#################### FUNCTIONS 2-way combinatorics ################################

# Function that transforms matrix Theta_hat into vec_theta
def get_theta_vec_2way4(Theta_hat, l1=21, l2=14, l3=2, l4=3):
    """
    Transform a Theta_hat matrix into vector form for 2-way interactions
    """
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    range3 = list(range(l1 + l2, l1 + l2 + l3))
    range4 = list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))
    counter = 0
    vec_theta = np.zeros(l1 * (l2 + l3 + l4) + l2 * (l3 + l4) + l3 * l4)
    
    # Case 1: a with b, c, or d
    for i in range1:
        for j in range2 + range3 + range4:
            vec_theta[counter] = (Theta_hat[i, j] + Theta_hat[j, i]) / 2
            counter += 1
    
    # Case 2: b with c or d
    for i in range2:
        for j in range3 + range4:
            vec_theta[counter] = (Theta_hat[i, j] + Theta_hat[j, i]) / 2
            counter += 1
    
    # Case 3: c with d
    for i in range3:
        for j in range4:
            vec_theta[counter] = (Theta_hat[i, j] + Theta_hat[j, i]) / 2
            counter += 1
    
    return vec_theta


# Get theta matrix from vector notation
def get_theta_from_theta_vec_2way4(vec_theta, l1=21, l2=14, l3=2, l4=3):
    """
    Get Theta_hat matrix from vec_theta for 2-way interactions
    """
    counter = 0
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    range3 = list(range(l1 + l2, l1 + l2 + l3))
    range4 = list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))
    Theta_hat = np.zeros((l1 + l2 + l3 + l4, l1 + l2 + l3 + l4))
    
    # Case 1: a with b, c, or d
    for i in range1:
        for j in range2 + range3 + range4:
            Theta_hat[i, j] = vec_theta[counter]
            Theta_hat[j, i] = vec_theta[counter]
            counter += 1
    
    # Case 2: b with c or d
    for i in range2:
        for j in range3 + range4:
            Theta_hat[i, j] = vec_theta[counter]
            Theta_hat[j, i] = vec_theta[counter]
            counter += 1
    
    # Case 3: c with d
    for i in range3:
        for j in range4:
            Theta_hat[i, j] = vec_theta[counter]
            Theta_hat[j, i] = vec_theta[counter]
            counter += 1
    
    expected_len = l1*(l2 + l3 + l4) + l2*(l3 + l4) + l3*l4
    assert counter == expected_len, "Something wrong with counter"
    return Theta_hat


# Position in vector form from matrix indices for 2-way
def get_position_vec_from_theta_matrix4(position_tuple, l1=21, l2=14, l3=2, l4=3):
    """
    Returns position in vector form from a pair of matrix indices
    """
    x, y = position_tuple
    # 0-based indices
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    range3 = list(range(l1 + l2, l1 + l2 + l3))
    range4 = list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))
    
    assert x < y, "x should be less than y"
    
    if x in range1:
        # ab, ac, ad
        return (x) * (l2 + l3 + l4) + (y - l1)
    elif x in range2:
        # bc, bd
        return l1 * (l2 + l3 + l4) + (x - l1) * (l3 + l4) + (y - (l1 + l2))
    elif x in range3:
        # cd
        return l1 * (l2 + l3 + l4) + l2 * (l3 + l4) + (x - (l1 + l2)) * l4 + (y - (l1 + l2 + l3))
    else:
        raise ValueError("Invalid x for 2-way interaction")


# Get pairwise products of main effects from different factors
def get_beta_vec_2way4(beta, l1, l2, l3, l4, gamma, only_beta=False):
    """
    Returns vector of pairwise products in the same order as columns of X
    """
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    range3 = list(range(l1 + l2, l1 + l2 + l3))
    range4 = list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))
    
    beta_vec2way = np.zeros(l1*(l2 + l3 + l4) + l2*l3 + l2*l4 + l3*l4)
    counter = 0
    
    for i in range1:
        for j in range2 + range3 + range4:
            beta_vec2way[counter] = beta[i] * beta[j]
            counter += 1
    
    for i in range2:
        for j in range3 + range4:
            beta_vec2way[counter] = beta[i] * beta[j]
            counter += 1
    
    for i in range3:
        for j in range4:
            beta_vec2way[counter] = beta[i] * beta[j]
            counter += 1
    
    expected_len = l1*(l2 + l3 + l4) + l2*l3 + l2*l4 + l3*l4
    assert counter == expected_len, "Counter mismatch in get_beta_vec_2way4"
    
    if not only_beta:
        beta_vec2way *= gamma
    
    return beta_vec2way


# Get list of positions in vector form from list of positions in matrix
def get_positions_2way4(ls_positions, l1=21, l2=14, l3=2, l4=3):
    all_positions = []
    for pos_tuple in ls_positions:
        pos_vec = get_position_vec_from_theta_matrix4(
            position_tuple=pos_tuple,
            l1=l1, l2=l2, l3=l3, l4=l4
        )
        all_positions.append(pos_vec)
    return all_positions



import numpy as np

# ------------------- FUNCTIONS 3-way combinatorics -------------------

# Value of psi_vec at table[i,j,k] (0-based indexing)
def psi_value_from_table_position(table, i, j, k):
    return (table[i, j, k] + table[i, k, j] + table[j, i, k] +
            table[j, k, i] + table[k, i, j] + table[k, j, i]) / 6



# Map 3D table position to vector index (0-based)
def psi_table_position_to_vector_index4(position_tuple, l1=21, l2=14, l3=2, l4=3):
    x, y, z = position_tuple
    assert x < l1 + l2 + l3 + l4, "x should be in correct range"
    assert x<y, "x should be smaller than y"
    assert y>=l1, "y should be >=l1"
    assert y<z, "y should be smaller than z"
    # ranges for reference
    range_x = get_range4(x, l1, l2, l3, l4)
    range_y = get_range4(y, l1, l2, l3, l4)
    range_z = get_range4(z, l1, l2, l3, l4)

    if x < l1:
        if np.array_equal(range_y, np.arange(l1, l1+l2)):
            return (x)*(l2*(l3+l4)+l3*l4) + (y-l1)*(l3+l4) + (z-l1-l2)
        elif np.array_equal(range_y, np.arange(l1+l2, l1+l2+l3)):
            return (x)*(l2*(l3+l4)+l3*l4) + l2*(l3+l4) + (y-l1-l2)*l4 + (z-l1-l2-l3)
    elif x < l1+l2:
        return l1*(l2*(l3+l4)+l3*l4) + (x-l1)*(l3*l4) + (y-(l1+l2))*l4 + (z-(l1+l2+l3))
    return None

# Get psi vector from 3D table (0-based)
def get_psi_vec4(psi, l1=21, l2=14, l3=2, l4=3):
    vec_len = l1*l2*(l3+l4) + l1*l3*l4 + l2*l3*l4
    psi_vec = np.zeros(vec_len, dtype=float)
    assert psi.shape == (l1+l2+l3+l4, l1+l2+l3+l4, l1+l2+l3+l4), "Dimensions are not ok or psi is not 3D"

    range1 = np.arange(0, l1)
    range2 = np.arange(l1, l1+l2)
    range3 = np.arange(l1+l2, l1+l2+l3)
    range4 = np.arange(l1+l2+l3, l1+l2+l3+l4)

    # Fill psi_vec with values from 3D table
    for i in range1:
        for j in range2:
            for k in np.concatenate([range3, range4]):
                idx = psi_table_position_to_vector_index4((i, j, k), l1, l2, l3, l4)
                psi_vec[idx] = psi_value_from_table_position(psi, i, j, k)
        for j in range3:
            for k in range4:
                idx = psi_table_position_to_vector_index4((i, j, k), l1, l2, l3, l4)
                psi_vec[idx] = psi_value_from_table_position(psi, i, j, k)

    for i in range2:
        for j in range3:
            for k in range4:
                idx = psi_table_position_to_vector_index4((i, j, k), l1, l2, l3, l4)
                psi_vec[idx] = psi_value_from_table_position(psi, i, j, k)
    
    return psi_vec



def get_psi_from_psi_vec4(psi_vec, l1=21, l2=14, l3=2, l4=3):
    range1 = np.arange(l1)
    range2 = np.arange(l1, l1 + l2)
    range3 = np.arange(l1 + l2, l1 + l2 + l3)
    range4 = np.arange(l1 + l2 + l3, l1 + l2 + l3 + l4)

    psi = np.zeros((l1 + l2 + l3 + l4, l1 + l2 + l3 + l4, l1 + l2 + l3 + l4))
    counter = 0

    for i in range1:
        for j in range2:
            for k in np.concatenate([range3, range4]):
                for a, b, c in [(i,j,k),(i,k,j),(k,i,j),(k,j,i),(j,i,k),(j,k,i)]:
                    psi[a, b, c] = psi_vec[counter]
                counter += 1

        for j in range3:
            for k in range4:
                for a, b, c in [(i,j,k),(i,k,j),(k,i,j),(k,j,i),(j,i,k),(j,k,i)]:
                    psi[a, b, c] = psi_vec[counter]
                counter += 1

    for i in range2:
        for j in range3:
            for k in range4:
                for a, b, c in [(i,j,k),(i,k,j),(k,i,j),(k,j,i),(j,i,k),(j,k,i)]:
                    psi[a, b, c] = psi_vec[counter]
                counter += 1

    assert counter == l1 * l2 * (l3 + l4) + l3 * l4 * (l1 + l2), "Counter mismatch"
    return psi


def get_beta_vec_3way4(beta_2way, delta, l1=21, l2=14, l3=2, l4=3, only_beta=False):
    beta_vec3way = np.zeros(l1 * l2 * (l3 + l4) + l3 * l4 * (l1 + l2))
    counter = 0
    range1 = np.arange(l1)
    range2 = np.arange(l1, l1 + l2)
    range3 = np.arange(l1 + l2, l1 + l2 + l3)
    range4 = np.arange(l1 + l2 + l3, l1 + l2 + l3 + l4)

    for i in range1:
        for j in range2:
            for k in np.concatenate([range3, range4]):
                beta_vec3way[counter] = (
                    beta_2way[get_position_vec_from_theta_matrix4((i,j), l1,l2,l3,l4)] *
                    beta_2way[get_position_vec_from_theta_matrix4((i,k), l1,l2,l3,l4)] *
                    beta_2way[get_position_vec_from_theta_matrix4((j,k), l1,l2,l3,l4)]
                )
                counter += 1

        for j in range3:
            for k in range4:
                beta_vec3way[counter] = (
                    beta_2way[get_position_vec_from_theta_matrix4((i,j), l1,l2,l3,l4)] *
                    beta_2way[get_position_vec_from_theta_matrix4((i,k), l1,l2,l3,l4)] *
                    beta_2way[get_position_vec_from_theta_matrix4((j,k), l1,l2,l3,l4)]
                )
                counter += 1

    for i in range2:
        for j in range3:
            for k in range4:
                beta_vec3way[counter] = (
                    beta_2way[get_position_vec_from_theta_matrix4((i,j), l1,l2,l3,l4)] *
                    beta_2way[get_position_vec_from_theta_matrix4((i,k), l1,l2,l3,l4)] *
                    beta_2way[get_position_vec_from_theta_matrix4((j,k), l1,l2,l3,l4)]
                )
                counter += 1

    if not only_beta:
        beta_vec3way *= delta

    assert counter == l1 * l2 * (l3 + l4) + l3 * l4 * (l1 + l2), "Counter mismatch"
    return beta_vec3way


def get_positions_3way4(ls_positions, l1=21, l2=14, l3=2, l4=3):
    all_positions = []
    for position_tuple in ls_positions:
        pos = psi_table_position_to_vector_index4(position_tuple, l1,l2,l3,l4)
        all_positions.append(pos)
    return all_positions





import numpy as np
from itertools import permutations


def phi_value_from_table_position(table, i, j, k, l):
    """Get average over all permutations of four indices in a 4D table."""
    indices = [i, j, k, l]
    sum_val = 0
    for perm in set(permutations(indices)):
        sum_val += table[perm]
    return sum_val / 24


def phi_table_position_to_vector_index4(position_tuple, l1=21, l2=14, l3=2, l4=3):
    """Convert a 4D table position to a vector index for phi_vec."""
    x, y, z, t = position_tuple

    assert 0 <= x < l1
    assert l1 <= y < l1 + l2
    assert l1 + l2 <= z < l1 + l2 + l3
    assert l1 + l2 + l3 <= t < l1 + l2 + l3 + l4

    return x * l2 * l3 * l4 + (y - l1) * l3 * l4 + (z - l1 - l2) * l4 + (t - l1 - l2 - l3)


def get_phi_vec4(phi, l1=21, l2=14, l3=2, l4=3):
    """Convert a 4D phi table to a 1D phi vector."""
    assert phi.shape == (l1 + l2 + l3 + l4,) * 4
    phi_vec = np.zeros(l1 * l2 * l3 * l4)
    range1 = range(l1)
    range2 = range(l1, l1 + l2)
    range3 = range(l1 + l2, l1 + l2 + l3)
    range4 = range(l1 + l2 + l3, l1 + l2 + l3 + l4)

    for i in range1:
        for j in range2:
            for k in range3:
                for l in range4:
                    idx = phi_table_position_to_vector_index4((i , j , k , l ), l1, l2, l3, l4)
                    phi_vec[idx] = phi_value_from_table_position(phi, i, j, k, l)
    return phi_vec


def get_phi_from_phi_vec4(phi_vec, l1=21, l2=14, l3=2, l4=3):
    """Convert a phi vector back to a 4D table."""
    phi = np.zeros((l1 + l2 + l3 + l4,) * 4)
    range1 = range(l1)
    range2 = range(l1, l1 + l2)
    range3 = range(l1 + l2, l1 + l2 + l3)
    range4 = range(l1 + l2 + l3, l1 + l2 + l3 + l4)

    counter = 0
    for i in range1:
        for j in range2:
            for k in range3:
                for l in range4:
                    val = phi_vec[counter]
                    indices = [(i, j, k, l), (i, j, l, k), (i, k, j, l), (i, k, l, j),
                               (i, l, j, k), (i, l, k, j), (j, i, k, l), (j, i, l, k),
                               (j, k, i, l), (j, k, l, i), (j, l, i, k), (j, l, k, i),
                               (k, i, j, l), (k, i, l, j), (k, j, i, l), (k, j, l, i),
                               (k, l, i, j), (k, l, j, i), (l, i, j, k), (l, i, k, j),
                               (l, j, i, k), (l, j, k, i), (l, k, i, j), (l, k, j, i)]
                    for idx in indices:
                        phi[idx] = val
                    counter += 1
    assert counter == l1 * l2 * l3 * l4, "counter error in phi"
    return phi


def get_beta_vec_4way4(beta_3way, tau, l1=21, l2=14, l3=2, l4=3, only_beta=False):
    """Get 4-way beta vector from 3-way beta vector."""
    beta_vec4way = np.zeros(l1 * l2 * l3 * l4)
    range1 = range(l1)
    range2 = range(l1, l1 + l2)
    range3 = range(l1 + l2, l1 + l2 + l3)
    range4 = range(l1 + l2 + l3, l1 + l2 + l3 + l4)

    counter = 0
    for i in range1:
        for j in range2:
            for k in range3:
                for l in range4:
                    beta_vec4way[counter] = (beta_3way[psi_table_position_to_vector_index4((j, k, l), l1, l2, l3, l4)]
                                             * beta_3way[psi_table_position_to_vector_index4((i, j, k), l1, l2, l3, l4)]
                                             * beta_3way[psi_table_position_to_vector_index4((i, j, l), l1, l2, l3, l4)]
                                             * beta_3way[psi_table_position_to_vector_index4((i, k, l), l1, l2, l3, l4)])
                    counter += 1
    if not only_beta:
        beta_vec4way *= tau
    assert counter == l1 * l2 * l3 * l4
    return beta_vec4way


def get_positions_4way4(ls_positions, l1=21, l2=14, l3=2, l4=3):
    """Get vector indices for list of 4D positions."""
    return [phi_table_position_to_vector_index4(pos, l1, l2, l3, l4) for pos in ls_positions]




import numpy as np




def get_ranges5(l1, l2, l3, l4, l5):
    # Function to get ranges of main effects, two-way interactions, three-way, and four-way
    
    l_main = l1 + l2 + l3 + l4 + l5
    l_theta = l1 * (l2 + l3 + l4 +l5) + l2 * (l3 + l4 +l5) + l3 * (l4 +l5) + l4*l5
    l_psi = l1 * l2 * (l3 + l4 +l5) + l1*l3*(l4+l5) +l1*l4*l5 + l2*l3*(l4+l5) +  l2*l4*l5 + l3*l4*l5
    l_phi = l1 * l2 * l3 * (l4+l5) + l1*l2*l4*l5 + l1*l3*l4*l5 + l2*l3*l4*l5
    l_omega=l1*l2*l3*l4*l5
    
    
    range_main = list(range(0, l_main))
    range_theta = list(range(l_main, l_main + l_theta))
    range_psi = list(range(l_main + l_theta, l_main + l_theta + l_psi))
    range_phi = list(range(l_main + l_theta + l_psi, l_main + l_theta + l_psi + l_phi))
    range_omega = list(range(l_main + l_theta + l_psi + l_phi, l_main + l_theta + l_psi + l_phi+ l_omega))
    
    return [range_main, range_theta, range_psi, range_phi, range_omega]





#################### FUNCTIONS 2-way combinatorics ################################




# Get pairwise products of main effects from different factors
def get_X_2way5(X, l1, l2, l3, l4, l5):
    """
    Returns vector of pairwise products in the same order as columns of X
    """
    range1 = list(range(l1))
    range2 = list(range(l1, l1 + l2))
    range3 = list(range(l1 + l2, l1 + l2 + l3))
    range4 = list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))
    range5 = list(range(l1 + l2 + l3 +l4, l1 + l2 + l3 + l4 +l5))
    X_2way = np.zeros( (X.shape[0], l1*(l2 + l3 + l4 +l5) + l2*(l3+l4+l5) + l3*(l4+l5) +l4*l5) )
    counter = 0
    
    for i in range1:
        for j in range2 + range3 + range4 +range5:
            X_2way[:,counter] = X[:,i]*X[:,j]
            counter += 1
    
    for i in range2:
        for j in range3 + range4 +range5:
            X_2way[:,counter] = X[:,i]*X[:,j]
            counter += 1
    
    for i in range3:
        for j in range4 +range5:
            X_2way[:,counter] = X[:,i]*X[:,j]
            counter += 1
            
    for i in range4:
        for j in range5:
            X_2way[:,counter] = X[:,i]*X[:,j]
            counter += 1
    
    assert counter == l1*(l2 + l3 + l4 +l5) + l2*(l3+l4+l5) + l3*(l4+l5) +l4*l5, "Counter mismatch in X_2way"
    
    
    return X_2way



def get_X_3way5(X, l1, l2, l3, l4, l5):
    """
    Returns matrix of 3-way interaction products of columns of X,
    grouped as:
      - i in r1, j in r2, k in r3+r4+r5
      - i in r1, j in r3, k in r4+r5
      - i in r1, j in r4, k in r5
      - i in r2, j in r3, k in r4+r5
      - i in r2, j in r4, k in r5
      - i in r3, j in r4, k in r5
    """

    # Define ranges
    r1 = list(range(l1))
    r2 = list(range(l1, l1 + l2))
    r3 = list(range(l1 + l2, l1 + l2 + l3))
    r4 = list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))
    r5 = list(range(l1 + l2 + l3 + l4, l1 + l2 + l3 + l4 + l5))

    n = X.shape[0]

    # total number of columns
    n_cols = (
        l1*l2*(l3+l4+l5) +
        l1*l3*(l4+l5) +
        l1*l4*l5 +
        l2*l3*(l4+l5) +
        l2*l4*l5 +
        l3*l4*l5
    )

    X_3way = np.zeros((n, n_cols))
    counter = 0

    # i in r1, j in r2, k in r3+r4+r5
    for i in r1:
        for j in r2:
            for k in r3 + r4 + r5:
                X_3way[:, counter] = X[:, i] * X[:, j] * X[:, k]
                counter += 1

    # i in r1, j in r3, k in r4+r5
    for i in r1:
        for j in r3:
            for k in r4 + r5:
                X_3way[:, counter] = X[:, i] * X[:, j] * X[:, k]
                counter += 1

    # i in r1, j in r4, k in r5
    for i in r1:
        for j in r4:
            for k in r5:
                X_3way[:, counter] = X[:, i] * X[:, j] * X[:, k]
                counter += 1

    # i in r2, j in r3, k in r4+r5
    for i in r2:
        for j in r3:
            for k in r4 + r5:
                X_3way[:, counter] = X[:, i] * X[:, j] * X[:, k]
                counter += 1

    # i in r2, j in r4, k in r5
    for i in r2:
        for j in r4:
            for k in r5:
                X_3way[:, counter] = X[:, i] * X[:, j] * X[:, k]
                counter += 1

    # i in r3, j in r4, k in r5
    for i in r3:
        for j in r4:
            for k in r5:
                X_3way[:, counter] = X[:, i] * X[:, j] * X[:, k]
                counter += 1

    # final sanity check
    assert counter == n_cols, f"Counter mismatch: got {counter}, expected {n_cols}"

    return X_3way


import numpy as np

def get_X_4way5(X, l1, l2, l3, l4, l5):
    """
    Returns matrix of 4-way interaction products of columns of X,
    grouped as:
      - i in r1, j in r2, k in r3, l in r4+r5
      - i in r1, j in r2, k in r4, l in r5
      - i in r1, j in r3, k in r4, l in r5
      - i in r2, j in r3, k in r4, l in r5
    """

    # Define ranges
    r1 = list(range(l1))
    r2 = list(range(l1, l1 + l2))
    r3 = list(range(l1 + l2, l1 + l2 + l3))
    r4 = list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))
    r5 = list(range(l1 + l2 + l3 + l4, l1 + l2 + l3 + l4 + l5))

    n = X.shape[0]

    # total number of 4-way columns
    n_cols = (
        l1*l2*l3*(l4+l5) +
        l1*l2*l4*l5 +
        l1*l3*l4*l5 +
        l2*l3*l4*l5
    )

    X_4way = np.zeros((n, n_cols))
    counter = 0

    # i in r1, j in r2, k in r3, l in r4+r5
    for i in r1:
        for j in r2:
            for k in r3:
                for l in r4 + r5:
                    X_4way[:, counter] = X[:, i] * X[:, j] * X[:, k] * X[:, l]
                    counter += 1

    # i in r1, j in r2, k in r4, l in r5
    for i in r1:
        for j in r2:
            for k in r4:
                for l in r5:
                    X_4way[:, counter] = X[:, i] * X[:, j] * X[:, k] * X[:, l]
                    counter += 1

    # i in r1, j in r3, k in r4, l in r5
    for i in r1:
        for j in r3:
            for k in r4:
                for l in r5:
                    X_4way[:, counter] = X[:, i] * X[:, j] * X[:, k] * X[:, l]
                    counter += 1

    # i in r2, j in r3, k in r4, l in r5
    for i in r2:
        for j in r3:
            for k in r4:
                for l in r5:
                    X_4way[:, counter] = X[:, i] * X[:, j] * X[:, k] * X[:, l]
                    counter += 1

    # sanity check
    assert counter == n_cols, f"Counter mismatch in get_X_4way5: got {counter}, expected {n_cols}"

    return X_4way

import numpy as np

def get_X_5way5(X, l1, l2, l3, l4, l5):
    """
    Returns matrix of 5-way interaction products of columns of X:
      i in r1, j in r2, k in r3, l in r4, m in r5
    """

    # Define ranges
    r1 = list(range(l1))
    r2 = list(range(l1, l1 + l2))
    r3 = list(range(l1 + l2, l1 + l2 + l3))
    r4 = list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))
    r5 = list(range(l1 + l2 + l3 + l4, l1 + l2 + l3 + l4 + l5))

    n = X.shape[0]

    # total number of 5-way columns
    n_cols = l1 * l2 * l3 * l4 * l5
    X_5way = np.zeros((n, n_cols))
    counter = 0

    for i in r1:
        for j in r2:
            for k in r3:
                for l in r4:
                    for m in r5:
                        X_5way[:, counter] = (
                            X[:, i] * X[:, j] * X[:, k] * X[:, l] * X[:, m]
                        )
                        counter += 1

    # sanity check
    assert counter == n_cols, f"Counter mismatch in get_X_5way5: got {counter}, expected {n_cols}"

    return X_5way


import numpy as np

def drop_zero_columns(X):
    """
    Given a matrix X, return a reduced version without all-zero columns
    and the indices of the kept columns.
    """
    # Find columns where at least one element is nonzero
    nonzero_mask = np.any(X != 0, axis=0)
    indices = np.where(nonzero_mask)[0]

    # Reduce the matrix
    X_reduced = X[:, indices]

    return X_reduced, indices


def get_2way_indices(l1, l2, l3, l4, l5):
    r1 = list(range(l1))
    r2 = list(range(l1, l1+l2))
    r3 = list(range(l1+l2, l1+l2+l3))
    r4 = list(range(l1+l2+l3, l1+l2+l3+l4))
    r5 = list(range(l1+l2+l3+l4, l1+l2+l3+l4+l5))
    
    # All 2-way combinations across groups
    two_way_pairs = []
    # i in r1, j in r2+r3+r4+r5
    for i in r1:
        for j in r2 + r3 + r4 + r5:
            two_way_pairs.append((i,j))
    for i in r2:
        for j in r3 + r4 + r5:
            two_way_pairs.append((i,j))
    for i in r3:
        for j in r4 + r5:
            two_way_pairs.append((i,j))
    for i in r4:
        for j in r5:
            two_way_pairs.append((i,j))
    
    return two_way_pairs  # list of tuples (col_i, col_j)


def get_mapping_3way_to_2way(l1,l2,l3,l4,l5):
    r1 = list(range(l1))
    r2 = list(range(l1, l1+l2))
    r3 = list(range(l1+l2, l1+l2+l3))
    r4 = list(range(l1+l2+l3, l1+l2+l3+l4))
    r5 = list(range(l1+l2+l3+l4, l1+l2+l3+l4+l5))

    two_way_pairs = get_2way_indices(l1,l2,l3,l4,l5)
    n_2way = len(two_way_pairs)
    
    three_way_cols = []

    # i ∈ r1
    for i in r1:
        for j in r2:
            for k in r3+r4+r5:
                three_way_cols.append((i,j,k))
        for j in r3:
            for k in r4+r5:
                three_way_cols.append((i,j,k))
        for j in r4:
            for k in r5:
                three_way_cols.append((i,j,k))

    # i ∈ r2
    for i in r2:
        for j in r3:
            for k in r4+r5:
                three_way_cols.append((i,j,k))
        for j in r4:
            for k in r5:
                three_way_cols.append((i,j,k))

    # i ∈ r3
    for i in r3:
        for j in r4:
            for k in r5:
                three_way_cols.append((i,j,k))

    n_3way = len(three_way_cols)
    mapping = np.zeros((n_3way, n_2way), dtype=int)

    # Fill mapping: 3-way involves 2-way if two indices match
    for tw_idx, (i,j) in enumerate(two_way_pairs):
        for th_idx, (a,b,c) in enumerate(three_way_cols):
            if (i,j) == (a,b) or (i,j) == (a,c) or (i,j) == (b,c):
                mapping[th_idx, tw_idx] = 1

    return mapping

def get_mapping_4way_to_3way(l1, l2, l3, l4, l5):
    # first generate all 3-way and 4-way column indices
    r1 = list(range(l1))
    r2 = list(range(l1, l1+l2))
    r3 = list(range(l1+l2, l1+l2+l3))
    r4 = list(range(l1+l2+l3, l1+l2+l3+l4))
    r5 = list(range(l1+l2+l3+l4, l1+l2+l3+l4+l5))

    # generate 3-way columns
    three_way_cols = []

    # i ∈ r1
    for i in r1:
        for j in r2:
            for k in r3+r4+r5:
                three_way_cols.append((i,j,k))
        for j in r3:
            for k in r4+r5:
                three_way_cols.append((i,j,k))
        for j in r4:
            for k in r5:
                three_way_cols.append((i,j,k))

    # i ∈ r2
    for i in r2:
        for j in r3:
            for k in r4+r5:
                three_way_cols.append((i,j,k))
        for j in r4:
            for k in r5:
                three_way_cols.append((i,j,k))

    # i ∈ r3
    for i in r3:
        for j in r4:
            for k in r5:
                three_way_cols.append((i,j,k))
    n_3way=len(three_way_cols)

    # generate 4-way columns
    four_way_cols = []
    for i in r1:
        for j in r2 :
            for k in r3:
                for l in  r4+r5 :
                    four_way_cols.append((i,j,k,l))
    for i in r1:
        for j in r2:
            for k in r4:
                for l in r5:          
                    four_way_cols.append((i,j,k,l))
    for i in r1:
        for j in r3:
            for k in r4:
                for l in r5:          
                    four_way_cols.append((i,j,k,l))
    for i in r2:
        for j in r3:
            for k in r4:
                for l in r5:  # nothing left
                    four_way_cols.append((i,j,k,l))   


    n_4way = len(four_way_cols)

    # initialize mapping
    mapping = np.zeros((n_4way, n_3way), dtype=int)

    # fill mapping: 4-way column involves 3-way col if any triple matches
    for th_idx, (a,b,c) in enumerate(three_way_cols):
        for fr_idx, (i,j,k,l) in enumerate(four_way_cols):
            if (a,b,c) == (i,j,k) or (a,b,c) == (i,j,l) or (a,b,c) == (i,k,l) or (a,b,c) == (j,k,l):
                mapping[fr_idx, th_idx] = 1

    return mapping














