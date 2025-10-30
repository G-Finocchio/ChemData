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
