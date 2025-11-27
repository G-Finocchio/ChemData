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








