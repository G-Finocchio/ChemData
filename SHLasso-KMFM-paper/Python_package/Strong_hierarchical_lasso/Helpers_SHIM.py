########### get_ranges functions that are used ######################
import numpy as np
import pandas as pd

def get_ranges2(l1, l2):
    # Function to get ranges of main effects, two-way interactions
    
    l_main = l1 + l2 
    l_theta = l1 * l2 
    range_main = list(range(0, l_main))
    range_theta = list(range(l_main, l_main + l_theta))
    
    return [range_main, range_theta]


def get_ranges3(l1=36, l2=3, l3=4):
    # Function to get ranges of main effects, two-way interactions, three-way
    
    l_main = l1 + l2 + l3
    l_theta = l1 * (l2 + l3 ) + l2 * l3
    l_psi = l1 * l2 * l3
    
    range_main = list(range(0, l_main))
    range_theta = list(range(l_main, l_main + l_theta))
    range_psi = list(range(l_main + l_theta, l_main + l_theta + l_psi))
    
    return [range_main, range_theta, range_psi]




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

def get_ranges5(l1, l2, l3, l4, l5):
    # Function to get ranges of main effects, two-way interactions, three-way, four-way and five-way
    
    l_main = l1 + l2 + l3 + l4 + l5
    l_theta = l1 * (l2 + l3 + l4 +l5) + l2 * (l3 + l4 +l5) + l3 * (l4 +l5) + l4*l5
    l_psi = l1 * l2 * (l3 + l4 +l5) + l1*l3*(l4+l5) +l1*l4*l5 + l2*l3*(l4+l5) +  l2*l4*l5 + l3*l4*l5
    l_phi = l1 * l2 * l3 * (l4+l5) + l1*l2*l4*l5 + l1*l3*l4*l5 + l2*l3*l4*l5
    l_omega=l1*l2*l3*l4*l5
    
    
    range_main = list(range(0, l_main))
    range_theta = list(range(l_main, l_main + l_theta))
    range_psi = list(range(l_main + l_theta, l_main + l_theta + l_psi))
    range_phi = list(range(l_main + l_theta + l_psi, l_main + l_theta + l_psi + l_phi))
    range_omega = list(range(l_main + l_theta + l_psi + l_phi, l_main + l_theta + l_psi + l_phi + l_omega))
    
    return [range_main, range_theta, range_psi, range_phi, range_omega]




    
    
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

    


def get_ranges_main( l1, l2, l3, l4, l5):
    
    r1=list(range(0, l1))
    r2=list(range(l1, l1 + l2))
    r3=list(range(l1 + l2, l1 + l2 + l3))
    r4 = list(range(l1 + l2 + l3, l1 + l2 + l3 + l4))
    r5= list(range(l1 + l2 + l3 + l4, l1 + l2 + l3 + l4 + l5))

    return r1,r2,r3,r4,r5


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

    
    
    
    
def recover_main(main):
    """
    Append one last "pseudo-level" for a main effect such that:
        sum(all levels) = 0.
    This enforces the effect-coding identifiability condition.
    """
    last = -np.sum(main)                      # compute last = -sum(others)
    return np.concatenate([main, [last]])     # append to main effects

    
    
def recover_tensor(tensor, shape):
    """
    Recover an N-way tensor enforcing sum-to-zero along all axes.

    Parameters
    ----------
    tensor : np.ndarray
        Flattened coefficients for N-way interaction.
    shape : tuple
        Original shape of the N-way interaction (without the last extra level).

    Returns
    -------
    full : np.ndarray
        Tensor with one extra level along each axis so that sums along
        each axis = 0.
    """
    ndim = len(shape)
    full_shape = tuple(s + 1 for s in shape)
    full = np.zeros(full_shape)

    # Fill the original tensor
    full[tuple(slice(0, s) for s in shape)] = tensor.reshape(shape, order="C")

    # Enforce sum-to-zero along each axis
    for axis in range(ndim):
        other_axes = [i for i in range(ndim) if i != axis]
        other_shape = [full_shape[i] for i in other_axes]

        for idx in np.ndindex(*other_shape):
            slc = [slice(None)] * ndim
            for k, ax in enumerate(other_axes):
                slc[ax] = idx[k]

            sub = full[tuple(slc)]
            # The current axis is always axis 0 in sub after slicing other axes
            full_idx = list(slc)
            full_idx[axis] = shape[axis]
            full[tuple(full_idx)] = -np.sum(sub[:shape[axis]], axis=0)

    return full




def recover_coefs_2way(beta_all, l1, l2):
    """
    Recover identifiable coefficients for a 2-way model.

    Enforces sum-to-zero (effect coding) on main and 2-way effects.

    Parameters
    ----------
    beta_all : np.ndarray
        Coefficients for a 2-way model before identifiability recovery.
    l1, l2 : int
        Number of modeled levels for factors 1 and 2 (before recovery).

    Returns
    -------
    dict
        {
            "main": np.ndarray,
            "theta": np.ndarray,
            "all": np.ndarray
        }
    """


    # ---------------- helper functions ----------------
    def recover_main(main):
        """Recover last main effect level by enforcing sum=0."""
        return np.concatenate([main, [-np.sum(main)]])




    # Get index ranges
    range_main, range_theta = get_ranges2(l1, l2)

    # --- Main effects ---
    main1 = recover_main(beta_all[range_main[:l1]])
    main2 = recover_main(beta_all[range_main[l1:l1 + l2]])
    main_full = np.concatenate([main1, main2])

    # --- 2-way interactions ---
    theta_flat = beta_all[range_theta]
    theta_full = recover_tensor(theta_flat, (l1, l2))

    # Flatten in C order (row-major)
    theta_flat_full = theta_full.ravel(order="C")

    return {
        "main": main_full,
        "theta": theta_flat_full,
        "all": np.concatenate([main_full, theta_flat_full]),
    }





def recover_coefs_3way(beta_all, l1, l2, l3):
    """
    Recover identifiable coefficients for a 3-way model.

    Enforces the sum-to-zero (effect coding) constraint on
    main, 2-way, and 3-way effects.

    Flattening order is row-major (C-style):
        (f1 outermost, f2 middle, f3 innermost).

    Parameters
    ----------
    beta_all : np.ndarray
        Coefficients for a 3-way model before identifiability recovery.
        Must include main, 2-way, and 3-way coefficients.
    l1, l2, l3 : int
        Number of modeled levels for factors 1, 2, and 3 (before recovery).

    Returns
    -------
    dict
        {
            "main": np.ndarray,
            "two_way": np.ndarray,
            "three_way": np.ndarray,
            "all": np.ndarray
        }
    """




    range_main, range_2way, range_3way = get_ranges3(l1, l2, l3)

    main1 = recover_main(beta_all[range_main[:l1]])
    main2 = recover_main(beta_all[range_main[l1:l1 + l2]])
    main3 = recover_main(beta_all[range_main[l1 + l2:l1 + l2 + l3]])

    # Concatenate all main effects together
    main_full = np.concatenate([main1, main2, main3])

    # ==========================================================
    #   TWO-WAY INTERACTIONS RECOVERY
    # ==========================================================
    # F1 x F2
    beta_2way = beta_all[range_2way]
    ls_pos_th12 = [
    matrix_position_to_vector_index_2way3((i, j), l1, l2, l3)
    for i in range(l1)
    for j in range(l1, l1+l2)
    ]
    # F1 x F3
    ls_pos_th13 = [
    matrix_position_to_vector_index_2way3((i, j), l1, l2, l3)
    for i in range(l1)
    for j in range(l1+l2, l1+l2+l3)
    ]

    # F2 x F3
    ls_pos_th23 = [
    matrix_position_to_vector_index_2way3((i, j), l1, l2, l3,)
    for i in range(l1,l1+l2)
    for j in range(l1+l2, l1+l2+l3)
    ]

    th_12=beta_2way[ls_pos_th12]
    th_13=beta_2way[ls_pos_th13]
    th_23=beta_2way[ls_pos_th23]

    th_12_full=recover_tensor(th_12, (l1,l2)).ravel(order="C")
    th_13_full=recover_tensor(th_13, (l1,l3)).ravel(order="C")
    th_23_full=recover_tensor(th_23, (l2,l3)).ravel(order="C")
    
    ### put Theta full in good order
    L1, L2, L3 =l1+1, l2+1, l3+1
    two_way_flat=np.zeros(L1*(L2+L3)+L2*L3)

    it=0
    for i in range(L1):
        for j in range(L2):
            pos_i=i
            pos_j=L1+j
            two_way_flat[matrix_position_to_vector_index_2way3((pos_i,pos_j), L1, L2, L3)]=th_12_full[it]
            it+=1
    assert it==L1*L2, "it should be L1L2"
    
    it=0
    for i in range(L1):
        for j in range(L3):
            pos_i=i
            pos_j=L1+L2+j
            two_way_flat[matrix_position_to_vector_index_2way3((pos_i,pos_j), L1, L2, L3)]=th_13_full[it]
            it+=1
    assert it==L1*L3, "it should be L1L3"
    

    it=0
    for i in range(L2):
        for j in range(L3):
            pos_i=L1 +i
            pos_j=L1+L2+j
            two_way_flat[matrix_position_to_vector_index_2way3((pos_i,pos_j), L1, L2, L3)]=th_23_full[it]
            it+=1
    assert it==L2*L3, "it should be L2L3"
    

    



    # ==========================================================
    #   THREE-WAY INTERACTIONS RECOVERY
    # ==========================================================
    arr3 = beta_all[range_3way]               # extract 3-way flat coefficients
    arr3_full = recover_tensor(arr3, (l1, l2, l3))  # recover full 3D tensor
    three_way_flat = arr3_full.ravel(order="C")      # flatten in C order

    # ==========================================================
    #   CONCATENATE EVERYTHING TO A SINGLE VECTOR
    # ==========================================================
    all_flat = np.concatenate([main_full, two_way_flat, three_way_flat])

    # ==========================================================
    # ✅  Return everything clearly labeled
    # ==========================================================
    return {
        "main": main_full,           # main effects (with recovered last levels)
        "two_way": two_way_flat,     # flattened 2-way interactions
        "three_way": three_way_flat, # flattened 3-way interactions
        "all": all_flat,             # everything concatenated
    }


def recover_coefs_4way(beta_all, l1, l2, l3, l4):
    """
    Recover identifiable coefficients for a 4-way model.

    Enforces sum-to-zero (effect coding) on main, 2-way, 3-way, and 4-way effects.

    Flattening order is C-style: f1 outermost, f2, f3, f4 innermost.

    Parameters
    ----------
    beta_all : np.ndarray
        Coefficients before identifiability recovery.
    l1,l2,l3,l4 : int
        Number of modeled levels for factors 1,2,3,4.

    Returns
    -------
    dict with keys "main", "two_way", "three_way", "four_way", "all"
    """

    # --- Get ranges (implement get_ranges4 separately) ---
    range_main, range_2way, range_3way, range_4way = get_ranges4(l1,l2,l3,l4)

    # --- Main effects ---
    main1 = recover_main(beta_all[range_main[:l1]])
    main2 = recover_main(beta_all[range_main[l1:l1+l2]])
    main3 = recover_main(beta_all[range_main[l1+l2:l1+l2+l3]])
    main4 = recover_main(beta_all[range_main[l1+l2+l3:l1+l2+l3+l4]])
    main_full = np.concatenate([main1, main2, main3, main4])

    # ==========================================================
    # 2-WAY INTERACTIONS
    # ==========================================================
    L1, L2, L3, L4 = l1+1, l2+1, l3+1, l4+1
    beta_2way = beta_all[range_2way]
    
    
    ls_pos_th12 = [
    get_position_vec_from_theta_matrix4((i, j), l1, l2, l3,l4)
    for i in range(l1)
    for j in range(l1, l1+l2)
    ]

    # F1 x F3
    ls_pos_th13 = [
    get_position_vec_from_theta_matrix4((i, j), l1, l2, l3,l4)
    for i in range(l1)
    for j in range(l1+l2, l1+l2+l3)
    ]
    
    # F1 x F4
    ls_pos_th14 = [
    get_position_vec_from_theta_matrix4((i, j), l1, l2, l3,l4)
    for i in range(l1)
    for j in range(l1+l2+l3, l1+l2+l3+l4)
    ]

    # F2 x F3
    ls_pos_th23 = [
    get_position_vec_from_theta_matrix4((i, j), l1, l2, l3,l4)
    for i in range(l1,l1+l2)
    for j in range(l1+l2, l1+l2+l3)
    ]
    
    # F2 x F4
    ls_pos_th24 = [
    get_position_vec_from_theta_matrix4((i, j), l1, l2, l3,l4)
    for i in range(l1,l1+l2)
    for j in range(l1+l2+l3, l1+l2+l3+l4)
    ]

    # F3 x F4
    ls_pos_th34 = [
    get_position_vec_from_theta_matrix4((i, j), l1, l2, l3,l4)
    for i in range(l1+l2,l1+l2+l3)
    for j in range(l1+l2+l3, l1+l2+l3+l4)
    ]
    
    th_12, th_13, th_14, th_23, th_24, th_34 = [
        beta_2way[idx] for idx in (ls_pos_th12, ls_pos_th13, ls_pos_th14, ls_pos_th23, ls_pos_th24, ls_pos_th34)]

    
    
    th_12_full=recover_tensor(th_12, (l1,l2)).ravel(order="C")
    th_13_full=recover_tensor(th_13, (l1,l3)).ravel(order="C")
    th_14_full=recover_tensor(th_14, (l1,l4)).ravel(order="C")
    th_23_full=recover_tensor(th_23, (l2,l3)).ravel(order="C")
    th_24_full=recover_tensor(th_24, (l2,l4)).ravel(order="C")
    th_34_full=recover_tensor(th_34, (l3,l4)).ravel(order="C")
    
    two_way_flat = np.zeros(L1*L2 + L1*L3 + L1*L4 + L2*L3 + L2*L4 + L3*L4)

    it = 0
    # F1 x F2
    for i in range(L1):
        for j in range(L2):
            pos_i = i
            pos_j = L1 + j
            two_way_flat[get_position_vec_from_theta_matrix4((pos_i,pos_j),L1,L2,L3,L4)] = th_12_full[it]
            it += 1
    assert it==L1*L2, "it should be L1L2"
            
    it = 0
    # F1 x F3
    for i in range(L1):
        for j in range(L3):
            pos_i = i
            pos_j = L1 + L2 + j
            two_way_flat[get_position_vec_from_theta_matrix4((pos_i,pos_j),L1,L2,L3,L4)] = th_13_full[it]
            it += 1
    assert it==L1*L3, "it should be L1L3"
    
    it=0
    # F1 x F4
    for i in range(L1):
        for j in range(L4):
            pos_i = i
            pos_j = L1 + L2 + L3 + j
            two_way_flat[get_position_vec_from_theta_matrix4((pos_i,pos_j),L1,L2,L3,L4)] = th_14_full[it]
            it += 1
    assert it==L1*L4, "it should be L1L4"
    
    it=0
    # F2 x F3
    for i in range(L2):
        for j in range(L3):
            pos_i = L1 + i
            pos_j = L1 + L2 + j
            two_way_flat[get_position_vec_from_theta_matrix4((pos_i,pos_j),L1,L2,L3,L4)] = th_23_full[it]
            it += 1
    assert it==L2*L3, "it should be L2L3"        
    
    it=0
    # F2 x F4
    for i in range(L2):
        for j in range(L4):
            pos_i = L1 + i
            pos_j = L1 + L2 + L3 + j
            two_way_flat[get_position_vec_from_theta_matrix4((pos_i,pos_j),L1,L2,L3,L4)] = th_24_full[it]
            it += 1
    assert it==L2*L4, "it should be L2L4"
    
    it=0        
    # F3 x F4
    for i in range(L3):
        for j in range(L4):
            pos_i = L1 + L2 + i
            pos_j = L1 + L2 + L3 + j
            two_way_flat[get_position_vec_from_theta_matrix4((pos_i,pos_j),L1,L2,L3,L4)] = th_34_full[it]
            it += 1
    assert it==L3*L4, "it should be L3L4"
    



    # ==========================================================
    # 3-WAY INTERACTIONS (empty for now)
    # ==========================================================
    L1, L2, L3, L4 = l1+1, l2+1, l3+1, l4+1
    beta_3way = beta_all[range_3way]
    
    
    ls_pos_th123 = [
    psi_table_position_to_vector_index4((i, j,k), l1, l2, l3,l4)
    for i in range(l1)
    for j in range(l1, l1+l2)
    for k in range(l1+l2,l1+l2+l3)
    ]
    
    ls_pos_th124 = [
    psi_table_position_to_vector_index4((i, j,k), l1, l2, l3,l4)
    for i in range(l1)
    for j in range(l1, l1+l2)
    for k in range(l1+l2+l3,l1+l2+l3+l4)
    ]
    
    ls_pos_th134 = [
    psi_table_position_to_vector_index4((i, j,k), l1, l2, l3,l4)
    for i in range(l1)
    for j in range(l1+l2, l1+l2+l3)
    for k in range(l1+l2+l3,l1+l2+l3+l4)
    ]
    
    ls_pos_th234 = [
    psi_table_position_to_vector_index4((i, j,k), l1, l2, l3,l4)
    for i in range(l1,l1+l2)
    for j in range(l1+l2, l1+l2+l3)
    for k in range(l1+l2+l3,l1+l2+l3+l4)
    ]
    
    th_123, th_124, th_134, th_234 = [beta_3way[idx] for idx in (ls_pos_th123, ls_pos_th124, ls_pos_th134, ls_pos_th234)]
    th_123_full=recover_tensor(th_123, (l1,l2,l3)).ravel(order="C")
    th_124_full=recover_tensor(th_124, (l1,l2,l4)).ravel(order="C")
    th_134_full=recover_tensor(th_134, (l1,l3,l4)).ravel(order="C")
    th_234_full=recover_tensor(th_234, (l2,l3,l4)).ravel(order="C")
    
    three_way_flat = np.zeros(L1*L2*L3 + L1*L2*L4 + L1*L3*L4 + L2*L3*L4)
    
    # F1 x F2 x F3
    it = 0
    for i in range(L1):
        for j in range(L2):
            for k in range(L3):
                pos_i = i
                pos_j = L1 + j
                pos_k = L1 + L2 + k
                three_way_flat[psi_table_position_to_vector_index4((pos_i, pos_j, pos_k), L1, L2, L3, L4)] = th_123_full[it]
                it += 1
    assert it == L1 * L2 * L3, "it should be L1L2L3"

    # F1 x F2 x F4
    it = 0
    for i in range(L1):
        for j in range(L2):
            for k in range(L4):
                pos_i = i
                pos_j = L1 + j
                pos_k = L1 + L2 + L3 + k
                three_way_flat[psi_table_position_to_vector_index4((pos_i, pos_j, pos_k), L1, L2, L3, L4)] = th_124_full[it]
                it += 1
    assert it == L1 * L2 * L4, "it should be L1L2L4"

    # F1 x F3 x F4
    it = 0
    for i in range(L1):
        for j in range(L3):
            for k in range(L4):
                pos_i = i
                pos_j = L1 + L2 + j
                pos_k = L1 + L2 + L3 + k
                three_way_flat[psi_table_position_to_vector_index4((pos_i, pos_j, pos_k), L1, L2, L3, L4)] = th_134_full[it]
                it += 1
    assert it == L1 * L3 * L4, "it should be L1L3L4"

    # F2 x F3 x F4
    it = 0
    for i in range(L2):
        for j in range(L3):
            for k in range(L4):
                pos_i = L1 + i
                pos_j = L1 + L2 + j
                pos_k = L1 + L2 + L3 + k
                three_way_flat[psi_table_position_to_vector_index4((pos_i, pos_j, pos_k), L1, L2, L3, L4)] = th_234_full[it]
                it += 1
    assert it == L2 * L3 * L4, "it should be L2L3L4"



    # ==========================================================
    # 4-WAY INTERACTIONS
    # ==========================================================
    beta_4way = beta_all[range_4way]
    four_way_full = recover_tensor(beta_4way, (l1,l2,l3,l4)).ravel(order="C")

    # ==========================================================
    # Concatenate everything
    # ==========================================================
    all_flat = np.concatenate([main_full, two_way_flat, three_way_flat, four_way_full])

    return {
        "main": main_full,
        "two_way": two_way_flat,
        "three_way": three_way_flat,
        "four_way": four_way_full,
        "all": all_flat
    }






def recover_coefs_5way(beta_all, l1, l2, l3,l4,l5):
    """
    Recover identifiable coefficients for a 5-way model.


    """


    range_main, range_2way, range_3way, range_4way, range_5way = get_ranges5(l1, l2, l3,l4,l5)


    main1 = recover_main(beta_all[range_main[:l1]])
    main2 = recover_main(beta_all[range_main[l1:l1 + l2]])
    main3 = recover_main(beta_all[range_main[l1 + l2 : l1 + l2 + l3]])
    main4 = recover_main(beta_all[range_main[l1 + l2 + l3 : l1 + l2 + l3 +l4]])
    main5 = recover_main(beta_all[range_main[l1 + l2 + l3 +l4 : l1 + l2 + l3 +l4 +l5]])

    # Concatenate all main effects together
    main_full = np.concatenate([main1, main2, main3,main4,main5])

    # ==========================================================
    #   TWO-WAY INTERACTIONS RECOVERY
    # ==========================================================
    # F1 x F2
    beta_2way = beta_all[range_2way]

    th12=np.zeros((l1,l2))
    th13=np.zeros((l1,l3))
    th14=np.zeros((l1,l4))
    th15=np.zeros((l1,l5))
    th23=np.zeros((l2,l3))
    th24=np.zeros((l2,l4))
    th25=np.zeros((l2,l5))
    th34=np.zeros((l3,l4))
    th35=np.zeros((l3,l5))
    th45=np.zeros((l4,l5))

    range1, range2, range3, range4, range5 = get_ranges_main(l1,l2,l3,l4,l5)
    # Assuming range1, range2, range3, range4, range5, and beta_2way are defined

    it = 0  # iterator over beta_2way

    # ---- i in range1 ----
    for i in range1:
        for j in range2 + range3 + range4 + range5:
            if j in range2:
                th12[i, j - l1] = beta_2way[it]
            elif j in range3:
                th13[i, j - l1 - l2] = beta_2way[it]
            elif j in range4:
                th14[i, j - l1 - l2 - l3] = beta_2way[it]
            elif j in range5:
                th15[i, j - l1 - l2 - l3 - l4] = beta_2way[it]
            it += 1

    # ---- i in range2 ----
    for i in range2:
        for j in range3 + range4 + range5:
            if j in range3:
                th23[i - l1, j - l1 - l2] = beta_2way[it]
            elif j in range4:
                th24[i - l1, j - l1 - l2 - l3] = beta_2way[it]
            elif j in range5:
                th25[i - l1, j - l1 - l2 - l3 - l4] = beta_2way[it]
            it += 1

    # ---- i in range3 ----
    for i in range3:
        for j in range4 + range5:
            if j in range4:
                th34[i - l1 - l2, j - l1 - l2 - l3] = beta_2way[it]
            elif j in range5:
                th35[i - l1 - l2, j - l1 - l2 - l3 - l4] = beta_2way[it]
            it += 1

    # ---- i in range4 ----
    for i in range4:
        for j in range5:
            th45[i - l1 - l2 - l3, j - l1 - l2 - l3 - l4] = beta_2way[it]
            it += 1



    th_12_full = recover_tensor(th12, (l1, l2))
    th_13_full = recover_tensor(th13, (l1, l3))
    th_14_full = recover_tensor(th14, (l1, l4))
    th_15_full = recover_tensor(th15, (l1, l5))
    th_23_full = recover_tensor(th23, (l2, l3))
    th_24_full = recover_tensor(th24, (l2, l4))
    th_25_full = recover_tensor(th25, (l2, l5))
    th_34_full = recover_tensor(th34, (l3, l4))
    th_35_full = recover_tensor(th35, (l3, l5))
    th_45_full = recover_tensor(th45, (l4, l5))


    
    ### put Theta full in good order
    L1, L2, L3, L4, L5 =l1+1, l2+1, l3+1, l4+1, l5+1
    two_way_flat=np.zeros(L1*(L2+L3+L4+L5)+L2*(L3+L4+L5) +L3*(L4+L5) +L4*L5)

    range1, range2, range3, range4, range5 = get_ranges_main(L1,L2,L3,L4,L5)

    it = 0  # iterator over two_way_flat

    # ---- i in range1 ----
    for i in range1:
        for j in range2 + range3 + range4 + range5:
            if j in range2:
                two_way_flat[it] = th_12_full[i, j - L1]
            elif j in range3:
                two_way_flat[it] = th_13_full[i, j - L1 - L2]
            elif j in range4:
                two_way_flat[it] = th_14_full[i, j - L1 - L2 - L3]
            elif j in range5:
                two_way_flat[it] = th_15_full[i, j - L1 - L2 - L3 - L4]
            it += 1

    assert it == L1 * (L2 + L3 + L4 + L5), "Iterator check failed after range1 group"

    # ---- i in range2 ----
    for i in range2:
        for j in range3 + range4 + range5:
            if j in range3:
                two_way_flat[it] = th_23_full[i - L1, j - L1 - L2]
            elif j in range4:
                two_way_flat[it] = th_24_full[i - L1, j - L1 - L2 - L3]
            elif j in range5:
                two_way_flat[it] = th_25_full[i - L1, j - L1 - L2 - L3 - L4]
            it += 1

    assert it == (
        L1 * (L2 + L3 + L4 + L5)
        + L2 * (L3 + L4 + L5)
    ), "Iterator check failed after range2 group"

    # ---- i in range3 ----
    for i in range3:
        for j in range4 + range5:
            if j in range4:
                two_way_flat[it] = th_34_full[i - L1 - L2, j - L1 - L2 - L3]
            elif j in range5:
                two_way_flat[it] = th_35_full[i - L1 - L2, j - L1 - L2 - L3 - L4]
            it += 1

    assert it == (
        L1 * (L2 + L3 + L4 + L5)
        + L2 * (L3 + L4 + L5)
        + L3 * (L4 + L5)
    ), "Iterator check failed after range3 group"

    # ---- i in range4 ----
    for i in range4:
        for j in range5:
            two_way_flat[it] = th_45_full[i - L1 - L2 - L3, j - L1 - L2 - L3 - L4]
            it += 1

    assert it == (
        L1 * (L2 + L3 + L4 + L5)
        + L2 * (L3 + L4 + L5)
        + L3 * (L4 + L5)
        + L4 * L5
    ), "Iterator check failed after range4 group"

    



    # ==========================================================
    #   THREE-WAY INTERACTIONS RECOVERY
    # ==========================================================

    beta_3way = beta_all[range_3way]
    

    th123=np.zeros((l1,l2,l3))
    th124=np.zeros((l1,l2,l4))
    th125=np.zeros((l1,l2,l5))
    th134=np.zeros((l1,l3,l4))
    th135=np.zeros((l1,l3,l5))
    th145=np.zeros((l1,l4,l5))
    th234=np.zeros((l2,l3,l4))
    th235=np.zeros((l2,l3,l5))
    th245=np.zeros((l2,l4,l5))
    th345=np.zeros((l3,l4,l5))

    range1, range2, range3, range4, range5 = get_ranges_main(l1, l2, l3, l4, l5)

    it = 0  # iterator over beta_3way

    # ---- i in range1 ----
    for i in range1:
        for j in range2 + range3 + range4 :
            for k in range(j + 1, l1 + l2 + l3 + l4 + l5):
                if j in range2 and k in range3:
                    th123[i, j - l1, k - l1 - l2] = beta_3way[it]
                    it += 1
                elif j in range2 and k in range4:
                    th124[i, j - l1, k - l1 - l2 - l3] = beta_3way[it]
                    it += 1
                elif j in range2 and k in range5:
                    th125[i, j - l1, k - l1 - l2 - l3 - l4] = beta_3way[it]
                    it += 1
                elif j in range3 and k in range4:
                    th134[i, j - l1 - l2, k - l1 - l2 - l3] = beta_3way[it]
                    it += 1
                elif j in range3 and k in range5:
                    th135[i, j - l1 - l2, k - l1 - l2 - l3 - l4] = beta_3way[it]
                    it += 1
                elif j in range4 and k in range5:
                    th145[i, j - l1 - l2 - l3, k - l1 - l2 - l3 - l4] = beta_3way[it]
                    it += 1

    assert it == (
        l1 * (l2 * (l3 + l4 + l5) + l3 * (l4 + l5) + l4 * l5)
    ), "Iterator check failed after range1 group"

    # ---- i in range2 ----
    for i in range2:
        for j in range3 + range4:
            for k in range(j + 1, l1 + l2 + l3 + l4 + l5):
                if j in range3 and k in range4:
                    th234[i - l1, j - l1 - l2, k - l1 - l2 - l3] = beta_3way[it]
                    it += 1
                elif j in range3 and k in range5:
                    th235[i - l1, j - l1 - l2, k - l1 - l2 - l3 - l4] = beta_3way[it]
                    it += 1
                elif j in range4 and k in range5:
                    th245[i - l1, j - l1 - l2 - l3, k - l1 - l2 - l3 - l4] = beta_3way[it]
                    it += 1

    assert it == (
        l1 * (l2 * (l3 + l4 + l5) + l3 * (l4 + l5) + l4 * l5)
        + l2 * (l3 * (l4 + l5) + l4 * l5)
    ), "Iterator check failed after range2 group"

    # ---- i in range3 ----
    for i in range3:
        for j in range4 + range5:
            for k in range(j + 1, l1 + l2 + l3 + l4 + l5):
                if j in range4 and k in range5:
                    th345[i - l1 - l2, j - l1 - l2 - l3, k - l1 - l2 - l3 - l4] = beta_3way[it]
                    it += 1

    assert it == (
        l1 * (l2 * (l3 + l4 + l5) + l3 * (l4 + l5) + l4 * l5)
        + l2 * (l3 * (l4 + l5) + l4 * l5)
        + l3 * (l4 * l5)
    ), "Iterator check failed after range3 group"


    th_123_full = recover_tensor(th123, (l1, l2, l3))
    th_124_full = recover_tensor(th124, (l1, l2, l4))
    th_125_full = recover_tensor(th125, (l1, l2, l5))
    th_134_full = recover_tensor(th134, (l1, l3, l4))
    th_135_full = recover_tensor(th135, (l1, l3, l5))
    th_145_full = recover_tensor(th145, (l1, l4, l5))
    th_234_full = recover_tensor(th234, (l2, l3, l4))
    th_235_full = recover_tensor(th235, (l2, l3, l5))
    th_245_full = recover_tensor(th245, (l2, l4, l5))
    th_345_full = recover_tensor(th345, (l3, l4, l5))


    # ---- put 3-way Theta full in good order ----
    L1, L2, L3, L4, L5 = l1 + 1, l2 + 1, l3 + 1, l4 + 1, l5 + 1

    three_way_flat = np.zeros(
        L1 * (L2 * (L3 + L4 + L5) + L3 * (L4 + L5) + L4 * L5)
        + L2 * (L3 * (L4 + L5) + L4 * L5)
        + L3 * (L4 * L5)
    )

    range1, range2, range3, range4, range5 = get_ranges_main(L1, L2, L3, L4, L5)

    it = 0  # iterator over three_way_flat

    # ---- i in range1 ----
    for i in range1:
        for j in range2 + range3 + range4:
            for k in range(j + 1, L1 + L2 + L3 + L4 + L5):
                if j in range2 and k in range3:
                    three_way_flat[it] = th_123_full[i, j - L1, k - L1 - L2]
                    it += 1
                elif j in range2 and k in range4:
                    three_way_flat[it] = th_124_full[i, j - L1, k - L1 - L2 - L3]
                    it += 1
                elif j in range2 and k in range5:
                    three_way_flat[it] = th_125_full[i, j - L1, k - L1 - L2 - L3 - L4]
                    it += 1
                elif j in range3 and k in range4:
                    three_way_flat[it] = th_134_full[i, j - L1 - L2, k - L1 - L2 - L3]
                    it += 1
                elif j in range3 and k in range5:
                    three_way_flat[it] = th_135_full[i, j - L1 - L2, k - L1 - L2 - L3 - L4]
                    it += 1
                elif j in range4 and k in range5:
                    three_way_flat[it] = th_145_full[i, j - L1 - L2 - L3, k - L1 - L2 - L3 - L4]
                    it += 1

    assert it == L1 * (L2 * (L3 + L4 + L5) + L3 * (L4 + L5) + L4 * L5), \
        "Iterator check failed after range1 group"

    # ---- i in range2 ----
    for i in range2:
        for j in range3 + range4:
            for k in range(j + 1, L1 + L2 + L3 + L4 + L5):
                if j in range3 and k in range4:
                    three_way_flat[it] = th_234_full[i - L1, j - L1 - L2, k - L1 - L2 - L3]
                    it += 1
                elif j in range3 and k in range5:
                    three_way_flat[it] = th_235_full[i - L1, j - L1 - L2, k - L1 - L2 - L3 - L4]
                    it += 1
                elif j in range4 and k in range5:
                    three_way_flat[it] = th_245_full[i - L1, j - L1 - L2 - L3, k - L1 - L2 - L3 - L4]
                    it += 1

    assert it == L1 * (L2 * (L3 + L4 + L5) + L3 * (L4 + L5) + L4 * L5) \
        + L2 * (L3 * (L4 + L5) + L4 * L5), "Iterator check failed after range2 group"

    # ---- i in range3 ----
    for i in range3:
        for j in range4:
            for k in range(j + 1, L1 + L2 + L3 + L4 + L5):
                if j in range4 and k in range5:
                    three_way_flat[it] = th_345_full[i - L1 - L2, j - L1 - L2 - L3, k - L1 - L2 - L3 - L4]
                    it += 1

    assert it == L1 * (L2 * (L3 + L4 + L5) + L3 * (L4 + L5) + L4 * L5) \
        + L2 * (L3 * (L4 + L5) + L4 * L5) \
        + L3 * (L4 * L5), "Iterator check failed after range3 group"
    



    # ==========================================================
    #   FOUR-WAY INTERACTIONS RECOVERY & FLATTEN
    # ==========================================================
    beta_4way = beta_all[range_4way]

    # ---- create 4-way tensors ----
    th1234 = np.zeros((l1, l2, l3, l4))
    th1235 = np.zeros((l1, l2, l3, l5))
    th1245 = np.zeros((l1, l2, l4, l5))
    th1345 = np.zeros((l1, l3, l4, l5))
    th2345 = np.zeros((l2, l3, l4, l5))

    # ---- get ranges ----
    range1, range2, range3, range4, range5 = get_ranges_main(l1, l2, l3, l4, l5)

    it = 0  # iterator over beta_4way

    # ---- i in range1 ----
    for i in range1:
        for j in range2 + range3 :
            for k in range(j + 1, l1 + l2 + l3 + l4):
                for m in range(k + 1, l1 + l2 + l3 + l4 + l5):
                    if j in range2 and k in range3 and m in range4:
                        th1234[i, j - l1, k - l1 - l2, m - l1 - l2 - l3] = beta_4way[it]
                        it += 1
                    elif j in range2 and k in range3 and m in range5:
                        th1235[i, j - l1, k - l1 - l2, m - l1 - l2 - l3 - l4] = beta_4way[it]
                        it += 1
                    elif j in range2 and k in range4 and m in range5:
                        th1245[i, j - l1, k - l1 - l2 - l3, m - l1 - l2 - l3 - l4] = beta_4way[it]
                        it += 1
                    elif j in range3 and k in range4 and m in range5:
                        th1345[i, j - l1 - l2, k - l1 - l2 - l3, m - l1 - l2 - l3 - l4] = beta_4way[it]
                        it += 1

    assert it == l1 * (l2 * l3 * l4 + l2 * l3 * l5 + l2 * l4 * l5 + l3 * l4 * l5), \
        "Iterator check failed after range1 group"

    # ---- i in range2 ----
    for i in range2:
        for j in range3 :
            for k in range(j + 1, l1 + l2 + l3 + l4 ):
                for m in range(k + 1, l1 + l2 + l3 + l4 + l5):
                    if j in range3 and k in range4 and m in range5:
                        th2345[i - l1, j - l1 - l2, k - l1 - l2 - l3, m - l1 - l2 - l3 - l4] = beta_4way[it]
                        it += 1

    assert it == l1 * (l2 * l3 * l4 + l2 * l3 * l5 + l2 * l4 * l5 + l3 * l4 * l5) \
                 + l2 * (l3 * l4 * l5), "Iterator check failed after range2 group"


    assert it == l1 * (l2 * l3 * l4 + l2 * l3 * l5 + l2 * l4 * l5 + l3 * l4 * l5) \
                 + l2 * (l3 * l4 * l5) , "Iterator check failed after range3 group"

    # ---- recover tensors to flat arrays ----
    th_1234_full = recover_tensor(th1234, (l1, l2, l3, l4))
    th_1235_full = recover_tensor(th1235, (l1, l2, l3, l5))
    th_1245_full = recover_tensor(th1245, (l1, l2, l4, l5))
    th_1345_full = recover_tensor(th1345, (l1, l3, l4, l5))
    th_2345_full = recover_tensor(th2345, (l2, l3, l4, l5))

    # ---- put 4-way Theta full in good order ----

    range1, range2, range3, range4, range5 = get_ranges_main(L1, L2, L3, L4, L5)
    four_way_flat = np.zeros(
        L1 * (L2 * L3 * L4 + L2 * L3 * L5 + L2 * L4 * L5 + L3 * L4 * L5)
        + L2 * (L3 * L4 * L5)
   
    )

    it = 0  # iterator over four_way_flat

    # ---- i in range1 ----
    for i in range1:
        for j in range2 + range3 :
            for k in range(j + 1, L1 + L2 + L3 + L4 ):
                for m in range(k + 1, L1 + L2 + L3 + L4 + L5):
                    if j in range2 and k in range3 and m in range4:
                        four_way_flat[it] = th_1234_full[i, j - L1, k - L1 - L2, m - L1 - L2 - L3]
                        it += 1
                    elif j in range2 and k in range3 and m in range5:
                        four_way_flat[it] = th_1235_full[i, j - L1, k - L1 - L2, m - L1 - L2 - L3 - L4]
                        it += 1
                    elif j in range2 and k in range4 and m in range5:
                        four_way_flat[it] = th_1245_full[i, j - L1, k - L1 - L2 - L3, m - L1 - L2 - L3 - L4]
                        it += 1
                    elif j in range3 and k in range4 and m in range5:
                        four_way_flat[it] = th_1345_full[i, j - L1 - L2, k - L1 - L2 - L3, m - L1 - L2 - L3 - L4]
                        it += 1

    assert it == L1 * (L2 * L3 * L4 + L2 * L3 * L5 + L2 * L4 * L5 + L3 * L4 * L5), \
        "Iterator check failed after range1 group"

    # ---- i in range2 ----
    for i in range2:
        for j in range3 :
            for k in range4:
                for m in range5:
                    if j in range3 and k in range4 and m in range5:
                        four_way_flat[it] = th_2345_full[i - L1, j - L1 - L2, k - L1 - L2 - L3, m - L1 - L2 - L3 - L4]
                        it += 1

    assert it == L1 * (L2 * L3 * L4 + L2 * L3 * L5 + L2 * L4 * L5 + L3 * L4 * L5) +L2*L3*L4*L5 , "Iterator check failed after range2 group"




    # ==========================================================
    # 6️⃣  CONCATENATE EVERYTHING TO A SINGLE VECTOR
    # ==========================================================
    all_flat = np.concatenate([main_full, two_way_flat, three_way_flat, four_way_flat])

    # ==========================================================
    # ✅  Return everything clearly labeled
    # ==========================================================
    return {
        "main": main_full,           # main effects (with recovered last levels)
        "two_way": two_way_flat,     # flattened 2-way interactions
        "three_way": three_way_flat, # flattened 3-way interactions
        "four_way": four_way_flat,
        "all": all_flat,             # everything concatenated
    }

    

def recover_coefs(beta_all, levels):
    """
    Wrapper to recover identifiable coefficients for 2-5 factor models.

    Parameters
    ----------
    beta_all : np.ndarray
        Flattened coefficients for the model (main, 2-way, etc.).
    levels : list of int
        Number of levels for each factor. Length determines model type (2-5 factors).

    Returns
    -------
    dict
        Output from the appropriate recover_coefs_*way function.
    """
    n_factors = len(levels)
    assert 2 <= n_factors <= 5, "Only 2-5 factors are supported."

    if n_factors == 2:
        return recover_coefs_2way(beta_all, levels[0], levels[1])
    elif n_factors == 3:
        return recover_coefs_3way(beta_all, levels[0], levels[1], levels[2])
    elif n_factors == 4:
        return recover_coefs_4way(beta_all, levels[0], levels[1], levels[2], levels[3])
    elif n_factors == 5:
        return recover_coefs_5way(beta_all, levels[0], levels[1], levels[2], levels[3], levels[4])




################################## GENERATE COEFS NAME #####################################



def generate_names_2factors(main_effects, n_levels):
    f1 = main_effects[:n_levels[0]]
    f2 = main_effects[n_levels[0]:n_levels[0]+n_levels[1]]

    names = f1 + f2
    # 2-way interactions
    for a in f1:
        for b in f2:
            names.append(f"{a}:{b}")
    return names


def generate_names_3factors(main_effects, n_levels):
    idx = 0
    f1 = main_effects[idx: idx + n_levels[0]]
    idx += n_levels[0]
    f2 = main_effects[idx: idx + n_levels[1]]
    idx += n_levels[1]
    f3 = main_effects[idx: idx + n_levels[2]]

    names = f1 + f2 + f3

    # 2-way interactions
    for a in f1:
        for b in f2 + f3:
            names.append(f"{a}:{b}")
    for a in f2:
        for b in f3:
            names.append(f"{a}:{b}")

    # 3-way interactions
    for a in f1:
        for b in f2:
            for c in f3:
                names.append(f"{a}:{b}:{c}")
    return names


def generate_names_4factors(main_effects, n_levels):
    idx = 0
    f1 = main_effects[idx: idx + n_levels[0]]
    idx += n_levels[0]
    f2 = main_effects[idx: idx + n_levels[1]]
    idx += n_levels[1]
    f3 = main_effects[idx: idx + n_levels[2]]
    idx += n_levels[2]
    f4 = main_effects[idx: idx + n_levels[3]]

    names = f1 + f2 + f3 + f4

    # 2-way interactions
    for a in f1:
        for b in f2 + f3 + f4:
            names.append(f"{a}:{b}")
    for a in f2:
        for b in f3 + f4:
            names.append(f"{a}:{b}")
    for a in f3:
        for b in f4:
            names.append(f"{a}:{b}")

    # 3-way interactions
    for a in f1:
        for b in f2:
            for c in f3 + f4:
                names.append(f"{a}:{b}:{c}")
        for b in f3:
            for c in f4:
                names.append(f"{a}:{b}:{c}")
    for a in f2:
        for b in f3:
            for c in f4:
                names.append(f"{a}:{b}:{c}")

    # 4-way interactions (explicit C-order)
    for a in f1:
        for b in f2:
            for c in f3:
                for d in f4:
                    names.append(f"{a}:{b}:{c}:{d}")

    return names


def generate_names_5factors(main_effects, n_levels):
    idx = 0
    f1 = main_effects[idx: idx + n_levels[0]]
    idx += n_levels[0]
    f2 = main_effects[idx: idx + n_levels[1]]
    idx += n_levels[1]
    f3 = main_effects[idx: idx + n_levels[2]]
    idx += n_levels[2]
    f4 = main_effects[idx: idx + n_levels[3]]
    idx += n_levels[3]
    f5 = main_effects[idx: idx + n_levels[4]]

    names = f1 + f2 + f3 + f4 + f5

    # 2-way interactions
    for a in f1:
        for b in f2 + f3 + f4 + f5:
            names.append(f"{a}:{b}")
    for a in f2:
        for b in f3 + f4 + f5:
            names.append(f"{a}:{b}")
    for a in f3:
        for b in f4 + f5:
            names.append(f"{a}:{b}")
    for a in f4:
        for b in f5:
            names.append(f"{a}:{b}")

    # 3-way interactions
    for a in f1:
        for b in f2:
            for c in f3 + f4 + f5:
                names.append(f"{a}:{b}:{c}")
        for b in f3:
            for c in f4 + f5:
                names.append(f"{a}:{b}:{c}")
        for b in f4:
            for c in f5:
                names.append(f"{a}:{b}:{c}")
    for a in f2:
        for b in f3:
            for c in f4 + f5:
                names.append(f"{a}:{b}:{c}")
        for b in f4:
            for c in f5:
                names.append(f"{a}:{b}:{c}")
    for a in f3:
        for b in f4:
            for c in f5:
                names.append(f"{a}:{b}:{c}")

    # 4-way interactions (explicit C-order)
    for a in f1:
        for b in f2:
            for c in f3:
                for d in f4 + f5:
                    names.append(f"{a}:{b}:{c}:{d}")
            for c in f4:
                for d in f5:
                    names.append(f"{a}:{b}:{c}:{d}")
        for b in f3:
            for c in f4:
                for d in f5:
                    names.append(f"{a}:{b}:{c}:{d}")
    for a in f2:
        for b in f3:
            for c in f4:
                for d in f5:
                    names.append(f"{a}:{b}:{c}:{d}")

    # 5-way interactions  - FOR the moment are not used
    #for a in f1:
        #for b in f2:
            #for c in f3:
                #for d in f4:
                    #for e in f5:
                        #names.append(f"{a}:{b}:{c}:{d}:{e}")

    return names






############### Get identif from X_main ######


# ----------------------
# 2-factor model
# ----------------------
def get_X_2way2(X, l1, l2):
    """
    2-way interactions for 2 factors.
    """
    r1 = list(range(l1))
    r2 = list(range(l1, l1+l2))

    n = X.shape[0]
    n_cols = l1*l2
    X_2way = np.zeros((n, n_cols))
    counter = 0

    for i in r1:
        for j in r2:
            X_2way[:, counter] = X[:, i] * X[:, j]
            counter += 1

    assert counter == n_cols
    return X_2way

# ----------------------
# 3-factor model
# ----------------------
def get_X_2way3(X, l1, l2, l3):
    """
    2-way interactions for 3 factors.
    """
    r1 = list(range(l1))
    r2 = list(range(l1, l1+l2))
    r3 = list(range(l1+l2, l1+l2+l3))

    n = X.shape[0]
    n_cols = l1*(l2+l3) + l2*l3
    X_2way = np.zeros((n, n_cols))
    counter = 0

    for i in r1:
        for j in r2+r3:
            X_2way[:, counter] = X[:, i]*X[:, j]
            counter += 1
    for i in r2:
        for j in r3:
            X_2way[:, counter] = X[:, i]*X[:, j]
            counter += 1

    assert counter == n_cols
    return X_2way


def get_X_3way3(X, l1, l2, l3):
    """
    3-way interactions for 3 factors.
    """
    r1 = list(range(l1))
    r2 = list(range(l1, l1+l2))
    r3 = list(range(l1+l2, l1+l2+l3))

    n = X.shape[0]
    n_cols = l1*l2*l3
    X_3way = np.zeros((n, n_cols))
    counter = 0

    for i in r1:
        for j in r2:
            for k in r3:
                X_3way[:, counter] = X[:, i]*X[:, j]*X[:, k]
                counter += 1

    assert counter == n_cols
    return X_3way

# ----------------------
# 4-factor model
# ----------------------
def get_X_2way4(X, l1, l2, l3, l4):
    """
    2-way interactions for 4 factors.
    """
    r1 = list(range(l1))
    r2 = list(range(l1, l1+l2))
    r3 = list(range(l1+l2, l1+l2+l3))
    r4 = list(range(l1+l2+l3, l1+l2+l3+l4))

    n = X.shape[0]
    n_cols = l1*(l2+l3+l4) + l2*(l3+l4) + l3*l4
    X_2way = np.zeros((n, n_cols))
    counter = 0

    for i in r1:
        for j in r2+r3+r4:
            X_2way[:, counter] = X[:, i]*X[:, j]
            counter += 1
    for i in r2:
        for j in r3+r4:
            X_2way[:, counter] = X[:, i]*X[:, j]
            counter += 1
    for i in r3:
        for j in r4:
            X_2way[:, counter] = X[:, i]*X[:, j]
            counter += 1

    assert counter == n_cols
    return X_2way


def get_X_3way4(X, l1, l2, l3, l4):
    """
    3-way interactions for 4 factors.
    """
    r1 = list(range(l1))
    r2 = list(range(l1, l1+l2))
    r3 = list(range(l1+l2, l1+l2+l3))
    r4 = list(range(l1+l2+l3, l1+l2+l3+l4))

    n = X.shape[0]
    n_cols = l1*l2*(l3+l4) + l1*l3*l4 + l2*l3*l4
    X_3way = np.zeros((n, n_cols))
    counter = 0

    # i in r1, j in r2, k in r3+r4
    for i in r1:
        for j in r2:
            for k in r3+r4:
                X_3way[:, counter] = X[:, i]*X[:, j]*X[:, k]
                counter += 1
    # i in r1, j in r3, k in r4
    for i in r1:
        for j in r3:
            for k in r4:
                X_3way[:, counter] = X[:, i]*X[:, j]*X[:, k]
                counter += 1
    # i in r2, j in r3, k in r4
    for i in r2:
        for j in r3:
            for k in r4:
                X_3way[:, counter] = X[:, i]*X[:, j]*X[:, k]
                counter += 1

    assert counter == n_cols
    return X_3way


def get_X_4way4(X, l1, l2, l3, l4):
    """
    4-way interactions for 4 factors.
    """
    r1 = list(range(l1))
    r2 = list(range(l1, l1+l2))
    r3 = list(range(l1+l2, l1+l2+l3))
    r4 = list(range(l1+l2+l3, l1+l2+l3+l4))

    n = X.shape[0]
    n_cols = l1*l2*l3*l4
    X_4way = np.zeros((n, n_cols))
    counter = 0

    for i in r1:
        for j in r2:
            for k in r3:
                for l in r4:
                    X_4way[:, counter] = X[:, i]*X[:, j]*X[:, k]*X[:, l]
                    counter += 1

    assert counter == n_cols
    return X_4way


def get_identifiable_mains(X_main, levels):
    """
    Converts full dummy-coded factors into effect coding , 
    where each factor has l_i dummy columns. For each factor, if the 
    'last level' dummy column equals 1, the other columns of the same 
    factor are set to -1. After all replacements, the last-level columns 
    are dropped.

    Parameters
    ----------
    X_main : pd.DataFrame
        Fully dummy-encoded dataframe.
    levels : list of int
        Number of dummy columns for each factor.
    Returns
    -------
    pd.DataFrame
        DataFrame with effect-coded factors.
    """
    drop_cols = []
    col_start = 0

    for l in levels:
        col_end = col_start + l-1  # one past the end
        
        # Replace -1 where last_col is 1
        rows_to_replace = X_main.iloc[:, col_end] == 1
        X_main.iloc[rows_to_replace, col_start:col_end] = -1
        
        # Mark last column for dropping later
        drop_cols.append(col_end)

        # Move to next factor block
        col_start = col_end+1

    # Now drop all last dummy columns at once (no shifting issues)
    X_main = X_main.drop(X_main.columns[drop_cols], axis=1)

    return X_main



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


def get_X_all_identifiable(X, levels):
    """
    Returns the full design matrix including main effects and all interaction terms
    for a 2–5 factor model.

    Parameters
    ----------
    X : np.ndarray or pd.DataFrame
        Dummy-coded MAIN effects matrix.
    levels : list of int
        Number of levels for each factor (length determines number of factors).

    Returns
    -------
    np.ndarray
        Design matrix with main effects, 2-way, 3-way, ..., up to n-way interactions.
    """
    n_factors = len(levels)
    assert 2 <= n_factors <= 5, f"Only 2–5 factors supported, got {n_factors}."
    assert X.shape[1] == np.sum(np.array(levels))
    # Convert levels to l_i - 1 for identifiable coding
    l = [lvl - 1 for lvl in levels]

    # Step 1: main effects in identifiable coding
    X_main = get_identifiable_mains(X, levels)

    # Step 2: 2-way interactions
    if n_factors == 2:
        X_2way = get_X_2way2(X_main, l[0], l[1])
        return np.hstack([X_main, X_2way])

    elif n_factors == 3:
        X_2way = get_X_2way3(X_main, l[0], l[1], l[2])
        X_3way = get_X_3way3(X_main, l[0], l[1], l[2])
        return np.hstack([X_main, X_2way, X_3way])

    elif n_factors == 4:
        X_2way = get_X_2way4(X_main, l[0], l[1], l[2], l[3])
        X_3way = get_X_3way4(X_main, l[0], l[1], l[2], l[3])
        X_4way = get_X_4way4(X_main, l[0], l[1], l[2], l[3])
        return np.hstack([X_main, X_2way, X_3way, X_4way])

    elif n_factors == 5:
        X_2way = get_X_2way5(X_main, l[0], l[1], l[2], l[3], l[4])
        X_3way = get_X_3way5(X_main, l[0], l[1], l[2], l[3], l[4])
        X_4way = get_X_4way5(X_main, l[0], l[1], l[2], l[3], l[4])
        X_5way = get_X_5way5(X_main, l[0], l[1], l[2], l[3], l[4])
        return np.hstack([X_main, X_2way, X_3way, X_4way, X_5way])


def recover_coefs_with_names(main_names, beta, intercept, levels):
    
    n_levels=len(levels)
    func_name = f"generate_names_{n_levels}factors"
    all_names = ["intercept"]+ globals()[func_name](main_names, [l+1 for l in levels])
    beta_all = recover_coefs(beta_all=beta, levels=levels)["all"]
    coefs_all=np.insert(beta_all,0,intercept)
    coefs_named = pd.Series(coefs_all, index=all_names, name="Values").to_frame()
    return coefs_named
