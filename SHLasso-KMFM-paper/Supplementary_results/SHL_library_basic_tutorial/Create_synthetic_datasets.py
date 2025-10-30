import numpy as np
import pandas as pd
from itertools import product


# Kappa functions
def kappa0(x):
    x = np.asarray(x, dtype=float)
    k0 = np.zeros_like(x)
    
    tt = (x <= 700) & (x != 0)
    k0[tt] = np.log((np.exp(x[tt]) - 1) / x[tt])
    
    tt = x > 700
    k0[tt] = x[tt] - np.log(x[tt])
    
    return k0


def kappa1(x):
    x = np.asarray(x, dtype=float)
    k1 = np.full_like(x, 0.5)
    
    tt = np.abs(x) <= 0.0001
    k1[tt] = 0.5 + x[tt]/12 - x[tt]**3/720 + x[tt]**5/30240
    
    tt = np.abs(x) > 0.0001
    k1[tt] = 1 - 1/x[tt] - 1/(1 - np.exp(x[tt]))
    
    return k1


#def kappa2(x):
#    x = np.asarray(x, dtype=float)
#    k2 = np.full_like(x, 1/12)

#    tt = np.abs(x) <= 0.015
#    k2[tt] = 1/12 - x[tt]**2/240 + x[tt]**4/6048
    
#    tt = np.abs(x) > 0.015
#    k2[tt] = 1/x[tt]**2 + 1/(2 - 2 * np.cosh(x[tt]))
    
#    return k2

def kappa2(x):
    x = np.asarray(x, dtype=float)
    k2 = np.full_like(x, 1/12)
    
    # small x: Taylor expansion
    tt = np.abs(x) <= 0.015
    k2[tt] = 1/12 - x[tt]**2/240 + x[tt]**4/6048
    
    # clip large x to a safe range to avoid overflow
    tt = np.abs(x) > 0.015
    x_tt = np.clip(x[tt], -100, 100)  # prevents cosh overflow
    k2[tt] = 1/x_tt**2 + 1/(2 - 2 * np.cosh(x_tt))
    
    return k2


def _outer_paste(xs, ys, sep=":"):
    """Replicates sort(as.vector(outer(xs, ys, paste, sep=':'))) when followed by sort."""
    return [f"{x}{sep}{y}" for x in xs for y in ys]

def _sorted_outer_paste(xs, ys, sep=":"):
    return sorted(_outer_paste(xs, ys, sep=sep))

def _factor_labels(tag, L):
    """Return ['A.1', ..., 'A.L'] for tag='A' and L levels."""
    return [f"{tag}.{i}" for i in range(1, L + 1)]

def _split_row_components(row_labels, NF):
    """
    Given row labels like 'A.1:B.2' (NF=2) or 'A.1:B.2:C.3' (NF=3), etc.,
    return a list of lists, one per position.
    E.g., for NF=3 -> [first_components, second_components, third_components]
    """
    comps = [r.split(":") for r in row_labels]
    pos_lists = [[] for _ in range(NF)]
    for c in comps:
        for k in range(NF):
            pos_lists[k].append(c[k])
    return pos_lists  # list length = NF, each length = len(row_labels)

# -------------------------------
# Main function: dummy_matrix
# -------------------------------

def dummy_matrix(NF=2, NL=None):
    """
    Python translation of the R function `dummy.matrix`.
    - NF: integer in {2, 3, 4}
    - NL: list/array of length NF with integers in [2, 100]
    Returns a pandas.DataFrame with the same row/column labels and coding as the R version.
    """
    if NL is None:
        NL = [2] * NF
    if NF not in (2, 3, 4):
        raise ValueError("NF must be 2, 3, or 4.")
    if len(NL) != NF:
        raise ValueError("NL must have length NF.")
    if any((n < 2 or n > 100) for n in NL):
        raise ValueError("Each entry of NL must be between 2 and 100.")

    # Factors and Levels (mirror R tags)
    fac_tags = ["A", "B", "C", "D"]
    # fac_levs = as.character(1:100) in R; we generate as needed via _factor_labels

    if NF == 2:
        # One-way
        L1, L2 = NL
        fac1 = _factor_labels(fac_tags[0], L1)  # A.1 ... A.L1
        fac2 = _factor_labels(fac_tags[1], L2)  # B.1 ... B.L2

        # Two-ways (full grid)
        fac12 = _sorted_outer_paste(fac1, fac2, sep=":")

        # Dummy design matrix 2-way
        n_2w = L1 * L2
        p_2w = n_2w - 1

        # Column names: A (excluding last), B (excluding last), then interactions among those (sorted)
        main_A = fac1[:-1]
        main_B = fac2[:-1]
        inter_AB = _sorted_outer_paste(main_A, main_B, sep=":")

        colnames = list(main_A) + list(main_B) + list(inter_AB)

        # Initialize DataFrame
        X = pd.DataFrame(0, index=fac12, columns=colnames, dtype=int)

        # Precompute row components
        first_comp, second_comp = _split_row_components(X.index.tolist(), 2)

        # Fill columns
        for col in X.columns:
            col_tags = col.split(":")
            if len(col_tags) == 1:
                # main effect
                tag = col_tags[0]                # e.g., 'A.1' or 'B.2'
                fac_letter = tag.split(".")[0]   # 'A' or 'B'

                # Set +1 wherever the tag appears in either position (grepl on first and second)
                mask_first = [c == tag for c in first_comp]
                mask_second = [c == tag for c in second_comp]
                X.loc[np.array(mask_first) | np.array(mask_second), col] = 1

                # Set -1 on the last level of the *same factor letter* in each position
                # This reproduces the exact R logic (it writes to both positions; one will match)
                last_tag_first = f"{fac_letter}.{L1}"
                last_tag_second = f"{fac_letter}.{L2}"

                X.loc[np.array([c == last_tag_first for c in first_comp]), col] = -1
                X.loc[np.array([c == last_tag_second for c in second_comp]), col] = -1

            elif len(col_tags) == 2:
                # interaction: product of the two corresponding main-effect columns
                c1, c2 = col_tags
                X[col] = X[c1] * X[c2]

        return X

    if NF == 3:
        # One-way
        L1, L2, L3 = NL
        fac1 = _factor_labels(fac_tags[0], L1)  # A
        fac2 = _factor_labels(fac_tags[1], L2)  # B
        fac3 = _factor_labels(fac_tags[2], L3)  # C

        # Two-ways (full)
        fac12 = _sorted_outer_paste(fac1, fac2, sep=":")
        fac13 = _sorted_outer_paste(fac1, fac3, sep=":")
        fac23 = _sorted_outer_paste(fac2, fac3, sep=":")

        # Three-ways: outer(A, fac23)
        fac123 = _sorted_outer_paste(fac1, fac23, sep=":")

        # Dummy design matrix 3-way
        n_3w = L1 * L2 * L3
        p_3w = n_3w - 1

        main_A = fac1[:-1]
        main_B = fac2[:-1]
        main_C = fac3[:-1]

        inter_AB = _sorted_outer_paste(main_A, main_B, sep=":")
        inter_AC = _sorted_outer_paste(main_A, main_C, sep=":")
        inter_BC = _sorted_outer_paste(main_B, main_C, sep=":")

        # Three-way interactions with A first, then all BC pairs
        inter_ABC = _sorted_outer_paste(main_A, _sorted_outer_paste(main_B, main_C, sep=":"), sep=":")

        colnames = list(main_A) + list(main_B) + list(main_C) \
                 + list(inter_AB) + list(inter_AC) + list(inter_BC) \
                 + list(inter_ABC)

        X = pd.DataFrame(0, index=fac123, columns=colnames, dtype=int)

        # Precompute row components (A, B, C by position)
        compA, compB, compC = _split_row_components(X.index.tolist(), 3)

        for col in X.columns:
            col_tags = col.split(":")
            if len(col_tags) == 1:
                tag = col_tags[0]
                fac_letter = tag.split(".")[0]

                # +1 wherever tag appears in any position
                maskA = np.array([c == tag for c in compA])
                maskB = np.array([c == tag for c in compB])
                maskC = np.array([c == tag for c in compC])
                X.loc[maskA | maskB | maskC, col] = 1

                # -1 on last level of same factor letter at each position
                X.loc[np.array([c == f"{fac_letter}.{L1}" for c in compA]), col] = -1
                X.loc[np.array([c == f"{fac_letter}.{L2}" for c in compB]), col] = -1
                X.loc[np.array([c == f"{fac_letter}.{L3}" for c in compC]), col] = -1

            elif len(col_tags) == 2:
                c1, c2 = col_tags
                X[col] = X[c1] * X[c2]

            elif len(col_tags) == 3:
                c1, c2, c3 = col_tags
                X[col] = X[c1] * X[c2] * X[c3]

        return X

    if NF == 4:
        # One-way
        L1, L2, L3, L4 = NL
        fac1 = _factor_labels(fac_tags[0], L1)  # A
        fac2 = _factor_labels(fac_tags[1], L2)  # B
        fac3 = _factor_labels(fac_tags[2], L3)  # C
        fac4 = _factor_labels(fac_tags[3], L4)  # D

        # Two-ways (full)
        fac12 = _sorted_outer_paste(fac1, fac2, sep=":")
        fac13 = _sorted_outer_paste(fac1, fac3, sep=":")
        fac14 = _sorted_outer_paste(fac1, fac4, sep=":")
        fac23 = _sorted_outer_paste(fac2, fac3, sep=":")
        fac24 = _sorted_outer_paste(fac2, fac4, sep=":")
        fac34 = _sorted_outer_paste(fac3, fac4, sep=":")

        # Three-ways (full)
        fac123 = _sorted_outer_paste(fac1, fac23, sep=":")
        fac124 = _sorted_outer_paste(fac1, fac24, sep=":")
        fac134 = _sorted_outer_paste(fac1, fac34, sep=":")
        fac234 = _sorted_outer_paste(fac2, fac34, sep=":")

        # Four-ways
        fac1234 = _sorted_outer_paste(fac1, fac234, sep=":")

        # Dummy design matrix 4-way
        n_4w = L1 * L2 * L3 * L4
        p_4w = n_4w - 1

        main_A = fac1[:-1]
        main_B = fac2[:-1]
        main_C = fac3[:-1]
        main_D = fac4[:-1]

        inter_AB = _sorted_outer_paste(main_A, main_B, sep=":")
        inter_AC = _sorted_outer_paste(main_A, main_C, sep=":")
        inter_AD = _sorted_outer_paste(main_A, main_D, sep=":")
        inter_BC = _sorted_outer_paste(main_B, main_C, sep=":")
        inter_BD = _sorted_outer_paste(main_B, main_D, sep=":")
        inter_CD = _sorted_outer_paste(main_C, main_D, sep=":")

        inter_ABC = _sorted_outer_paste(main_A, _sorted_outer_paste(main_B, main_C, sep=":"), sep=":")
        inter_ABD = _sorted_outer_paste(main_A, _sorted_outer_paste(main_B, main_D, sep=":"), sep=":")
        inter_ACD = _sorted_outer_paste(main_A, _sorted_outer_paste(main_C, main_D, sep=":"), sep=":")
        inter_BCD = _sorted_outer_paste(main_B, _sorted_outer_paste(main_C, main_D, sep=":"), sep=":")

        inter_ABCD = _sorted_outer_paste(
            main_A,
            _sorted_outer_paste(
                main_B,
                _sorted_outer_paste(main_C, main_D, sep=":"),
                sep=":"
            ),
            sep=":"
        )

        colnames = (
            list(main_A) + list(main_B) + list(main_C) + list(main_D)
            + list(inter_AB) + list(inter_AC) + list(inter_AD)
            + list(inter_BC) + list(inter_BD) + list(inter_CD)
            + list(inter_ABC) + list(inter_ABD) + list(inter_ACD) + list(inter_BCD)
            + list(inter_ABCD)
        )

        X = pd.DataFrame(0, index=fac1234, columns=colnames, dtype=int)

        # Precompute row components A, B, C, D
        compA, compB, compC, compD = _split_row_components(X.index.tolist(), 4)

        for col in X.columns:
            col_tags = col.split(":")
            if len(col_tags) == 1:
                tag = col_tags[0]
                fac_letter = tag.split(".")[0]

                maskA = np.array([c == tag for c in compA])
                maskB = np.array([c == tag for c in compB])
                maskC = np.array([c == tag for c in compC])
                maskD = np.array([c == tag for c in compD])
                X.loc[maskA | maskB | maskC | maskD, col] = 1

                X.loc[np.array([c == f"{fac_letter}.{L1}" for c in compA]), col] = -1
                X.loc[np.array([c == f"{fac_letter}.{L2}" for c in compB]), col] = -1
                X.loc[np.array([c == f"{fac_letter}.{L3}" for c in compC]), col] = -1
                X.loc[np.array([c == f"{fac_letter}.{L4}" for c in compD]), col] = -1

            elif len(col_tags) == 2:
                c1, c2 = col_tags
                X[col] = X[c1] * X[c2]

            elif len(col_tags) == 3:
                c1, c2, c3 = col_tags
                X[col] = X[c1] * X[c2] * X[c3]

            elif len(col_tags) == 4:
                c1, c2, c3, c4 = col_tags
                X[col] = X[c1] * X[c2] * X[c3] * X[c4]

        return X


# Continuous Bernoulli CDF
def cdf_contbern(x, p):
    if p == 0.5:
        return x
    Bp = 2 * p - 1
    return (p**x * (1 - p) ** (1 - x) + p - 1) / Bp


# Inverse CDF
def inverse_cdf_contbern(x, p):
    if p == 0.5:
        return x
    return (np.log(x * (2 * p - 1) + 1 - p) - np.log(1 - p)) / (np.log(p) - np.log(1 - p))


# Sampling function
def sample_contbern(p):
    p = np.array(p, ndmin=1)
    n = len(p)
    u = np.random.rand(n)
    samples = [inverse_cdf_contbern(ui, pi) for ui, pi in zip(u, p)]
    return np.array(samples)




def create_basic_GLM_4way():
    np.random.seed(124)
    x4w = dummy_matrix(NF=3, NL=[8,7,9])

    l1, l2, l3= 7,6,8

    p4w = x4w.shape[1]
    n4w = p4w + 1
    beta_min, beta_max = 2,5

    beta_true = pd.DataFrame(np.zeros((n4w, 1)), columns=["coeffs"])
    beta_true.index = ["interc"] + list(x4w.columns)

    beta_true["coeffs"] = np.random.uniform(beta_min, beta_max, n4w) * np.random.choice(
        [1, -1], size=n4w, replace=True
    )
    
    beta_true["coeffs"][1:(l1+l2+l3+1)] = beta_true["coeffs"][1:(l1+l2+l3+1)] *2# bigger scale for main

    # Main effects
    col_mains = list(x4w.columns[: l1 + l2 + l3 ])

    # Two-way interactions
    col_theta_good = []
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            col_theta_good.append(f"A.{i}:B.{j}")
        for k in range(1, l3 + 1):
            col_theta_good.append(f"A.{i}:C.{k}")
    for j in range(1, l2 + 1):
        for k in range(1, l3 + 1):
            col_theta_good.append(f"B.{j}:C.{k}")


    # Three-way interactions
    col_psi_good = []
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            for k in range(1, l3 + 1):
                col_psi_good.append(f"A.{i}:B.{j}:C.{k}")


    col_all_good = col_mains + col_theta_good + col_psi_good 
    x4w = x4w[col_all_good]

    beta_true.index = ["interc"] + list(x4w.columns)

    levs_true = [
        "interc",
        "A.1",
        "A.2",
        "B.1",
        "B.2",
        "B.3",
        "C.1",
        "C.2",
        "A.1:B.1",
        "A.1:B.2",
        "A.1:B.3",
        "A.2:B.2",
        "A.1:C.1",
        "A.1:C.2",
        "A.2:C.1"
    ]

    beta_true.loc[~beta_true.index.isin(levs_true), "coeffs"] = 0

    eta = beta_true["coeffs"].iloc[0] + np.dot(x4w.values, beta_true["coeffs"].iloc[1:])
    prob = np.exp(eta) / (np.exp(eta) + 1)
    prob[prob>0.9999]=0.9999
    
    y4w = pd.DataFrame(index=x4w.index)
    y4w["obs"] = sample_contbern(prob)
    y4w["true"] = kappa1(beta_true["coeffs"].iloc[0] + np.dot(x4w.values, beta_true["coeffs"].iloc[1:]))

    if all(beta_true.index[1:] == x4w.columns):
        print("Data loaded properly")

    return {
        "X": x4w.values,
        "y": y4w,
        "beta": beta_true
    }





def create_basic_normal_4way():
    np.random.seed(124)
    x4w = dummy_matrix(NF=4, NL=[7, 6, 4, 7])

    l1, l2, l3, l4 = 6,5,3,6

    p4w = x4w.shape[1]
    n4w = p4w + 1
    beta_min, beta_max = 2, 5

    beta_true = pd.DataFrame(np.zeros((n4w, 1)), columns=["coeffs"])
    beta_true.index = ["interc"] + list(x4w.columns)

    beta_true["coeffs"] = np.random.uniform(beta_min, beta_max, n4w) * np.random.choice(
        [1, -1], size=n4w, replace=True
    )
    beta_true["coeffs"][1:(l1+l2+l3+l4+1)] = beta_true["coeffs"][1:(l1+l2+l3+l4+1)]*2# bigger scale for main

    # Main effects
    col_mains = list(x4w.columns[: l1 + l2 + l3 + l4])

    # Two-way interactions
    col_theta_good = []
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            col_theta_good.append(f"A.{i}:B.{j}")
        for k in range(1, l3 + 1):
            col_theta_good.append(f"A.{i}:C.{k}")
        for l in range(1, l4 + 1):
            col_theta_good.append(f"A.{i}:D.{l}")
    for j in range(1, l2 + 1):
        for k in range(1, l3 + 1):
            col_theta_good.append(f"B.{j}:C.{k}")
        for l in range(1, l4 + 1):
            col_theta_good.append(f"B.{j}:D.{l}")
    for k in range(1, l3 + 1):
        for l in range(1, l4 + 1):
            col_theta_good.append(f"C.{k}:D.{l}")

    # Three-way interactions
    col_psi_good = []
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            for k in range(1, l3 + 1):
                col_psi_good.append(f"A.{i}:B.{j}:C.{k}")
            for l in range(1, l4 + 1):
                col_psi_good.append(f"A.{i}:B.{j}:D.{l}")
        for k in range(1, l3 + 1):
            for l in range(1, l4 + 1):
                col_psi_good.append(f"A.{i}:C.{k}:D.{l}")
    for j in range(1, l2 + 1):
        for k in range(1, l3 + 1):
            for l in range(1, l4 + 1):
                col_psi_good.append(f"B.{j}:C.{k}:D.{l}")

    # Four-way interactions
    col_phi_good = []
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            for k in range(1, l3 + 1):
                for l in range(1, l4 + 1):
                    col_phi_good.append(f"A.{i}:B.{j}:C.{k}:D.{l}")

    col_all_good = col_mains + col_theta_good + col_psi_good + col_phi_good
    x4w = x4w[col_all_good]

    beta_true.index = ["interc"] + list(x4w.columns)

    levs_true = [
        "interc",
        "A.1",
        "A.2",
        "A.5",
        "B.1",
        "B.2",
        "C.1",
        "C.2",
        "D.1",
        "A.1:B.1",
        "A.1:B.2",
        "A.2:B.1",
        "A.1:C.1",
        "A.1:D.1",
        "B.1:C.1",
        "B.1:D.1",
        "B.2:C.1",
        "C.1:D.1",
        "A.1:B.1:C.1",
        "A.1:B.1:D.1",
        "A.1:C.1:D.1"
    ]

    beta_true.loc[~beta_true.index.isin(levs_true), "coeffs"] = 0

    eta = beta_true["coeffs"].iloc[0] + np.dot(x4w.values, beta_true["coeffs"].iloc[1:])

    y4w = pd.DataFrame(index=x4w.index)
    y4w["obs"] = eta + np.random.normal(loc=0.0, scale=2.0, size=len(eta))
    y4w["true"] = eta

    if all(beta_true.index[1:] == x4w.columns):
        print("Normal data loaded properly")

    return {
        "X": x4w.values,
        "y": y4w,
        "beta": beta_true
    }