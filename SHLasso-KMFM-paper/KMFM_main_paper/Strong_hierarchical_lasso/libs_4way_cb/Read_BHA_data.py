import numpy as np
import pandas as pd

def read_load_BHA_data():
    # --- Read data ---

    dx_train = pd.read_csv("../Analysis/Data/XTrainNScaled.txt", sep=r"\s+", header=0)
    dx_test  = pd.read_csv("../Analysis/Data/XTestNScaled.txt", sep=r"\s+", header=0)
    dy_train = pd.read_csv("../Analysis/Data/YTrain.txt", sep=r"\s+", header=None)
    dy_test  = pd.read_csv("../Analysis/Data/YTest.txt", sep=r"\s+", header=None)


    # Response variable
    y_all = pd.concat([dy_train, dy_test], axis=0).to_numpy()
    n_all = y_all.shape[0]
    tt = np.arange(1, n_all+1)/n_all

    # Add three artificial descriptors for additives
    xc = pd.concat([dx_train, dx_test], axis=0).reset_index(drop=True)
    np.random.seed(2134)
    n_c1 = xc.iloc[:,0].nunique()
    n_c2 = xc.iloc[:,2].nunique()
    n_c3 = xc.iloc[:,3].nunique()
    c1 = pd.Categorical(xc.iloc[:,0]).codes
    c2 = pd.Categorical(xc.iloc[:,2]).codes
    c3 = pd.Categorical(xc.iloc[:,3]).codes
    c1 = c1 + np.random.uniform(size=n_c1)[c1]
    c2 = c2 + np.random.normal(size=n_c2)[c2]
    c3 = c3 + np.random.uniform(size=n_c3)[c3]
    xc = pd.concat([pd.DataFrame({'add_new1': c1, 'add_new2': c2, 'add_new3': c3}), xc], axis=1)

    xs = xc.iloc[:, [3, 22, 49, 59]]  # 0-based indexing
    xs.columns = ["additive", "aryl_halide", "base", "ligand"]

    # Number of levels
    a = np.array([xs.iloc[:, i].nunique() for i in range(xs.shape[1])])
    xcf = np.full((xs.shape[0], np.sum(a)), np.nan)
    b = np.cumsum(a)
    
    # Column names
    colnames = []
    for i in range(len(a)):
        colnames.extend([xs.columns[i]] * a[i])
    colnum = list(np.argsort(xs.iloc[:,0].unique()))
    for i in range(1, xs.shape[1]):
        colnum.extend(list(np.argsort(xs.iloc[:,i].unique())))
    colnames = [f"{col}_{num+1}" for col, num in zip(colnames, colnum)]
    
    # One-hot encoding
    for i in range(xs.shape[0]):
        for j in range(len(a)):
            res = np.zeros(a[j])
            where = np.where(xs.iloc[i,j] == xs.iloc[:,j].unique())[0][0]
            res[where] = 1
            start_idx = 0 if j==0 else b[j-1]
            end_idx = b[j]
            xcf[i, start_idx:end_idx] = res
    for i in range(len(b)):
        ind = xcf[:, b[i]-1] == 1
        start_idx = 0 if i==0 else b[i-1]
        xcf[ind, start_idx:b[i]] = -1
    xcf = np.delete(xcf, b-1, axis=1)
    x_all = xcf.copy()

    # Label for yield
    y_labels = []
    for ijkl in range(n_all):
        nz_idx = np.where(x_all[ijkl,:]!=0)[0]
        add_I = sum(["additive" in colnames[idx] for idx in nz_idx])
        add_i = "additive 22" if add_I>1 else colnames[nz_idx[0]]
        ary_J = sum(["aryl_halide" in colnames[idx] for idx in nz_idx])
        ary_j = "aryl_halide 15" if ary_J>1 else colnames[nz_idx[0+add_I-1]]
        bas_K = sum(["base" in colnames[idx] for idx in nz_idx])
        bas_k = "base 1" if bas_K>1 else colnames[nz_idx[0+add_I+ary_J-1]]
        lig_L = sum(["ligand" in colnames[idx] for idx in nz_idx])
        lig_l = "ligand 3" if lig_L>1 else colnames[nz_idx[0+add_I+ary_J+bas_K-1]]
        y_labels.append(f"{add_i}:{ary_j}:{bas_k}:{lig_l}")
    y_labels = np.array(y_labels)

    # --- Mixed terms 2-way ---
    xx = np.ones((xcf.shape[0],1))
    bb = np.cumsum(a-1)
    for j in range(3):
        start_j = 0 if j==0 else bb[j-1]
        end_j = bb[j]
        for i in range(start_j, end_j):
            xxp = xcf[:,i].reshape(-1,1) * np.delete(xcf, np.s_[:end_j], axis=1)
            xx = np.hstack((xx, xxp))
    xx = np.hstack((xcf, xx[:,1:]))
    xx_all = xx.copy()

    # --- Mixed terms 3-way ---
    xcf1 = x_all.copy()
    # Dynamic assignment for column names
    n_levels = [xs.iloc[:, i].nunique() for i in range(xs.shape[1])]
    xcf1_cols = ["additive"]*n_levels[0] + ["aryl_halide"]*n_levels[1] + ["base"]*n_levels[2] + ["ligand"]*n_levels[3]
    xcf1 = pd.DataFrame(xcf1, columns=xcf1_cols)
    xx1 = xx[:, xcf1.shape[1]:]
    xxx = np.ones((xcf.shape[0],1))
    ind = np.ones(xx1.shape[1], dtype=bool)
    for j in range(2):
        start_j = 0 if j==0 else bb[j-1]
        end_j = bb[j]
        cols_to_check = xcf1.columns[start_j:end_j]
        ind = np.array([not any(c in cols_to_check for c in xcf1.columns)] * len(ind))
        for i in range(start_j, end_j):
            xxxp = xcf[:,i].reshape(-1,1) * xx1[:, ind]
            xxx = np.hstack((xxx, xxxp))
    xxx = np.hstack((xx, xxx[:,1:]))
    xxx_all = xxx.copy()

    # --- Mixed terms 4-way ---
    xxx1 = xxx[:, 40:]  # adjust based on actual 3-way columns
    xxxx = np.ones((xcf.shape[0],1))
    for i in range(21):  # additive levels
        xxxxp = xcf[:,i].reshape(-1,1) * xxx1[:,1597:1680]  # same index as R
        xxxx = np.hstack((xxxx, xxxxp))
    xxxx = np.hstack((xxx, xxxx[:,1:]))
    xxxx_all = xxxx.copy()

    # --- Train/Test split ---
    fac_1234 = y_labels
    fac_1234_train = fac_1234[:dy_train.shape[0]]
    fac_1234_test = fac_1234[dy_train.shape[0]:]

    # Train data
    x_train = xxxx_all[np.isin(fac_1234, fac_1234_train), :]
    y_train = y_all[np.isin(fac_1234, fac_1234_train), :]
    x_test = xxxx_all[np.isin(fac_1234, fac_1234_test), :]
    y_test = y_all[np.isin(fac_1234, fac_1234_test), :]

    print("BHA data is loaded.")
    return {
        "xxxx_all": xxxx_all,
        "y_all": y_all,
        "x_train": x_train,
        "x_test": x_test,
        "y_train": y_train,
        "y_test": y_test
    }
