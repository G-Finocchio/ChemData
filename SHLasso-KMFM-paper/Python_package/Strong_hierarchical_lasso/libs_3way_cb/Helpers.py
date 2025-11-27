import pandas as pd
import os 
from sklearn.metrics import r2_score as r2
import numpy as np



def l1_loss_beta(beta, beta_hat, scale=True):
    beta = np.array(beta)       
    beta_hat = np.array(beta_hat)  
    result = np.sum(np.abs(beta - beta_hat))   
    if scale:                 
        result = result / len(beta)   
    return result               



def MSE_beta(beta, beta_hat):
    beta = np.array(beta)       # beta<-c(beta)
    beta_hat = np.array(beta_hat)
    return np.linalg.norm(beta - beta_hat, 2)**2 / len(beta)



def TPR_zeros(beta, beta_hat):
    """TPR = sensitivity: correctly predicted zeros from true zeros"""
    beta = np.array(beta)
    beta_hat = np.array(beta_hat)
    idx_beta = np.where(beta == 0)[0]        
    idx_beta_hat = np.where(beta_hat == 0)[0]  
    
    denom = len(idx_beta)
    if denom == 0:   # no true zeros -> undefined TPR
        return np.nan
    
    return len(np.intersect1d(idx_beta, idx_beta_hat)) / denom * 100


def FPR_zeros(beta, beta_hat):
    """FPR = predicted zero when true value != 0"""
    beta = np.array(beta)
    beta_hat = np.array(beta_hat)
    idx_beta = np.where(beta != 0)[0]         
    idx_beta_hat = np.where(beta_hat == 0)[0] 
    
    denom = len(idx_beta)
    if denom == 0:   # no true non-zeros -> undefined FPR
        return np.nan
    
    return len(np.intersect1d(idx_beta, idx_beta_hat)) / denom * 100

# Example:
# print(FPR_zeros([1,1,0,0,1], [1,0,0,0,0]))


def all_beta_functions(beta, beta_hat, scale=True):
    print("TPR zeros:", TPR_zeros(beta, beta_hat))
    print("FPR zeros:", FPR_zeros(beta, beta_hat))
    print("MSE beta:", MSE_beta(beta, beta_hat))
    print("L1 loss beta:", l1_loss_beta(beta, beta_hat, scale=scale))




def load_deoxi_flourination(data_path=os.path.join("../Analysis", "Data")):
    # --- File paths ---
    train_path = os.path.join(data_path, "fluorination_train.csv")
    test_path  = os.path.join(data_path, "fluorination_test.csv")

    agd_palette = ["#0063C3", "#00BCD0", "#CB7C85", "#E8C9BB",
                   "#252A5A", "#7D2C42", "#A9A9A9"]

    # --- Original data ---
    dat_1 = pd.read_csv(train_path)
    dat_2 = pd.read_csv(test_path)
    dat   = pd.concat([dat_1, dat_2], ignore_index=True)

    # Rename and clean
    dat = dat.rename(columns={"sulfonyl.fluoride": "s", "alcohol": "a", "base": "b"})
    dat["s"] = dat["s"].str.replace(".", "-", regex=False)
    dat["s"] = dat["s"].str.replace("3-", "4-", regex=False)
    dat["s"] = dat["s"].astype("category")

    # --- Reaction data ---
    rxns = dat[["a", "b", "s", "yield"]].copy()

    # --- Pre-process factor-like features ---
    da_orig = dat.loc[:, dat.columns.str.contains("alcohol") | (dat.columns == "a")].drop_duplicates()

    da_fact = da_orig.loc[:, ["a"] + da_orig.select_dtypes(include="integer").columns.tolist()].copy()
    da_fact = da_fact.rename(columns=lambda x: x.replace("alcohol_", "").replace(".", "_"))

    # --- Order: primary / secondary / tertiary ---
    order_df = da_fact[["primary", "secondary", "tertiary"]].astype(str)
    order_df["secondary"] = order_df["secondary"].replace({"1": "2"})
    order_df["tertiary"]  = order_df["tertiary"].replace({"1": "3"})
    order = order_df[["primary", "secondary", "tertiary"]].max(axis=1)
    order = order.replace({"1": "primary", "2": "secondary", "3": "tertiary"})

    # --- Ring size ---
    ring_df = da_fact[["4_membered_ring", "5_membered_ring",
                       "6_membered_ring", "7_membered_ring"]].astype(str)
    ring_df["r4"] = ring_df["4_membered_ring"].replace({"1": "4", "0": np.nan})
    ring_df["r5"] = ring_df["5_membered_ring"].replace({"1": "5", "0": np.nan})
    ring_df["r6"] = ring_df["6_membered_ring"].replace({"1": "6", "0": np.nan})
    ring_df["r7"] = ring_df["7_membered_ring"].replace({"1": "7", "0": np.nan})

    # Convert to numeric (np.nan safe) and take min across rows
    ring_size = ring_df.apply(pd.to_numeric, errors="coerce").min(axis=1)
    ring_size = ring_size.fillna(0).astype(int).astype(str)

    da_fact = da_fact.assign(ring_size=ring_size, order=order)
    da_fact = da_fact.drop(columns=[c for c in da_fact.columns if "membered_ring" in c] + 
                                    ["primary", "secondary", "tertiary", "cyclic"])
    da_fact = da_fact.astype("category")

    # --- Alcohol DFT features ---
    glb_feats = ["number_of_atoms", "molar_mass", "electronegativity",
                 "electronic_spatial_extent", "hardness", "homo_energy",
                 "lumo_energy", "molar_volume", "dipole", 
                 "OC_length", "OC_L", "OC_B1", "OC_B5"]
    c_feats = ["C_APT_charge", "C_Mulliken_charge", "C_NMR_shift",
               "C_NPA_charge", "C_VBur", "C_angle", "C_PVBur"]
    o_feats = ["O_APT_charge", "O_Mulliken_charge", "O_NMR_shift",
               "O_NPA_charge", "O_VBur"]

    da_dft_m062x = pd.read_csv(os.path.join(data_path,"alcohols_Boltzmann_M062X_THF.csv"))
    if "inchi" in da_dft_m062x.columns:
        da_dft_m062x = da_dft_m062x.drop(columns=["inchi"])
    da_dft_m062x = da_dft_m062x.rename(columns={"name": "a"})
    da_dft_m062x = da_dft_m062x[["a"] + glb_feats + c_feats + o_feats]

    da_dft = da_dft_m062x
    da = pd.merge(da_dft, da_fact, on="a")

    # --- Normalize yield ---
    def normalize_yield(y):
        y = np.where(y > 97, 100 - np.random.uniform(0.1, 3, size=len(y)),
            np.where(y < 3, np.random.uniform(0.1, 3, size=len(y)), y))
        return y / 100

    # --- Final dataframe ---
    data = rxns.copy()
    data["prob"] = normalize_yield(data["yield"].values)
    data = data.merge(da, on="a", how="left")
    data["a"] = data["a"].astype("category")

    print("Data loaded successfully")
    return data
