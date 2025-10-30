from sklearn.metrics import r2_score as r2
import numpy as np

def l1_loss_beta(beta, beta_hat, scale=True):
    beta = np.array(beta)       
    beta_hat = np.array(beta_hat)  
    result = np.sum(np.abs(beta - beta_hat))   
    if scale:                 
        result = result / len(beta)   
    return result               

# Example:
# print(l1_loss_beta([-1,2,3,5,0], [-1,0,0,-1,5]))


def MSE_beta(beta, beta_hat):
    beta = np.array(beta)       # beta<-c(beta)
    beta_hat = np.array(beta_hat)
    return np.linalg.norm(beta - beta_hat, 2)**2 / len(beta)

# Example:
# print(MSE_beta([1,2,3,4], [0,0,0,1]))


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