import numpy as np
from .Functions_for_updates import *

import numpy as np

class SHLasso_5way_CB_fast:
    def __init__(self, X, y, l1=7, l2=4, l3=12, l4=8, l5=4, scale=False, glm_final=True):
        self.self = {}
        self.self["means_X"] = np.mean(X, axis=0)
        self.self["stds_X"] = np.std(X, axis=0)

        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.l4 = l4
        self.l5 = l5
        self.scale = scale
        self.glm_final=glm_final

    def fit(self, X, y,
            lambda_beta, lambda_gamma, lambda_delta, lambda_tau, 
            w_beta=None, w_gamma=None, w_delta=None, w_tau=None,
            standardize=False, tol=1e-2, max_iter=5):

        # STEP 0: scale if needed
        if self.scale:
            print("was scaled")
            X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)

        # Ranges
        range_main, range_theta, range_psi, range_phi, range_omega = get_ranges5(
            l1=self.l1, l2=self.l2, l3=self.l3, l4=self.l4, l5=self.l5
        )

        # Initialize betas
        beta_main_hat = np.zeros(len(range_main))
        beta_2way_hat = np.zeros(len(range_theta))
        beta_3way_hat = np.zeros(len(range_psi))
        beta_4way_hat = np.zeros(len(range_phi))
        interc_init = 0

        X_main = X[:, range_main]
        X_2way = X[:, range_theta]
        X_3way = X[:, range_psi]
        X_4way = X[:, range_phi]

        beta_all_old = np.concatenate([beta_main_hat, beta_2way_hat, beta_3way_hat, beta_4way_hat])
            

        for it in range(max_iter):
            print(f"\n=== ITERATION {it+1} ===")

            # === UPDATE MAIN ===
            res_lasso = irlasso_cb(
                X=X_main, Y=y,
                lambdas=[lambda_beta], w_lambda=w_beta, beta0 =  None,
                centering=False, scaling=False, intercept=True, maxit=20,
                tol=0.5, sd_tol=1e-6, verbose=False,
                C=(  X_4way @ beta_4way_hat +
                  X_3way @ beta_3way_hat + X_2way @ beta_2way_hat).reshape(-1,1)
            )
            beta_main_hat = res_lasso["beta"][1:,0,0]
            interc_init = res_lasso["beta"][0,0,0]
            print("updated intercept and main")
            # === UPDATE 2-WAY ===
            X_2way_hier = get_X_2way5(X_main*np.sign(np.abs(beta_main_hat)), self.l1,self.l2,self.l3,self.l4,self.l5)
            if np.all(X_2way_hier == 0):
                beta_2way_hat = np.zeros_like(range_theta)
            else:
                C = (X_4way @ beta_4way_hat + X_3way @ beta_3way_hat + X_main @ beta_main_hat + interc_init).reshape(-1,1)
                X_2way_hier_reduced, nonzero_cols = drop_zero_columns(X_2way_hier)
                res_lasso = irlasso_cb(
                    X=X_2way_hier_reduced, Y=y,lambdas=[lambda_gamma],w_lambda=w_gamma,
                    beta0= None,
                    centering=False, scaling=False, intercept=False, maxit=20,
                    tol=0.5, sd_tol=1e-6, verbose= False,
                    C = C
                )
                beta_reduced = res_lasso["beta"][:,0,0] # we used interc separately
                beta_2way_hat = np.zeros(X_2way_hier.shape[1])
                beta_2way_hat[nonzero_cols] = beta_reduced
            print("updated 2 way")
            # === UPDATE 3-WAY ===
            mapping_3way_to_2way = get_mapping_3way_to_2way(self.l1, self.l2, self.l3, self.l4, self.l5)

            # hierarchical mask: 3-way column is allowed only if all its corresponding 2-way betas are non-zero
            mask_hier = np.array([np.all(beta_2way_hat[mapping_3way_to_2way[i, :] == 1] != 0)
                    if np.any(mapping_3way_to_2way[i, :] == 1) else False
                    for i in range(mapping_3way_to_2way.shape[0]) ])



            # mask for non-zero columns in X_3way
            mask_nonzero = np.any(X_3way != 0, axis=0) #just killed by alreday 0

            # combined mask: drop original zero columns and columns blocked by hierarchy
            combined_mask = mask_hier & mask_nonzero
            combined_mask = mask_hier.copy()
            #print("len max hier, X, combined: ", sum(mask_hier), sum(mask_nonzero), sum(combined_mask))

            #reduced X_3way and reduced beta0
            X_3way_hier_reduced = X_3way[:, combined_mask]
            beta0_reduced = beta_3way_hat[combined_mask]
            #print(sum(combined_mask), X_3way.shape)
    
            #print("interc init: ", interc_init)
            #C term including contributions from main, 2-way, 4-way, and intercept
            C = (X_4way @ beta_4way_hat + X_2way @ beta_2way_hat + X_main @ beta_main_hat + interc_init).reshape(-1, 1)

            # Check if anything remains
            if X_3way_hier_reduced.shape[1] == 0:
                beta_3way_hat = np.zeros_like(range_psi)
            else:
                res_lasso = irlasso_cb(
                    X=X_3way_hier_reduced, Y=y,
                    lambdas=[lambda_delta], w_lambda=w_delta, 
                    beta0= None,
                    centering=False, scaling=False, intercept=False, maxit=10,
                    tol=0.5, sd_tol=1e-6, verbose=False,
                    C=C)
                # expand reduced beta back to full size
                beta_3way_hat = np.zeros(X_3way.shape[1])
                beta_3way_hat[combined_mask] = res_lasso["beta"][:, 0, 0]

            print("updated 3-way")

            # === UPDATE 4-WAY ===
            mapping_4way_to_3way = get_mapping_4way_to_3way(self.l1,self.l2,self.l3,self.l4,self.l5)
            mask_hier = np.array([np.all(beta_3way_hat[mapping_4way_to_3way[i, :] == 1] != 0)
                    if np.any(mapping_4way_to_3way[i, :] == 1) else False
                    for i in range(mapping_4way_to_3way.shape[0]) ])
            mask_nonzero = np.any(X_4way != 0, axis=0) #just killed by alreday 0

            # combined mask: drop original zero columns and columns blocked by hierarchy
            combined_mask = mask_hier & mask_nonzero
            combined_mask = mask_hier.copy()
            #reduced X_3way and reduced beta0
            X_4way_hier_reduced = X_4way[:, combined_mask]
            beta0_reduced = beta_4way_hat[combined_mask]

            C = (X_3way @ beta_3way_hat + X_2way @ beta_2way_hat + X_main @ beta_main_hat + interc_init).reshape(-1, 1)

            if X_4way_hier_reduced.shape[1] == 0:
                beta_4way_hat = np.zeros_like(range_phi)
            else:
                res_lasso = irlasso_cb(
                    X=X_4way_hier_reduced, Y=y,
                    lambdas=[lambda_tau], w_lambda=w_tau, beta0=None,
                    centering=False, scaling=False, intercept=False, maxit=20,
                    tol=0.5, sd_tol=1e-6, verbose=False,
                    C= C
                )
                beta_4way_hat = np.zeros(X_4way.shape[1])
                beta_4way_hat[combined_mask] = res_lasso["beta"][:, 0, 0]
            print("updated 4way")
            # === 5-WAY ===Not included, assummed all 0


            # Combine all betas
            beta_all_new = np.concatenate([beta_main_hat, beta_2way_hat,
                                           beta_3way_hat, beta_4way_hat])

            # Check convergence
            diff = np.sum(np.abs(beta_all_old - beta_all_new))
            #print(f"Iteration {it+1}, sum(abs(beta_old - beta_new)) = {diff:.6f}")

            if diff <= tol*np.sum(abs(beta_all_old)):
                print(f"Converged at iteration {it+1}")
                break

            beta_all_old = beta_all_new.copy()
        
        
        ############# USE NO PENALTY ON SELECTED NON0 SET ################
        if self.glm_final == True :
            #print("did glm final")
            X_all_selected=X[:,beta_all_new!=0]
            res_glm =  irlasso_cb(X=X_all_selected, Y=y,
                lambdas=[0], w_lambda=None, beta0 = None,
                centering=False, scaling=False, intercept=True, maxit=10,
                tol=0.5, sd_tol=1e-6, verbose=False,
                C=0)
            beta_all_estimated = res_glm["beta"][1:,0,0]
            interc_init = res_glm["beta"][0,0,0]
            beta_all_new[beta_all_new!=0]=beta_all_estimated
            beta_main_hat=beta_all_new[range_main]
            beta_2way_hat=beta_all_new[range_theta]
            beta_3way_hat=beta_all_new[range_psi]
            beta_4way_hat=beta_all_new[range_phi]

        # Save final betas
        self.beta_all = beta_all_new
        self.intercept = interc_init
        self.beta_main_hat = beta_main_hat
        self.beta_2way_hat = beta_2way_hat
        self.beta_3way_hat = beta_3way_hat
        self.beta_4way_hat = beta_4way_hat
    

        return dict(
            intercept=interc_init,
            beta_main_hat=beta_main_hat,
            beta_2way_hat=beta_2way_hat,
            beta_3way_hat=beta_3way_hat,
            beta_4way_hat=beta_4way_hat,
            beta_all=beta_all_new
        )


    def predict(self, X_new, scale=False):
        if scale:
            X_new = (X_new - np.mean(X_new, axis=0)) / np.std(X_new, axis=0)
            print("Take care! X was scaled.")

        v = X_new @ self.beta_all + self.intercept
        y_pred = kappa1(v)  # same as in R
        return y_pred

    def R2_score(self, X_new, y_true, scale=False, verbose=True):
        y_pred = self.predict(X_new, scale=scale)
        score = r2(y_true, y_pred)
        if verbose:
            print("r2 score is ", score)
        return score
