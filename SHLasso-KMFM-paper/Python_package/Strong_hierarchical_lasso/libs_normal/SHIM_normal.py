import numpy as np
from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score
from .Combinatorics import *

# ------------------- 2-way SHIM -------------------
class SHIM_2way_normal:
    def __init__(self, X, y, l1=21, l2=14, scale=False, use_intercept=True, use_offset=True):
        self.means_X = np.mean(X, axis=0)
        self.stds_X = np.std(X, axis=0, ddof=1)
        self.intercept = np.mean(y)
        self.l1 = l1
        self.l2 = l2
        self.scale = scale
        self.use_intercept = use_intercept
        self.use_offset=use_offset

    def fit(self, X, y, lambda_beta, lambda_gamma):
        if self.scale:
            X = (X - self.means_X) / self.stds_X
            print("X was scaled")

        range_main, range_theta = get_ranges2(l1=self.l1, l2=self.l2)
        X_main = X[:, range_main]
        X_2way = X[:, range_theta]

        y_final = y if not self.use_intercept else y - np.mean(y)
        if self.use_intercept:
            self.intercept = np.mean(y)
            print(f"Intercept used: {self.intercept}")

        # Main effects
        lasso_main = Lasso(alpha=lambda_beta, fit_intercept=False, max_iter=10000)
        lasso_main.fit(X_main, y_final)
        self.beta_main_hat = lasso_main.coef_

        # 2-way
        beta_2way_without_gamma = get_beta_vec_2way2(beta=self.beta_main_hat, l1=self.l1, l2=self.l2, gamma=None, only_beta=True)
        X_2way_hier = X_2way * np.abs(np.sign(beta_2way_without_gamma))
        non_zero_cols = np.any(X_2way_hier != 0, axis=0)
        self.beta_2way_hat = np.zeros(len(range_theta))
        if self.use_offset ==True:  #get out contribution of mains
            y_final=y_final-X_main@self.beta_main_hat
        if np.any(non_zero_cols):
            lasso_2way = Lasso(alpha=lambda_gamma, fit_intercept=False, max_iter=10000)
            lasso_2way.fit(X_2way_hier[:, non_zero_cols], y_final)
            self.beta_2way_hat[non_zero_cols] = lasso_2way.coef_

        self.beta_all = np.concatenate([self.beta_main_hat, self.beta_2way_hat])

    def predict(self, X_new, scale=False):
        if scale:
            X_new = (X_new - self.means_X) / self.stds_X
            print("Take care! X was scaled.")
        return X_new @ self.beta_all + self.intercept

    def R2_score(self, X_new, y_true, scale=False, verbose=True):
        y_pred = self.predict(X_new, scale=scale)
        r2_val = r2_score(y_true, y_pred)
        if verbose:
            print("R2 score is", r2_val)
        return r2_val


# ------------------- 3-way SHIM -------------------
class SHIM_3way_normal:
    def __init__(self, X, y, l1, l2, l3, scale=False, use_intercept=True, use_offset=True):
        self.means_X = np.mean(X, axis=0)
        self.stds_X = np.std(X, axis=0, ddof=1)
        self.intercept = np.mean(y)
        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.scale = scale
        self.use_intercept = use_intercept
        self.use_offset=use_offset

    def fit(self, X, y, lambda_beta, lambda_gamma, lambda_delta):
        if self.scale:
            X = (X - self.means_X) / self.stds_X
            print("X was scaled")

        range_main, range_theta, range_psi= get_ranges3(l1=self.l1, l2=self.l2, l3=self.l3)
        X_main = X[:, range_main]
        X_2way = X[:, range_theta]
        X_3way = X[:, range_psi]

        y_final = y if not self.use_intercept else y - np.mean(y)
        if self.use_intercept:
            self.intercept = np.mean(y)
            print(f"Intercept used: {self.intercept}")

        # Main effects
        lasso_main = Lasso(alpha=lambda_beta, fit_intercept=False, max_iter=10000)
        lasso_main.fit(X_main, y_final)
        self.beta_main_hat = lasso_main.coef_

        # 2-way
        if self.use_offset ==True:  #get out contribution of mains
            y_final=y_final-X_main@self.beta_main_hat
        beta_2way_without_gamma = get_beta_vec_2way3(beta=self.beta_main_hat, l1=self.l1, l2=self.l2, l3=self.l3, gamma=None, only_beta=True)
        X_2way_hier = X_2way * np.abs(np.sign(beta_2way_without_gamma))
        non_zero_cols = np.any(X_2way_hier != 0, axis=0)
        self.beta_2way_hat = np.zeros(len(range_theta))
        if np.any(non_zero_cols):
            lasso_2way = Lasso(alpha=lambda_gamma, fit_intercept=False, max_iter=10000)
            lasso_2way.fit(X_2way_hier[:, non_zero_cols], y_final)
            self.beta_2way_hat[non_zero_cols] = lasso_2way.coef_

        # 3-way
        if self.use_offset ==True:  #get out contribution of 2ways
            y_final=y_final-X_2way@self.beta_2way_hat
        beta_3way_without_delta = get_beta_vec_3way3(beta_2way= self.beta_2way_hat,l1= self.l1, l2= self.l2,l3= self.l3, delta=None, only_beta=True)
        X_3way_hier = X_3way * np.abs(np.sign(beta_3way_without_delta))
        non_zero_cols = np.any(X_3way_hier != 0, axis=0)
        self.beta_3way_hat = np.zeros(len(range_psi))
        if np.any(non_zero_cols):
            lasso_3way = Lasso(alpha=lambda_delta, fit_intercept=False, max_iter=10000)
            lasso_3way.fit(X_3way_hier[:, non_zero_cols], y_final)
            self.beta_3way_hat[non_zero_cols] = lasso_3way.coef_

        self.beta_all = np.concatenate([self.beta_main_hat, self.beta_2way_hat, self.beta_3way_hat])

    def predict(self, X_new, scale=False):
        if scale:
            X_new = (X_new - self.means_X) / self.stds_X
            print("Take care! X was scaled.")
        return X_new @ self.beta_all + self.intercept

    def R2_score(self, X_new, y_true, scale=False, verbose=True):
        y_pred = self.predict(X_new, scale=scale)
        r2_val = r2_score(y_true, y_pred)
        if verbose:
            print("R2 score is", r2_val)
        return r2_val


# ------------------- 4-way SHIM -------------------
class SHIM_4way_normal:
    def __init__(self, X, y, l1=21, l2=14, l3=2, l4=3, scale=False, use_intercept=True, use_offset=True):
        self.means_X = np.mean(X, axis=0)
        self.stds_X = np.std(X, axis=0, ddof=1)
        self.intercept = np.mean(y)
        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.l4 = l4
        self.scale = scale
        self.use_intercept = use_intercept
        self.use_offset=use_offset

    def fit(self, X, y, lambda_beta, lambda_gamma, lambda_delta, lambda_tau):
        if self.scale:
            X = (X - self.means_X) / self.stds_X
            print("X was scaled")

        range_main, range_theta, range_psi, range_phi = get_ranges4(l1=self.l1, l2=self.l2, l3=self.l3, l4=self.l4)
        X_main = X[:, range_main]
        X_2way = X[:, range_theta]
        X_3way = X[:, range_psi]
        X_4way = X[:, range_phi]

        y_final = y if not self.use_intercept else y - np.mean(y)
        if self.use_intercept:
            self.intercept = np.mean(y)
            print(f"Intercept used: {self.intercept}")

        # Main effects
        lasso_main = Lasso(alpha=lambda_beta, fit_intercept=False, max_iter=10000)
        lasso_main.fit(X_main, y_final)
        self.beta_main_hat = lasso_main.coef_

        # 2-way
        if self.use_offset ==True:  #get out contribution of mains
            y_final=y_final-X_main@self.beta_main_hat
        beta_2way_without_gamma = get_beta_vec_2way4(beta= self.beta_main_hat, l1=self.l1, l2=self.l2, l3=self.l3, l4=self.l4, gamma=None, only_beta=True)
        X_2way_hier = X_2way * np.abs(np.sign(beta_2way_without_gamma))
        non_zero_cols = np.any(X_2way_hier != 0, axis=0)
        self.beta_2way_hat = np.zeros(len(range_theta))
        if np.any(non_zero_cols):
            lasso_2way = Lasso(alpha=lambda_gamma, fit_intercept=False, max_iter=10000)
            lasso_2way.fit(X_2way_hier[:, non_zero_cols], y_final)
            self.beta_2way_hat[non_zero_cols] = lasso_2way.coef_

        # 3-way
        if self.use_offset ==True:  #get out contribution of 2ways
            y_final=y_final-X_2way@self.beta_2way_hat
        beta_3way_without_delta = get_beta_vec_3way4(beta_2way=self.beta_2way_hat, l1=self.l1, l2=self.l2, l3=self.l3, l4=self.l4, delta=None, only_beta=True)
        X_3way_hier = X_3way * np.abs(np.sign(beta_3way_without_delta))
        non_zero_cols = np.any(X_3way_hier != 0, axis=0)
        self.beta_3way_hat = np.zeros(len(range_psi))
        if np.any(non_zero_cols):
            lasso_3way = Lasso(alpha=lambda_delta, fit_intercept=False, max_iter=10000)
            lasso_3way.fit(X_3way_hier[:, non_zero_cols], y_final)
            self.beta_3way_hat[non_zero_cols] = lasso_3way.coef_

        # 4-way
        if self.use_offset ==True:  #get out contribution of 3ways
            y_final=y_final-X_3way@self.beta_3way_hat
        beta_4way_without_tau = get_beta_vec_4way4(beta_3way=self.beta_3way_hat, tau=None, l1=self.l1, l2=self.l2, l3=self.l3, l4=self.l4, only_beta=True)
        X_4way_hier = X_4way * np.abs(np.sign(beta_4way_without_tau))
        non_zero_cols = np.any(X_4way_hier != 0, axis=0)
        self.beta_4way_hat = np.zeros(len(range_phi))
        if np.any(non_zero_cols):
            lasso_4way = Lasso(alpha=lambda_tau, fit_intercept=False, max_iter=10000)
            lasso_4way.fit(X_4way_hier[:, non_zero_cols], y_final)
            self.beta_4way_hat[non_zero_cols] = lasso_4way.coef_

        self.beta_all = np.concatenate([self.beta_main_hat, self.beta_2way_hat, self.beta_3way_hat, self.beta_4way_hat])
        
    def predict(self, X_new, scale=False):
        if scale:
            X_new = (X_new - self.means_X) / self.stds_X
            print("Take care! X was scaled.")
        return X_new @ self.beta_all + self.intercept

    def R2_score(self, X_new, y_true, scale=False, verbose=True):
        y_pred = self.predict(X_new, scale=scale)
        r2_val = r2_score(y_true, y_pred)
        if verbose:
            print("R2 score is", r2_val)
        return r2_val





class SHIM_5way_normal:
    def __init__(self, X, y, l1=7, l2=4, l3=12, l4=8, l5=4, scale=False, use_offset=True):
        self.self = {}
        self.self["means_X"] = np.mean(X, axis=0)
        self.self["stds_X"] = np.std(X, axis=0)

        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.l4 = l4
        self.l5 = l5
        self.scale = scale

        # placeholders
        self.beta_all = None
        self.intercept = 0.0
        self.beta_main_hat = None
        self.beta_2way_hat = None
        self.beta_3way_hat = None
        self.beta_4way_hat = None
        
        self.use_offset=use_offset

    def fit(self, X, y, lambda_beta, lambda_gamma, lambda_delta, lambda_tau, standardize=False):
        """One-pass hierarchical Lasso update for up to 4-way interactions (5-way ignored)."""

        # Scale if needed
        if self.scale:
            print("was scaled")
            X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)

        # Ranges for features by interaction order
        range_main, range_theta, range_psi, range_phi, range_omega = get_ranges5(
            l1=self.l1, l2=self.l2, l3=self.l3, l4=self.l4, l5=self.l5
        )

        # Extract blocks
        X_main = X[:, range_main]
        X_2way = X[:, range_theta]
        X_3way = X[:, range_psi]
        X_4way = X[:, range_phi]

        # --- MAIN EFFECTS ---
        lasso_main = Lasso(alpha=lambda_beta, fit_intercept=True, max_iter=10000)
        lasso_main.fit(X_main, y)
        beta_main_hat = lasso_main.coef_.copy()
        intercept = float(lasso_main.intercept_)
        print("updated main")

        # --- 2-WAY EFFECTS ---
        if self.use_offset ==True:  #get out contribution of mains
            y_final=y_final-X_main@self.beta_main_hat
        X_2way_hier = get_X_2way5(X_main * np.sign(np.abs(beta_main_hat)),
                                  self.l1, self.l2, self.l3, self.l4, self.l5)
        if np.all(X_2way_hier == 0):
            beta_2way_hat = np.zeros(X_2way.shape[1])
        else:
            X_2way_hier_reduced, nonzero_cols = drop_zero_columns(X_2way_hier)
            if X_2way_hier_reduced.shape[1] == 0:
                beta_2way_hat = np.zeros(X_2way.shape[1])
            else:
                lasso_2way = Lasso(alpha=lambda_gamma, fit_intercept=False, max_iter=10000)
                lasso_2way.fit(X_2way_hier_reduced, y)
                beta_2way_hat = np.zeros(X_2way.shape[1])
                beta_2way_hat[nonzero_cols] = lasso_2way.coef_
        print("updated 2-way")

        # --- 3-WAY EFFECTS ---
        if self.use_offset ==True:  #get out contribution of 2ways
            y_final=y_final-X_2way@self.beta_2way_hat
        mapping_3way_to_2way = get_mapping_3way_to_2way(self.l1, self.l2, self.l3, self.l4, self.l5)
        mask_hier = np.array([
            np.all(beta_2way_hat[mapping_3way_to_2way[i, :] == 1] != 0)
            if np.any(mapping_3way_to_2way[i, :] == 1) else False
            for i in range(mapping_3way_to_2way.shape[0])
        ])
        mask_nonzero = np.any(X_3way != 0, axis=0)
        combined_mask = mask_hier & mask_nonzero

        X_3way_hier_reduced = X_3way[:, combined_mask]
        if X_3way_hier_reduced.shape[1] == 0:
            beta_3way_hat = np.zeros(X_3way.shape[1])
        else:
            lasso_3way = Lasso(alpha=lambda_delta, fit_intercept=False, max_iter=10000)
            lasso_3way.fit(X_3way_hier_reduced, y)
            beta_3way_hat = np.zeros(X_3way.shape[1])
            beta_3way_hat[combined_mask] = lasso_3way.coef_
        print("updated 3-way")

        # --- 4-WAY EFFECTS ---
        if self.use_offset ==True:  #get out contribution of 3ways
            y_final=y_final-X_3way@self.beta_3way_hat
        mapping_4way_to_3way = get_mapping_4way_to_3way(self.l1, self.l2, self.l3, self.l4, self.l5)
        mask_hier = np.array([
            np.all(beta_3way_hat[mapping_4way_to_3way[i, :] == 1] != 0)
            if np.any(mapping_4way_to_3way[i, :] == 1) else False
            for i in range(mapping_4way_to_3way.shape[0])
        ])
        mask_nonzero = np.any(X_4way != 0, axis=0)
        combined_mask = mask_hier & mask_nonzero

        X_4way_hier_reduced = X_4way[:, combined_mask]
        if X_4way_hier_reduced.shape[1] == 0:
            beta_4way_hat = np.zeros(X_4way.shape[1])
        else:
            lasso_4way = Lasso(alpha=lambda_tau, fit_intercept=False, max_iter=10000)
            lasso_4way.fit(X_4way_hier_reduced, y)
            beta_4way_hat = np.zeros(X_4way.shape[1])
            beta_4way_hat[combined_mask] = lasso_4way.coef_
        print("updated 4-way")

        # --- Combine results ---
        beta_all = np.concatenate([beta_main_hat, beta_2way_hat, beta_3way_hat, beta_4way_hat])

        # Save results
        self.beta_all = beta_all
        self.intercept = intercept
        self.beta_main_hat = beta_main_hat
        self.beta_2way_hat = beta_2way_hat
        self.beta_3way_hat = beta_3way_hat
        self.beta_4way_hat = beta_4way_hat

        return dict(
            intercept=intercept,
            beta_main_hat=beta_main_hat,
            beta_2way_hat=beta_2way_hat,
            beta_3way_hat=beta_3way_hat,
            beta_4way_hat=beta_4way_hat,
            beta_all=beta_all
        )

    def predict(self, X_new, scale=False):
        if scale:
            X_new = (X_new - np.mean(X_new, axis=0)) / np.std(X_new, axis=0)
            print("Take care! X was scaled.")

        y_pred = X_new @ self.beta_all + self.intercept
        return y_pred

    def R2_score(self, X_new, y_true, scale=False, verbose=True):
        y_pred = self.predict(X_new, scale=scale)
        score = r2_score(y_true, y_pred)
        if verbose:
            print("r2 score is ", score)
        return score
