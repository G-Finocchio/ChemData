from sklearn.linear_model import Lasso as Lasso_sk
import numpy as np

# Criterion to stop the model
def compute_rel_dif(Q_old, Q_new):
    rel_dif = np.abs(Q_old - Q_new) / np.abs(Q_old)
    if Q_old<Q_new:
        print("Warning: Q_old<Q_new at this step!")
    return rel_dif


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



# Approximation of link function g
def g_link(x):
    x = np.asarray(x, dtype=float)
    # Clamp values between 0.001 and 0.999
    tt = np.clip(x, 0.001, 0.999)
    g = 3.5 * np.tan(np.pi * (2 * tt - 1) / 2)
    return g








def irlasso_cb(X, Y, lambdas, w_lambda=None, beta0=None, centering=True, scaling=True, intercept=True, maxit=10, tol=0.0545, sd_tol=1e-6, verbose=False, C=0):

    if verbose:
        print("Performing IRLASSO-NEW")

    # Convert to numpy arrays
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).reshape(-1, 1)
    
    # Center variables
    mu_X = np.zeros(X.shape[1])
    mu_Y = np.zeros(Y.shape[1])
    if centering:
        mu_X = X.mean(axis=0)
        mu_Y = Y.mean(axis=0)
        X = X - mu_X
        Y = Y - mu_Y

    # Scale variables
    sd_X = np.ones(X.shape[1])
    sd_Y = np.ones(Y.shape[1])
    if scaling:
        sd_X = np.maximum(X.std(axis=0, ddof=1), sd_tol)
        sd_Y = np.maximum(Y.std(axis=0, ddof=1), sd_tol)
        X = X / sd_X
        Y = Y / sd_Y

    # Add intercept column
    if intercept:
        X = np.hstack([np.ones((X.shape[0], 1)), X])

    # Dimensions
    n, p = X.shape
    
    # Offset handling
    if np.isscalar(C):  
        C_vec = np.zeros(n)  # scalar offset = 0 for all observations
    else:
        C_vec = np.asarray(C, dtype=float).reshape(-1)
    if C_vec.shape[0] != n:
        raise ValueError(f"Offset C must have length {n}, got {C_vec.shape[0]}")


    # Initialization
    if beta0 is None:
        beta0 = np.zeros(p)

    beta_old = np.repeat(beta0[:, None, None], len(lambdas), axis=2)
    Z_old = np.zeros((n, 1, len(lambdas)))
    W_old = np.zeros((n, n, len(lambdas)))

    beta = beta_old.copy()
    Z = Z_old.copy()
    W = W_old.copy()

    R = [None] * len(lambdas)
    nc = list(range(len(lambdas)))
    it_stop = np.zeros(len(lambdas), dtype=int)

    counter = 0
    while True:
        if verbose:
            print(f"Performing IRLASSO iter {counter}")

        beta_old = beta.copy()
        Z_old = Z.copy()
        W_old = W.copy()

        eta = np.zeros((n, 1, len(lambdas)))
        k_p = np.zeros_like(eta)
        k_pp = np.zeros_like(eta)

        if len(nc) == 0:
            if verbose:
                print("No lambda left...")
            break

        if verbose:
            print(f"Lambda left {nc}")

        for m in nc:

            eta[:, :, m] = X @ beta[:, :, m] + C_vec[:, None]

            k_p[:, :, m] = kappa1(eta[:, :, m])
            k_pp[:, :, m] = kappa2(eta[:, :, m])
            k_pp[:, :, m][k_pp[:, :, m] < 0.005] = 0.005

            W[:, :, m] = np.diag(k_pp[:, 0, m])
            Z[:, :, m] = eta[:, :, m] + np.diag(1.0 / k_pp[:, 0, m]) @ (Y - k_p[:, :, m])

            # Weighted X and Z
            sqrtW = np.sqrt(W[:, :, m])
            X_W = sqrtW @ X
            Z_W = sqrtW @ Z[:, :, m]

            # Penalty weights
            if w_lambda is None:
                w_lambda = np.ones(X_W.shape[1])
            if intercept:
                w_lambda[0]=1e-4  # Dont penalize intercept if it exists

            if lambdas[m] == 0:
                # === OLS ===
                # Closed-form solution: β = (XᵀX)⁻¹ Xᵀy
                beta[:, 0, m] = np.linalg.pinv(X_W.T @ X_W) @ (X_W.T @ Z_W.ravel())
            else:
                # === Lasso ===
                #model = Lasso(alpha=lambdas[m], fit_intercept=False, weights= w_lambda)
                #model.fit(X_W, Z_W.ravel())
                #beta[:, 0, m] = model.coef_
                model = Lasso_sk(alpha=lambdas[m], fit_intercept=False, max_iter=50000)
                model.fit(X_W, Z_W.ravel())
                beta[:, 0, m] = model.coef_

            # Model selection matrix
            if intercept:
                s_lasso = np.where(beta[:,0,m][1:] != 0)[0]
                R[m] = np.zeros((p, len(s_lasso)))
                if len(s_lasso) > 0:
                    for idx, s in enumerate(s_lasso):
                        R[m][s+1, idx] = 1
                if R[m].shape[1] > 0:
                    R[m][0, 0] = 1
            else:
                s_lasso = np.where(beta[:,0,m] != 0)[0]
                R[m] = np.zeros((p, len(s_lasso)))
                if len(s_lasso) > 0:
                    for idx, s in enumerate(s_lasso):
                        R[m][s, idx] = 1

        # Convergence checks
        epsilon = np.sqrt(((beta - beta_old) ** 2).sum(axis=(0, 1)) /
                          ((beta_old) ** 2).sum(axis=(0, 1)))
        if verbose:
            print(f"Min Divergence {epsilon[nc].min()}")

        log_like = np.array([
            (( X @ beta[:, :, m])  * Y).sum() - kappa0(X @ beta[:, :, m] ).sum()
            for m in range(len(lambdas))
        ])
        log_like_ratio = log_like - np.array([
            ((X @ beta_old[:, :, m])  * Y).sum() - kappa0(X @ beta_old[:, :, m] ).sum()
            for m in range(len(lambdas))
        ])
        if verbose:
            print(f"Min Loglike ratio {log_like_ratio[nc].min()}")

        # Remove NaNs
        nan_stop = np.where(np.isnan(epsilon))[0]
        if len(nan_stop) > 0:
            if verbose:
                print(f"Divergence NaN comps {nan_stop}")
            for m in nc:
                beta[:, :, m] = beta_old[:, :, m]
                Z[:, :, m] = Z_old[:, :, m]
                W[:, :, m] = W_old[:, :, m]
            nc = [m for m in nc if m not in nan_stop]

        # Stop conditions
        stop_mask = (epsilon < tol) | (log_like_ratio < tol)
        nc_stop = np.where(stop_mask)[0]
        it_stop[nc_stop] = counter
        nc = [m for m in nc if m not in nc_stop]

        # Max iterations
        if counter == maxit:
            if verbose:
                print("Maximum iteration, no convergence...")
            it_stop[it_stop == 0] = counter
            break
        else:
            counter += 1

    # Transform back to original variables
    beta_tilde = np.zeros_like(beta)
    for m in range(len(lambdas)):
        if intercept:
            beta_tilde[0, 0, m] = (
                sd_Y[0] * beta[0, 0, m] +
                mu_Y[0] -
                sd_Y[0] * ((mu_X / sd_X) @ beta[1:, 0, m])
            )
            beta_tilde[1:, 0, m] = sd_Y[0] * beta[1:, 0, m] / sd_X
        else:
            beta_tilde[:, 0, m] = sd_Y[0] * beta[:, 0, m] / sd_X

    return {
        "BETA": beta_tilde,
        "beta": beta,
        "Z": Z,
        "W": W,
        "R": R,
        "it": it_stop
    }


