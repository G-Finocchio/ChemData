import numpy as np
import itertools
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score as r2, root_mean_squared_error as rmse
from .cb_lasso import *

def grid_search_lasso(
    X_train, X_test, y_train, y_test,
    lambda_grid=None,
    maxit=20,
    tol=0.5,
    sd_tol=1e-6,
    val_size=0.3,
    random_state=42,
    plot=True
):
    """
    Perform grid search for irlasso_cb with an internal validation split.
    Select best lambda by validation R², then evaluate on the test set.

    Returns
    -------
    best_coefs : np.array
        Coefficients for best model (excluding intercept)
    best_intercept : float
        Intercept term
    best_lambda : float
        Lambda value giving best validation R²
    test_r2 : float
        Final R² on the test set
    results : list of (lambda, train_r2, val_r2)
    """

    # 1. Split train into (train, val)
    X_tr, X_val, y_tr, y_val = train_test_split(
        X_train, y_train, test_size=val_size, random_state=random_state
    )

    results = []
    best_val_r2 = -np.inf
    best_lambda = None
    best_coefs = None
    best_intercept = None

    # 2. Iterate over lambda grid
    for lmbd in lambda_grid:
        print(f"\nTrying lambda: {lmbd}")

        res_lasso_tr = irlasso_cb(
            X=X_tr, Y=y_tr,
            lambdas=[lmbd], w_lambda=None, beta0=None,
            centering=False, scaling=False, intercept=True,
            maxit=maxit, tol=tol, sd_tol=sd_tol, verbose=False
        )

        coefs = np.array(res_lasso_tr["beta"][1:, 0, 0])  # skip intercept
        intercept = res_lasso_tr["beta"][0, 0, 0]

        # Predictions
        y_pred_tr = kappa1(X_tr @ coefs + intercept)
        y_pred_val = kappa1(X_val @ coefs + intercept)

        # Compute R²
        train_r2 = r2(y_tr, y_pred_tr)
        val_r2 = r2(y_val, y_pred_val)
        results.append((lmbd, train_r2, val_r2))

        print(f"Train R² = {train_r2:.4f}, Val R² = {val_r2:.4f}")

        if val_r2 > best_val_r2:
            best_val_r2 = val_r2
            best_lambda = lmbd
            best_coefs = coefs
            best_intercept = intercept

    # 3. Retrain on full training data with best lambda
    print(f"\nBest λ on validation: {best_lambda} (Val R² = {best_val_r2:.4f})")

    res_best = irlasso_cb(
        X=X_train, Y=y_train,
        lambdas=[best_lambda], w_lambda=None, beta0=None,
        centering=False, scaling=False, intercept=True,
        maxit=maxit, tol=tol, sd_tol=sd_tol, verbose=False
    )

    best_coefs = np.array(res_best["beta"][1:, 0, 0])
    best_intercept = res_best["beta"][0, 0, 0]

    # 4. Evaluate on test set
    y_pred_test = kappa1(X_test @ best_coefs + best_intercept)
    test_r2 = r2(y_test, y_pred_test)
    print(f"\nFinal Test R² = {test_r2:.4f}")
    test_rmse=rmse(y_test, y_pred_test)
    print(f"\n Final Test RMSE: {test_rmse:.2f}")

    # Optional plot
    if plot:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6,6))
        plt.scatter(y_pred_test, y_test, alpha=0.6, edgecolors='k')
        # add red diagonal
        min_val = min(y_test.min(), y_pred_test.min())
        max_val = max(y_test.max(), y_pred_test.max())
        plt.plot([min_val, max_val], [min_val, max_val], 'r--', lw=2)
        plt.xlabel("Predicted values")
        plt.ylabel("True values")
        plt.title(f"True vs Predicted (R² = {test_r2:.2f})")
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.show()

    return best_coefs, best_intercept, best_lambda, test_r2, results


