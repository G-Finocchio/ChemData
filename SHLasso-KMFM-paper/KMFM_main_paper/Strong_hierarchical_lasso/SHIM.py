import importlib
from sklearn.metrics import r2_score, root_mean_squared_error as rmse
import matplotlib.pyplot as plt

import importlib

import itertools
import numpy as np
from sklearn.model_selection import train_test_split
from .Helpers_SHIM import *

class SHIM:
    def __init__(self, X, y, init_list, levels_list, lambdas_list, 
                 w_list=None, tol=5e-3, max_iter=3, use_init_main=True, use_init_main_w=True, lmd_lasso=1e-5, transform_to_identifiable=False,response='cb' ):
        self.X = X
        self.y = y
        self.init_list = init_list if init_list is not None else [None]*len(lambdas_list)
        self.levels_list = levels_list
        self.lambdas_list = lambdas_list
        self.w_list = [None] * len(lambdas_list) if w_list is None else w_list
        self.response=response.lower()
        if len(self.levels_list)==5:
            self.init_list=[]
            print("For the 5 way model initial values are not used/needed.")
        self.tol = tol
        self.max_iter = max_iter
        self.use_init_main = use_init_main
        self.use_init_main_w = use_init_main_w
        self.lmd_lasso=lmd_lasso
        if transform_to_identifiable ==True:
            self.X=get_X_all_identifiable(X, levels_list)

        if self.response=="cb":
            # Build model right away
            self.model = self._build_model_cb()
        elif  self.response=="normal":
            self.model = self._build_model_normal()
        else:
            raise ValueError(f"Unknown response_type: {response}")


    def _build_model_cb(self):
        """Internal: select the correct model class and instantiate it."""
        n_factors = len(self.levels_list)
        if n_factors not in (2, 3, 4, 5):
            raise ValueError("Only 2, 3, 4, or 5 factors are supported.")

        # --- Validation ---
        if n_factors == 2 and len(self.init_list) != 2:
            raise ValueError("2-way models require init_list of length 2.")
        elif n_factors == 3 and len(self.init_list) != 3:
            raise ValueError("3-way models require init_list of length 3.")
        elif n_factors == 4 and len(self.init_list) != 4:
            raise ValueError("4-way models require init_list of length 4.")
        elif n_factors == 5 and self.init_list:
            raise ValueError("5-way model does not use init_list, pass [] for init list!")

        expected_lambdas_len = {
            2: 2,
            3: 3,
            4: 4,
            5: 4,
        }[n_factors]
        if len(self.lambdas_list) != expected_lambdas_len:
            raise ValueError(f"{n_factors}-way model requires lambdas_list of length {expected_lambdas_len}.")

        # --- Import correct module ---
        model_mapping = {
            2: ("SHIM_2way_CB", "SHIM_2way_CB"),
            3: ("SHIM_3way_CB", "SHIM_3way_CB"),
            4: ("SHIM_4way_CB", "SHIM_4way_CB"),
            5: ("SHIM_5way_CB", "SHLasso_5way_CB_fast"),
        }
        module_name, class_name = model_mapping[n_factors]
        package_name = f".libs_{n_factors}way_cb"
        module = importlib.import_module(f"{package_name}.{module_name}", package="Strong_hierarchical_lasso")
        model_class = getattr(module, class_name)

        # --- kwargs ---
        levels_kwargs = {f"l{i+1}": self.levels_list[i] for i in range(n_factors)}
        init_kwargs = {}

        if n_factors in (2, 3, 4):
            if self.use_init_main:
                # run quick irlasso_cb on main effects only
                fct_module = importlib.import_module(f"{package_name}.Functions_for_updates", package="Strong_hierarchical_lasso")

                irlasso_cb = getattr(fct_module, "irlasso_cb")
                kappa0 = getattr(fct_module, "kappa0")
                kappa1 = getattr(fct_module, "kappa1")
                kappa2 = getattr(fct_module, "kappa2")

                res_lasso_tr = irlasso_cb(
                    X=self.X[:, :sum(self.levels_list)],
                    Y=self.y,
                    lambdas=[self.lmd_lasso],
                    w_lambda=None,
                    beta0=None,
                    centering=False,
                    scaling=False,
                    intercept=True,
                    maxit=10,
                    tol=0.05,
                    sd_tol=1e-6,
                    verbose=False
                )
                beta_init = np.array(res_lasso_tr["beta"][1:, 0, 0])  # includes intercept
                self.init_list[0] = beta_init  # overwrite main init
                if self.use_init_main_w:
                    self.w_list[0] = np.minimum(1e7, np.abs(1 / beta_init))

            # assign init kwargs
            if n_factors == 2:
                names = ["beta_init", "gamma_init"]
            elif n_factors == 3:
                names = ["beta_init", "gamma_init", "delta_init"]
            else:  # n_factors == 4
                names = ["beta_init", "gamma_init", "delta_init", "tau_init"]
            init_kwargs = {names[i]: self.init_list[i] for i in range(len(names))}
        if self.response =='normal':
            init_kwargs={} # normal resp does not use any init

        # --- instantiate ---
        return model_class(self.X, self.y, **init_kwargs, **levels_kwargs)
    
    
    
    def _build_model_normal(self):
        """Internal: select the correct model class and instantiate it."""
        n_factors = len(self.levels_list)
        if n_factors not in (2, 3, 4, 5):
            raise ValueError("Only 2, 3, 4, or 5 factors are supported.")

        # --- Validation ---
        if n_factors == 2 and len(self.init_list) != 2:
            raise ValueError("2-way models require init_list of length 2.")
        elif n_factors == 3 and len(self.init_list) != 3:
            raise ValueError("3-way models require init_list of length 3.")
        elif n_factors == 4 and len(self.init_list) != 4:
            raise ValueError("4-way models require init_list of length 4.")
        elif n_factors == 5 and self.init_list:
            raise ValueError("5-way model does not use init_list, pass [] for init list!")

        expected_lambdas_len = {
            2: 2,
            3: 3,
            4: 4,
            5: 4,
        }[n_factors]
        if len(self.lambdas_list) != expected_lambdas_len:
            raise ValueError(f"{n_factors}-way model requires lambdas_list of length {expected_lambdas_len}.")

        # --- Import correct module ---
        model_mapping = {
            2: ("SHIM_normal", "SHIM_2way_normal"),
            3: ("SHIM_normal", "SHIM_3way_normal"),
            4: ("SHIM_normal", "SHIM_4way_normal"),
            5: ("SHIM_normal", "SHIM_5way_normal"),
        }
        module_name, class_name = model_mapping[n_factors]
        package_name = f".libs_normal"
        module = importlib.import_module(f"{package_name}.{module_name}", package="Strong_hierarchical_lasso")
        model_class = getattr(module, class_name)

        # --- kwargs ---
        levels_kwargs = {f"l{i+1}": self.levels_list[i] for i in range(n_factors)}
        init_kwargs = {}

        # --- instantiate ---
        return model_class(self.X, self.y, **init_kwargs, **levels_kwargs)


    def fit_model(self):
        """Fit the model with the provided lambda values."""
        n_factors = len(self.levels_list)

        # --- Lambdas ---
        # Define lambda parameter names dynamically based on number of factors
        lambda_names = ["lambda_beta", "lambda_gamma"]
        if n_factors >= 3:
            lambda_names.append("lambda_delta")
        if n_factors >= 4:
            lambda_names.append("lambda_tau")

        lmd_kwargs = {lambda_names[i]: self.lambdas_list[i] for i in range(len(self.lambdas_list))}

        # --- Weights ---
        # Define weight parameter names dynamically 
        w_names = ["w_beta", "w_gamma"]
        if n_factors >= 3:
            w_names.append("w_delta")
        if n_factors >= 4:
            w_names.append("w_tau")

        w_kwargs = {w_names[i]: self.w_list[i] for i in range(len(self.w_list))}
        
        if self.response == 'cb':
            # --- Fit the model ---
            self.model.fit(
                self.X,self.y,
                **lmd_kwargs,
                **w_kwargs,
                tol=self.tol,
                max_iter=self.max_iter)
        else:
            self.model.fit(
                self.X,self.y,
                **lmd_kwargs)
            

        return self


    
    def predict(self, X_new, **kwargs):
        """Make predictions using the fitted model's own predict method."""
        if self.model is None:
            raise RuntimeError("Model is not fitted yet. Call fit_model() first.")
        return self.model.predict(X_new, **kwargs)

    def R2_score(self, X_new, y_true, scale=False, verbose=True, plot=True):
        """Compute R² score on new data, with optional true-vs-pred plot."""
       # Call predict safely
        if scale is not None and hasattr(self.model, "scale") and getattr(self.model, "scale"):
            y_pred = self.predict(X_new, scale=scale)
        else:
            y_pred = self.predict(X_new)
        score = r2_score(y_true, y_pred)

        if verbose:
            print("R² score:", score)

        if plot:
            plt.figure(figsize=(6,6))
            plt.scatter(y_pred, y_true, alpha=0.6, edgecolors='k')
            # add red diagonal
            min_val = min(y_true.min(), y_pred.min())
            max_val = max(y_true.max(), y_pred.max())
            plt.plot([min_val, max_val], [min_val, max_val], 'r--', lw=2)
            plt.xlabel("Predicted values")
            plt.ylabel("True values")
            plt.title(f"True vs Predicted (R² = {score:.2f})")
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.show()

        return score
    
    
    def rmse_score(self, X_new, y_true, scale=False):
        """Compute R² score on new data, with optional true-vs-pred plot."""
       # Call predict safely
        if scale is not None and hasattr(self.model, "scale") and getattr(self.model, "scale"):
            y_pred = self.predict(X_new, scale=scale)
        else:
            y_pred = self.predict(X_new)
        score = rmse(y_true, y_pred)
        print("RMSE:", score)
        return score




def grid_search_shim(
    X_train,X_test, y_train, y_test,
    levels_list,
    init_list=None,
    lambda_grid=None,
    max_iter=10,
    use_init_main=False,
    use_init_main_w=False,
    val_size=0.3,
    random_state=42,
    plot=True,
    w_beta=None,
    tol=5e-3,
    response='cb'
):
    """
    Perform grid search with an internal validation split, select best model
    by validation R², then evaluate on the test set.

    Returns
    -------
    best_model : SHIM object (fitted on train+val with best lambdas)
    best_config : tuple of best lambda values
    test_r2 : float, final R² on test set
    results : list of (lambdas, train_r2, val_r2)
    """

    # 1. Split train into (train, val)
    X_tr, X_val, y_tr, y_val = train_test_split(
        X_train, y_train, test_size=val_size, random_state=random_state
    )

    results = []
    best_val_r2 = -np.inf
    best_config = None
    best_model = None

    # 2. Loop over lambda grid
    for lambdas in itertools.product(*lambda_grid):
        print(f"\nTrying lambdas: {lambdas}")

        shim = SHIM(
            X=X_tr,
            y=y_tr,
            init_list=init_list,
            levels_list=levels_list,
            lambdas_list=list(lambdas),
            max_iter=max_iter,
            use_init_main=use_init_main,
            use_init_main_w=use_init_main_w,
            w_list=[w_beta, None, None],
            tol=tol,
            response=response
        )

        shim_fitted = shim.fit_model()
        train_r2 = shim_fitted.R2_score(X_tr, y_tr, plot=True, verbose=False)
        val_r2 = shim_fitted.R2_score(X_val, y_val, plot=True, verbose=False)

        results.append((lambdas, train_r2, val_r2))
        print(f"Train R² = {train_r2:.2f}, Val R² = {val_r2:.2f}")

        if val_r2 > best_val_r2:
            best_val_r2 = val_r2
            best_config = lambdas

    # 3. Retrain on full train set with best config
    print(f"\nBest config on validation: {best_config} (Val R² = {best_val_r2:.2f})")
    best_model = SHIM(
        X=X_train,
        y=y_train,
        init_list=init_list,
        levels_list=levels_list,
        lambdas_list=list(best_config),
        max_iter=max_iter,
        use_init_main=use_init_main,
        use_init_main_w=use_init_main_w,
        w_list=[w_beta, None, None], 
        tol=tol, response=response
    ).fit_model()

    # 4. Final evaluation on train set
    train_r2 = best_model.R2_score(X_train, y_train, plot=plot, verbose=True)
    print(f"\nFinal Train R² = {train_r2:.2f}")

    # 4. Final evaluation on test set
    test_r2 = best_model.R2_score(X_test, y_test, plot=plot, verbose=True)
    print(f"\nFinal Test R² = {test_r2:.2f}")
    
    test_rmse= best_model.rmse_score(X_test, y_test)
    print(f"\n Final Test RMSE: {test_rmse:.2f}")

    return best_model, best_config, test_r2, results