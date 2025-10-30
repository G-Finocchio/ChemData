#import sys
#import os

# Add ../libs to the Python path
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "libs")))

# Now you can import directly
from .Updates import *


class SHIM_3way_CB:
    def __init__(self, X, y, beta_init, gamma_init, delta_init, 
                 l1=36, l2=3, l3=4):


        self.l1, self.l2, self.l3 = l1, l2, l3
        
        range_main, range_theta, range_psi=get_ranges3(l1,l2,l3)
        
        if beta_init is None:
            beta_init = np.zeros(len(range_main))
        if gamma_init is None:
            gamma_init = np.zeros(len(range_theta))
        if delta_init is None:
            delta_init=np.zeros(len(range_psi)) 
         
        # store copies to avoid aliasing with caller
        self.beta_hat  = np.array(beta_init,  dtype=float).copy()
        self.gamma_hat = np.array(gamma_init, dtype=float).copy()
        self.delta_hat = np.array(delta_init, dtype=float).copy()

        self.means_X = np.mean(X, axis=0)
        self.stds_X  = np.std(X, axis=0, ddof=1)
        self.mean_y  = float(np.mean(y))

        self.intercept = 0.0
        self.beta_all  = None

    def fit(self, X, y,
            lambda_beta, lambda_gamma, lambda_delta,
            w_beta=1, w_gamma=1, w_delta=1, 
            tol=5e-3, max_iter=10, compute_Q=Q_bern,
            intercept=0.0, use_intercept=True):
        #Accept none as init for w
        if w_beta is None:
            w_beta=1
        if w_gamma is None:
            w_gamma=1
        if w_delta is None:
            w_delta=1
        # local working copies
        beta_hat  = self.beta_hat.copy()
        gamma_hat = self.gamma_hat.copy()
        delta_hat = self.delta_hat.copy()


        # ensure the intercept is a plain Python float (break possible 0-d np arrays)
        intercept_loc = float(intercept)

        # build initial beta_all
        beta_2way = get_beta_vec_2way3(beta=beta_hat,  l1=self.l1, l2=self.l2, l3=self.l3, 
                                       gamma=gamma_hat, only_beta=False)
        beta_3way = get_beta_vec_3way3(beta_2way=beta_2way, l1=self.l1, l2=self.l2, l3=self.l3,
                                       delta=delta_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way, beta_3way])
        # compute a real Q_old for a meaningful baseline
        Q_old = compute_Q(X=X, y=y, beta=beta_hat,
                          gamma_vec=gamma_hat, delta_vec=delta_hat, 
                          lambda_beta=lambda_beta, lambda_gamma=lambda_gamma,
                          lambda_delta=lambda_delta,
                          w_beta=w_beta, w_gamma=w_gamma, w_delta=w_delta,
                          l1=self.l1, l2=self.l2, l3=self.l3,
                          intercept=intercept_loc)
        #print("beta_hat: ",beta_hat)
        for it in range(1, max_iter + 1):
            print(f" Start iteration: {it};")

            if use_intercept:
                #print("intercept before update:", intercept_loc)
                intercept_loc = float(
                    update_intercept(X=X, y=y, beta_all=beta_all, intercept_old=intercept_loc)
                )
                #print("intercept after  update:", intercept_loc)

            # -- Update beta
            beta_hat = update_beta(X=X, y=y, beta_hat=beta_hat, gamma_hat=gamma_hat,
                                delta_hat=delta_hat,  lambda_beta=lambda_beta,
                                intercept=intercept_loc, l1=self.l1, l2=self.l2, l3=self.l3, w=w_beta)

            # -- Update gamma, delta
            gamma_hat = update_gamma(X=X, y=y, beta_hat=beta_hat, gamma_hat=gamma_hat,
                                    delta_hat=delta_hat,  lambda_gamma=lambda_gamma,
                                    l1=self.l1, l2=self.l2, l3=self.l3, 
                                    intercept=intercept_loc)

            delta_hat = update_delta(X=X, y=y, beta_hat=beta_hat, gamma_hat=gamma_hat,
                                    delta_hat=delta_hat, lambda_delta=lambda_delta,
                                    l1=self.l1, l2=self.l2, l3=self.l3, 
                                    intercept=intercept_loc)


            # rebuild beta_all
            beta_2way = get_beta_vec_2way3(beta=beta_hat,  l1=self.l1, l2=self.l2, l3=self.l3, 
                                           gamma=gamma_hat, only_beta=False)
            beta_3way = get_beta_vec_3way3(beta_2way=beta_2way, l1=self.l1, l2=self.l2, l3=self.l3, 
                                           delta=delta_hat, only_beta=False)
            beta_all = np.concatenate([beta_hat, beta_2way, beta_3way])

            # evaluate objective
            Q_new = compute_Q(X=X, y=y, beta=beta_hat,
                              gamma_vec=gamma_hat, delta_vec=delta_hat,
                              lambda_beta=lambda_beta, lambda_gamma=lambda_gamma,
                              lambda_delta=lambda_delta,
                              w_beta=w_beta, w_gamma=w_gamma, w_delta=w_delta,
                              l1=self.l1, l2=self.l2, l3=self.l3, 
                              intercept=intercept_loc)

            if Q_new == Q_old:
                self.beta_hat, self.gamma_hat, self.delta_hat = beta_hat, gamma_hat, delta_hat
                self.beta_all = beta_all
                self.intercept = intercept_loc
                print("The model has converged.")
                return {"beta_hat": self.beta_hat, "gamma_hat": self.gamma_hat, "delta_hat": self.delta_hat, "beta_all": beta_all, "intercept": self.intercept}

            rel_dif = compute_rel_dif(Q_old=Q_old, Q_new=Q_new)
            #print("  Relative difference is ", rel_dif)
            #print("Q_new old: ", Q_new, Q_old)

            if abs(rel_dif) <= tol:
                self.beta_hat, self.gamma_hat, self.delta_hat = beta_hat, gamma_hat, delta_hat
                self.beta_all = beta_all
                self.intercept = intercept_loc
                print("The model has converged.")
                #print("Last self.intercept", self.intercept, intercept_loc)
                return {"beta_hat": self.beta_hat, "gamma_hat": self.gamma_hat, "delta_hat": self.delta_hat, "beta_all": beta_all, "intercept": self.intercept}

            Q_old = Q_new  # advance baseline

        print("It has not converged. The max number of iterations can be increased.")
        self.beta_hat, self.gamma_hat, self.delta_hat = beta_hat, gamma_hat, delta_hat
        self.beta_all = beta_all
        self.intercept = intercept_loc
        return {"beta_hat": self.beta_hat, "gamma_hat": self.gamma_hat, "delta_hat": self.delta_hat, "beta_all": beta_all, "intercept": self.intercept}

    def predict(self, X_new):
        beta_all = self.beta_all
        v = X_new @ beta_all + self.intercept
        y_pred = kappa1(v)
        return y_pred

    def R2_score(self, X_new, y_true, verbose=True):
        y_pred = self.predict(X_new)
        if verbose:
            print("R2 score is", r2(y_true, y_pred))
        return r2(y_true, y_pred)
