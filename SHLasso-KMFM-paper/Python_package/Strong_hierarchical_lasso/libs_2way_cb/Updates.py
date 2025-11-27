from .Functions_for_updates import *
from .Helpers import *


####################### UPDATE INTERCEPT ######################


def update_intercept(X, y, beta_all, intercept_old):
    # function to find the minimum for intercept
    def fctint(intercept):
        v = X @ beta_all + intercept
        log_like = np.sum(y * v - kappa0(v))
        loss = -log_like  # scale does not matter
        return loss

    interval = (
        min(-intercept_old / 2 - 7, 5 * intercept_old / 2 - 7),
        max(-intercept_old / 2 + 7, 5 * intercept_old / 2 + 7)
    )

    res = minimize_scalar(fctint, bounds=interval, method="bounded")
    minimum = res.x

    f_0 = fctint(0)
    if f_0 <= fctint(minimum) and f_0 <= fctint(intercept_old):
        return 0

    if fctint(intercept_old) <= fctint(minimum):
        return intercept_old

    print("Updated intercept")
    return minimum





###################################### UPDATE GAMMA ###########################################

def update_gamma(X, y, beta_hat, gamma_hat,  lambda_gamma,
                l1=7, l2=8,  w=1, intercept=0):
    range1 = np.arange(0, l1 )
    range2 = np.arange(l1 , l1 + l2 )

    ranges = get_ranges2(l1, l2)
    X_main = X[:, np.ravel(ranges[0])]
    X_2way = X[:, np.ravel(ranges[1])]

    w = np.ones(len(gamma_hat)) if np.isscalar(w) or len(np.atleast_1d(w)) <= 1 else w


    #print('i','j')
    for i in range1:
        for j in range2:
            pos_ij = matrix_position_to_vector_index_2way2([i, j], l1, l2)
            if beta_hat[i] == 0 or beta_hat[j] == 0:
                #print("direct 0")
                gamma_hat[pos_ij] = 0
            else:
                
                discard_from_c_2way = [pos_ij]
                beta_2way = get_beta_vec_2way2(beta_hat, l1, l2, gamma_hat, only_beta=False)
                beta_all = np.concatenate([beta_hat, beta_2way]).reshape(-1, 1)
                eta = X @ beta_all + intercept
                
                # Contribution of 2-way interactions without gamma
                two_ways = X_2way[:, pos_ij] * beta_hat[i] * beta_hat[j]
                #print("toways: ", two_ways)


                X_tilde = two_ways 

                C = (
                    intercept
                    + X_main @ beta_hat
                    + X_2way[:, [k for k in range(X_2way.shape[1]) if k not in discard_from_c_2way]] @ beta_2way[[k for k in range(len(beta_2way)) if k not in discard_from_c_2way]]
                )

                Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, 
                                lambda_beta=0, lambda_gamma=lambda_gamma, 
                                w_beta=1, w_gamma=1,  
                                l1=l1, l2=l2,  intercept=intercept)

                gamma_old_value = gamma_hat[pos_ij].copy()
                new_val=minimizer_Q_bern_gamma(X=X_tilde, y=y, C=C,
                                                        lambda_ = lambda_gamma, beta_old = gamma_old_value,
                                                        weight=1, scaled=True)

                gamma_hat[pos_ij] = new_val


                # Recompute beta vectors after update
                beta_2way = get_beta_vec_2way2(beta_hat, l1, l2, gamma_hat, only_beta=False)

                beta_all = np.concatenate([beta_hat, beta_2way]).reshape(-1, 1)
                eta = X @ beta_all + intercept

                Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, 
                                lambda_beta=0, lambda_gamma=lambda_gamma,
                                w_beta=1, w_gamma=1, 
                                l1=l1, l2=l2, intercept=intercept)

                if Q_new - Q_old > abs(Q_old) * 1e-2:
                    print("There might be numerical instability in update gamma, which was taken care of.")

                if Q_new - Q_old >= 0:
                    gamma_hat[pos_ij] = gamma_old_value
                    


    print("Updated gamma")
    return gamma_hat




############# Update the vector of main effects  ##################### 


def update_beta(X, y, beta_hat, gamma_hat, lambda_beta,
                l1=7, l2=8, w=1, intercept=0):
    
    range1 = np.arange(l1)
    range2 = np.arange(l1, l1+l2)


    
    X_main, X_2way = (
        X[:, get_ranges2(l1=l1, l2=l2)[i]] for i in range(2)
    )
    
    w = np.ones(len(beta_hat)) if np.isscalar(w) or len(np.atleast_1d(w)) <= 1 else w
    beta_2way = get_beta_vec_2way2(beta=beta_hat, l1=l1, l2=l2, gamma=gamma_hat, only_beta=False)
    beta_all = np.concatenate([beta_hat, beta_2way]).reshape(-1, 1)
    eta = X @ beta_all + intercept
    #print("i")
    for i in range1:
        discard_from_c_main = [i]
        mains = X_main[:, i]
        #print('Mains: ', mains)
        discard_from_c_2way, two_ways = [], 0
        for j in range2:
            pos = matrix_position_to_vector_index_2way2([i, j], l1=l1, l2=l2)
            discard_from_c_2way.append(pos)
            two_ways += X_2way[:, pos] * beta_hat[j] * gamma_hat[pos]
        
        
        
        
        X_tilde= mains + two_ways
        C = intercept + X_main[:, np.arange(len(beta_hat)) != i] @ beta_hat[np.arange(len(beta_hat)) != i] + \
            X_2way[:, [x for x in range(X_2way.shape[1]) if x not in discard_from_c_2way]] @ beta_2way[[x for x in range(len(beta_2way)) if x not in discard_from_c_2way]] 
        
        Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat,  
                        lambda_beta=lambda_beta, lambda_gamma=0,  
                        w_beta=w, w_gamma=1, l1=l1, l2=l2,  intercept=intercept)
        
        beta_old_value = beta_hat[i]
        beta_hat[i] = minimizer_Q_bern_beta(X=X_tilde,  y=y, C=C,
                                            lambda_=lambda_beta, beta_old=beta_old_value, weight=w[i], scaled=True)
        
        beta_2way = get_beta_vec_2way2(beta=beta_hat, l1=l1, l2=l2, gamma=gamma_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way]).reshape(-1, 1)
        eta = X @ beta_all + intercept
        
        Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, 
                        lambda_beta=lambda_beta, lambda_gamma=0, 
                        w_beta=w, w_gamma=1,  l1=l1, l2=l2,   intercept=intercept)
        if Q_new - Q_old > abs(Q_old) * 1e-2:
            #print(Q_new, Q_old)
            print("There might be numerical instability in update beta which was taken care of.")
        if Q_new - Q_old >= 0:
            beta_hat[i] = beta_old_value

            
            
            
    # j
    #print("j")
    for j in range2:
        discard_from_c_main = [j]
        mains = np.array(X_main[:, j])

        discard_from_c_2way = []
        # Contribution of 2-ways without beta
        two_ways = 0
        for ik in range1 :
            small_idx = min(ik, j)
            big_idx = max(ik, j)

            pos = matrix_position_to_vector_index_2way2(
                [small_idx, big_idx], l1=l1, l2=l2
            )
            
            discard_from_c_2way.append(pos)
            two_ways = two_ways + X_2way[:, pos] * beta_hat[ik] * gamma_hat[pos]



        X_tilde = mains + two_ways

        C = intercept + X_main[:, np.arange(X_main.shape[1]) != j] @ beta_hat[np.arange(len(beta_hat)) != j] + \
            X_2way[:, [i for i in range(X_2way.shape[1]) if i not in discard_from_c_2way]] @ beta_2way[[i for i in range(len(beta_2way)) if i not in discard_from_c_2way]]  

        Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat,
                        lambda_beta=lambda_beta, lambda_gamma=0, 
                        w_beta=w, w_gamma=1,
                        l1=l1, l2=l2,  intercept=intercept)

        beta_old_value = beta_hat[j]
        # Use the minimizer
        beta_hat[j] = minimizer_Q_bern_beta(
            X=X_tilde, y=y, C=C,
            lambda_=lambda_beta, beta_old=beta_old_value, weight=w[j], scaled=True
        )

        beta_2way = get_beta_vec_2way2(beta=beta_hat, l1=l1, l2=l2,   gamma=gamma_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way]).reshape(-1, 1)
        eta = X @ beta_all + intercept

        Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat,
                        lambda_beta=lambda_beta, lambda_gamma=0, 
                        w_beta=w, w_gamma=1,
                        l1=l1, l2=l2,  intercept=intercept)

        if Q_new - Q_old > abs(Q_old) * 1e-2:
            print("There might be numerical instability in update beta which was taken care of.")

        # Keep the new estimator only if it is better than the old one
        if Q_new - Q_old >= 0:
            beta_hat[j] = beta_old_value





    print("Updated beta")
    return beta_hat







