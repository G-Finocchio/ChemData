from .Functions_for_updates import *
from .Helpers import *
def update_delta(X, y, beta_hat, gamma_hat, delta_hat,  lambda_delta,
                l1=36, l2=3, l3=4,  intercept=0):
    beta_2way = get_beta_vec_2way3(
        beta=beta_hat, l1=l1, l2=l2, l3=l3,
        gamma=gamma_hat, only_beta=False
    )

    beta_3way = get_beta_vec_3way3(
        beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, 
        delta=delta_hat, only_beta=True
    )



    X_3way = X[:, get_ranges3(l1, l2, l3)[2]]

    if np.var(beta_3way) == 0:
        print("Three ways are 0.")
        return beta_3way * 0

    X_tilde = np.tile(beta_3way, (X_3way.shape[0], 1)) * X_3way
    nonzero_cols = np.where(np.any(X_tilde != 0, axis=0))[0]
    X_tilde_nonzero=X_tilde[:, nonzero_cols]

    X_c = X[:, (
        get_ranges3(l1, l2, l3)[0]
        + get_ranges3(l1, l2, l3)[1]
    )]

    beta_c = np.concatenate((beta_hat, beta_2way))
    C = (X_c @ beta_c + intercept).reshape(-1,1)
    delta_hat_old = delta_hat.copy()
    
    Q_old = Q_bern(
        X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
            lambda_beta=0, lambda_gamma=0, lambda_delta=lambda_delta,
            w_beta=1, w_gamma=1, w_delta=1, l1=l1, l2=l2, l3=l3, intercept=intercept
    )
    delta_hat = 0 * delta_hat
    lasso_rez = irlasso_cb(
        X=X_tilde_nonzero, Y=y, lambdas=[lambda_delta], w_lambda=None,
        centering=False, scaling=False, intercept=False,
        maxit=10, tol=0.5, sd_tol=1e-6, verbose=False, C=C, beta0=None
    )
    delta_hat[nonzero_cols]= np.array(lasso_rez['beta'][:, 0, 0])


    Q_new = Q_bern(
        X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
        lambda_beta=0, lambda_gamma=0, lambda_delta=lambda_delta,
        w_beta=1, w_gamma=1, w_delta=1, 
        l1=l1, l2=l2, l3=l3, intercept=intercept
    )

    if Q_new - Q_old > abs(Q_old) * 1e-1:
        delta_hat = delta_hat_old
        print("There might be numerical instability in update delta "
              "which was taken care of by using old delta.")

    if Q_new - Q_old >= 0:
        delta_hat = delta_hat_old
       

    print("Updated delta")
    return delta_hat


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

def update_gamma(X, y, beta_hat, gamma_hat, delta_hat,  lambda_gamma,
                l1=36, l2=3, l3=4, w=1, intercept=0):
    range1 = np.arange(0, l1 )
    range2 = np.arange(l1 , l1 + l2 )
    range3 = np.arange(l1 + l2 , l1 + l2 + l3 )

    ranges = get_ranges3(l1, l2, l3)
    X_main = X[:, np.ravel(ranges[0])]
    X_2way = X[:, np.ravel(ranges[1])]
    X_3way = X[:, np.ravel(ranges[2])]

    w = np.ones(len(gamma_hat)) if np.isscalar(w) or len(np.atleast_1d(w)) <= 1 else w


    #print('i','j')
    for i in range1:
        for j in range2:
            pos_ij = matrix_position_to_vector_index_2way3([i, j], l1, l2, l3)
            if beta_hat[i] == 0 or beta_hat[j] == 0:
                #print("direct 0")
                gamma_hat[pos_ij] = 0
            else:
                
                discard_from_c_2way = [pos_ij]
                discard_from_c_3way = []
                beta_2way = get_beta_vec_2way3(beta_hat, l1, l2, l3, gamma_hat, only_beta=False)
                beta_3way = get_beta_vec_3way3(beta_2way=beta_2way, delta=delta_hat, l1=l1, l2=l2, l3=l3, only_beta=False)
                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way]).reshape(-1, 1)
                eta = X @ beta_all + intercept
                
                # Contribution of 2-way interactions without gamma
                two_ways = X_2way[:, pos_ij] * beta_hat[i  ] * beta_hat[j  ]
                #print("toways: ", two_ways)

                # Contribution of 3-way interactions without gamma
                three_ways = 0
                for k in range3:
                    idx_ijk = table_position_to_vector_index3([i, j, k], l1, l2, l3)
                    three_ways += (
                        X_3way[:, idx_ijk]
                        * (beta_hat[i ] * beta_hat[j ] * beta_hat[k ]) ** 2
                        * gamma_hat[matrix_position_to_vector_index_2way3([i, k], l1, l2, l3)]
                        * gamma_hat[matrix_position_to_vector_index_2way3([j, k], l1, l2, l3)]
                        * delta_hat[idx_ijk]
                    )
                    discard_from_c_3way.append(idx_ijk)



                X_tilde = two_ways + three_ways

                C = (
                    intercept
                    + X_main @ beta_hat
                    + X_2way[:, [k for k in range(X_2way.shape[1]) if k not in discard_from_c_2way]] @ beta_2way[[k for k in range(len(beta_2way)) if k not in discard_from_c_2way]]
                    + X_3way[:, [k for k in range(X_3way.shape[1]) if k not in discard_from_c_3way]] @ beta_3way[[k for k in range(len(beta_3way)) if k not in discard_from_c_3way]]
                )

                Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
                                lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, 
                                w_beta=1, w_gamma=1, w_delta=1, 
                                l1=l1, l2=l2, l3=l3,  intercept=intercept)

                gamma_old_value = gamma_hat[pos_ij].copy()
                new_val=minimizer_Q_bern_gamma(X=X_tilde, y=y, C=C,
                                                        lambda_ = lambda_gamma, beta_old = gamma_old_value,
                                                        weight=1, scaled=True)

                gamma_hat[pos_ij] = new_val


                # Recompute beta vectors after update
                beta_2way = get_beta_vec_2way3(beta_hat, l1, l2, l3, gamma_hat, only_beta=False)
                beta_3way = get_beta_vec_3way3(beta_2way, delta_hat, l1, l2, l3,  only_beta=False)

                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way]).reshape(-1, 1)
                eta = X @ beta_all + intercept

                Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
                                lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0,
                                w_beta=1, w_gamma=1, w_delta=1,
                                l1=l1, l2=l2, l3=l3,  intercept=intercept)

                if Q_new - Q_old > abs(Q_old) * 1e-2:
                    print("There might be numerical instability in update gamma, which was taken care of.")

                if Q_new - Q_old >= 0:
                    gamma_hat[pos_ij] = gamma_old_value
                    
    #print('i','k')               
    for i in range1:
        for k in range3:

            if beta_hat[i] == 0 or beta_hat[k] == 0:
                gamma_hat[matrix_position_to_vector_index_2way3([i, k], l1, l2, l3)] = 0
                continue

            discard_from_c_2way = [matrix_position_to_vector_index_2way3([i, k], l1, l2, l3)]
            discard_from_c_3way = []

            beta_2way = get_beta_vec_2way3(beta_hat, l1, l2, l3, gamma_hat, only_beta=False)
            beta_3way = get_beta_vec_3way3(beta_2way, delta_hat, l1, l2, l3, only_beta=False)
            beta_all = np.concatenate([beta_hat, beta_2way, beta_3way]).reshape(-1, 1)

            # Contribution of 2-ways without gamma
            two_ways = X_2way[:, matrix_position_to_vector_index_2way3([i, k], l1, l2, l3)] * beta_hat[i] * beta_hat[k]

            # Contribution of 3-ways without gamma
            three_ways = 0
            for j in range2:
                idx = table_position_to_vector_index3([i, j, k], l1, l2, l3)
                three_ways += (X_3way[:, idx] *
                               (beta_hat[i] * beta_hat[j] * beta_hat[k])**2 *
                               gamma_hat[matrix_position_to_vector_index_2way3([i, j], l1, l2, l3)] *
                               gamma_hat[matrix_position_to_vector_index_2way3([j, k], l1, l2, l3)] *delta_hat[idx])
                discard_from_c_3way.append(idx)


            X_tilde = two_ways + three_ways


            C = (intercept + X_main @ beta_hat +
                X_2way[:, [i for i in range(X_2way.shape[1]) if i not in discard_from_c_2way]] @
                beta_2way[[i for i in range(len(beta_2way)) if i not in discard_from_c_2way]] +
                X_3way[:, [i for i in range(X_3way.shape[1]) if i not in discard_from_c_3way]] @
                beta_3way[[i for i in range(len(beta_3way)) if i not in discard_from_c_3way]] )

            Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
                        lambda_beta=0, lambda_gamma=lambda_gamma,
                        lambda_delta=0,
                        w_beta=1, w_gamma=1, w_delta=1, 
                        l1=l1, l2=l2, l3=l3,
                        intercept=intercept)

            gamma_idx = matrix_position_to_vector_index_2way3([i, k], l1, l2, l3)
            gamma_old_value = gamma_hat[gamma_idx]

            gamma_hat[gamma_idx] =minimizer_Q_bern_gamma(X=X_tilde, y=y, C=C,
                                                        lambda_ = lambda_gamma, beta_old = gamma_old_value,
                                                        weight=1, scaled=True)
            # Recompute beta vectors after gamma update
            beta_2way = get_beta_vec_2way3(beta_hat, l1, l2, l3,  gamma_hat, only_beta=False)
            beta_3way = get_beta_vec_3way3(beta_2way, delta_hat, l1, l2, l3, only_beta=False)

            Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                            lambda_beta=0, lambda_gamma=lambda_gamma,
                            lambda_delta=0, 
                            w_beta=1, w_gamma=1, w_delta=1,
                            l1=l1, l2=l2, l3=l3, 
                            intercept=intercept)

            if Q_new - Q_old > abs(Q_old) * 1e-2:
                print("There might be numerical instability in update gamma, which was taken care of.")

            if Q_new - Q_old >= 0:
                gamma_hat[gamma_idx] = gamma_old_value

   
   
                    
    # j, k
    #print('j,k')
    for j in range2:
        for k in range3:
            if beta_hat[j] == 0 or beta_hat[k] == 0:
                gamma_hat[matrix_position_to_vector_index_2way3([j, k], l1=l1, l2=l2, l3=l3)] = 0
            else:
                discard_from_c_2way = [matrix_position_to_vector_index_2way3([j, k], l1=l1, l2=l2, l3=l3)]
                discard_from_c_3way = []

                beta_2way = get_beta_vec_2way3(
                    beta=beta_hat,
                    l1=l1,
                    l2=l2,
                    l3=l3,
                    gamma=gamma_hat,
                    only_beta=False
                )

                beta_3way = get_beta_vec_3way3(
                    beta_2way=beta_2way,
                    l1=l1,
                    l2=l2,
                    l3=l3,
                    delta=delta_hat,
                    only_beta=False
                )



                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way]).reshape(-1, 1)
                beta_all = beta_all.reshape(-1, 1)
                eta = X @ beta_all + intercept

                # Contribution of 2-ways without gamma
                two_ways = X_2way[:, matrix_position_to_vector_index_2way3([j, k], l1=l1, l2=l2, l3=l3)] * beta_hat[j] * beta_hat[k]

                # Contribution of 3-ways without gamma
                three_ways = 0
                for i in range1:
                    # Compute 3 ways ijk
                    three_ways += X_3way[:, table_position_to_vector_index3([i, j, k], l1=l1, l2=l2, l3=l3)] * \
                        (beta_hat[i] * beta_hat[j] * beta_hat[k])**2 * \
                        gamma_hat[matrix_position_to_vector_index_2way3([i, j], l1=l1, l2=l2, l3=l3)] * \
                        gamma_hat[matrix_position_to_vector_index_2way3([i, k], l1=l1, l2=l2, l3=l3)] * \
                        delta_hat[table_position_to_vector_index3([i, j, k], l1=l1, l2=l2, l3=l3)]
                    discard_from_c_3way.append(table_position_to_vector_index3([i, j, k], l1=l1, l2=l2, l3=l3))


                X_tilde = two_ways + three_ways

                C = intercept + X_main @ beta_hat + X_2way[:, np.setdiff1d(range(X_2way.shape[1]), discard_from_c_2way)] @ beta_2way[np.setdiff1d(range(beta_2way.shape[0]), discard_from_c_2way)] + \
                    X_3way[:, np.setdiff1d(range(X_3way.shape[1]), discard_from_c_3way)] @ beta_3way[np.setdiff1d(range(beta_3way.shape[0]), discard_from_c_3way)] 


                Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, 
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, w_beta =1, w_gamma=1,
                            w_delta=1,  l1=l1, l2=l2, l3=l3,  scaled=True, intercept=intercept)

                # Use the minimizer
                gamma_old_value = gamma_hat[matrix_position_to_vector_index_2way3([j, k], l1=l1, l2=l2, l3=l3)]
                gamma_hat[matrix_position_to_vector_index_2way3([j, k], l1=l1, l2=l2, l3=l3)] = minimizer_Q_bern_gamma(X=X_tilde, y=y, C=C,
                                                        lambda_ = lambda_gamma, beta_old = gamma_old_value, weight=1, scaled=True)

                # Update betas
                beta_2way = get_beta_vec_2way3(beta_hat, l1, l2, l3,  gamma_hat, only_beta=False)
                beta_3way = get_beta_vec_3way3(beta_2way, delta_hat, l1, l2, l3, only_beta=False)
                
                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way]).reshape(-1, 1)
                eta = X @ beta_all + intercept
                Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, 
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, w_beta =1, w_gamma=1,
                            w_delta=1,  l1=l1, l2=l2, l3=l3, scaled=True, intercept=intercept)

                if Q_new - Q_old > abs(Q_old) * 1e-2:
                    print("There might be numerical instability in update gamma.")

                #  Keep the new estimator only if better than old one
                if Q_new - Q_old >= 0:
                    gamma_hat[matrix_position_to_vector_index_2way3([j, k], l1, l2, l3)] = gamma_old_value


    print("Updated gamma")
    return gamma_hat




############# Update the vector of main effects  ##################### 


def update_beta(X, y, beta_hat, gamma_hat, delta_hat, lambda_beta,
                l1=36, l2=3, l3=4, w=1, intercept=0):
    
    range1 = np.arange(l1)
    range2 = np.arange(l1, l1+l2)
    range3 = np.arange(l1+l2, l1+l2+l3)

    
    X_main, X_2way, X_3way = (
        X[:, get_ranges3(l1=l1, l2=l2, l3=l3)[i]] for i in range(3)
    )
    
    w = np.ones(len(beta_hat)) if np.isscalar(w) or len(np.atleast_1d(w)) <= 1 else w
    beta_2way = get_beta_vec_2way3(beta=beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta=False)
    beta_3way = get_beta_vec_3way3(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3,  delta=delta_hat, only_beta=False)
    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way]).reshape(-1, 1)
    eta = X @ beta_all + intercept
    #print("i")
    for i in range1:
        discard_from_c_main = [i]
        mains = X_main[:, i]
        #print('Mains: ', mains)
        discard_from_c_2way, two_ways = [], 0
        for j in np.concatenate([range2, range3]):
            pos = matrix_position_to_vector_index_2way3([i, j], l1=l1, l2=l2, l3=l3)
            discard_from_c_2way.append(pos)
            two_ways += X_2way[:, pos] * beta_hat[j] * gamma_hat[pos]
        #print("2ways: ", two_ways)
        discard_from_c_3way, three_ways = [], 0
        for j in range2:
            for k in range3:
                idx = table_position_to_vector_index3([i, j, k], l1=l1, l2=l2, l3=l3)
                three_ways += X_3way[:, idx] * (beta_hat[j] * beta_hat[k])**2 * \
                              gamma_hat[matrix_position_to_vector_index_2way3([i, k], l1=l1, l2=l2, l3=l3)] * \
                              gamma_hat[matrix_position_to_vector_index_2way3([j, k], l1=l1, l2=l2, l3=l3)] * \
                              gamma_hat[matrix_position_to_vector_index_2way3([i, j], l1=l1, l2=l2, l3=l3)] * delta_hat[idx]
                discard_from_c_3way.append(idx)
        
        
        
        
        X_tilde, Z_tilde = mains + two_ways, three_ways
        C = intercept + X_main[:, np.arange(len(beta_hat)) != i] @ beta_hat[np.arange(len(beta_hat)) != i] + \
            X_2way[:, [x for x in range(X_2way.shape[1]) if x not in discard_from_c_2way]] @ beta_2way[[x for x in range(len(beta_2way)) if x not in discard_from_c_2way]] + \
            X_3way[:, [x for x in range(X_3way.shape[1]) if x not in discard_from_c_3way]] @ beta_3way[[x for x in range(len(beta_3way)) if x not in discard_from_c_3way]] 

        
        Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
                        w_beta=w, w_gamma=1, w_delta=1, l1=l1, l2=l2, l3=l3, intercept=intercept)
        
        beta_old_value = beta_hat[i]
        beta_hat[i] = minimizer_Q_bern_beta(X=X_tilde, Z=Z_tilde, y=y, C=C,
                                            lambda_=lambda_beta, beta_old=beta_old_value, weight=w[i], scaled=True)
        
        beta_2way = get_beta_vec_2way3(beta=beta_hat, l1=l1, l2=l2, l3=l3, gamma=gamma_hat, only_beta=False)
        beta_3way = get_beta_vec_3way3(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, delta=delta_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way, beta_3way]).reshape(-1, 1)
        eta = X @ beta_all + intercept
        
        Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0,
                        w_beta=w, w_gamma=1, w_delta=1, l1=l1, l2=l2, l3=l3,  intercept=intercept)
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
        for ik in list(range1) + list(range3) :
            small_idx = min(ik, j)
            big_idx = max(ik, j)

            pos = matrix_position_to_vector_index_2way3(
                [small_idx, big_idx], l1=l1, l2=l2, l3=l3
            )
            
            discard_from_c_2way.append(pos)
            two_ways = two_ways + X_2way[:, pos] * beta_hat[ik] * gamma_hat[pos]

        discard_from_c_3way = []
        # Contribution of 3-ways without beta
        three_ways = 0
        for i in range1:
            for k in list(range3) :
                pos = table_position_to_vector_index3([i, j, k], l1=l1, l2=l2, l3=l3)
                three_ways = three_ways + X_3way[:, pos] * (beta_hat[i] * beta_hat[k])**2 * \
                    gamma_hat[matrix_position_to_vector_index_2way3([i, k], l1=l1, l2=l2, l3=l3)] * \
                    gamma_hat[matrix_position_to_vector_index_2way3([j, k], l1=l1, l2=l2, l3=l3)] * \
                    gamma_hat[matrix_position_to_vector_index_2way3([i, j], l1=l1, l2=l2, l3=l3)] * \
                    delta_hat[pos]
                discard_from_c_3way.append(pos)


        X_tilde = mains + two_ways
        Z_tilde = three_ways

        C = intercept + X_main[:, np.arange(X_main.shape[1]) != j] @ beta_hat[np.arange(len(beta_hat)) != j] + \
            X_2way[:, [i for i in range(X_2way.shape[1]) if i not in discard_from_c_2way]] @ beta_2way[[i for i in range(len(beta_2way)) if i not in discard_from_c_2way]] + \
            X_3way[:, [i for i in range(X_3way.shape[1]) if i not in discard_from_c_3way]] @ beta_3way[[i for i in range(len(beta_3way)) if i not in discard_from_c_3way]] 

        Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, 
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0,
                        w_beta=w, w_gamma=1, w_delta=1,
                        l1=l1, l2=l2, l3=l3,  intercept=intercept)

        beta_old_value = beta_hat[j]
        # Use the minimizer
        beta_hat[j] = minimizer_Q_bern_beta(
            X=X_tilde, Z=Z_tilde, y=y, C=C,
            lambda_=lambda_beta, beta_old=beta_old_value, weight=w[j], scaled=True
        )

        beta_2way = get_beta_vec_2way3(beta=beta_hat, l1=l1, l2=l2, l3=l3,  gamma=gamma_hat, only_beta=False)
        beta_3way = get_beta_vec_3way3(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3,  delta=delta_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way, beta_3way]).reshape(-1, 1)
        eta = X @ beta_all + intercept

        Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0,
                        w_beta=w, w_gamma=1, w_delta=1, 
                        l1=l1, l2=l2, l3=l3,  intercept=intercept)

        if Q_new - Q_old > abs(Q_old) * 1e-2:
            print("There might be numerical instability in update beta which was taken care of.")

        # Keep the new estimator only if it is better than the old one
        if Q_new - Q_old >= 0:
            beta_hat[j] = beta_old_value

    #print("k")
    for k in range3:
        discard_from_c_main = [k]
        mains = X_main[:, k]

        discard_from_c_2way = []
        # Contribution of 2-ways without beta
        two_ways = 0
        for il in list(range1) + list(range2) :
            small_idx = min(il, k)
            big_idx = max(il, k)

            discard_from_c_2way.append(
                matrix_position_to_vector_index_2way3(
                    [small_idx, big_idx], l1=l1, l2=l2, l3=l3
                )
            )

            two_ways += (
                X_2way[:, matrix_position_to_vector_index_2way3([small_idx, big_idx],
                                                            l1=l1, l2=l2, l3=l3)]
                * beta_hat[il]
                * gamma_hat[matrix_position_to_vector_index_2way3([small_idx, big_idx],
                                                                l1=l1, l2=l2, l3=l3)]
            )

        discard_from_c_3way = []
        # Contribution of 3-ways without beta
        three_ways = 0

        for i in range1:
            for j in range2:
                three_ways += (
                    X_3way[:, table_position_to_vector_index3([i, j, k],
                                                                l1=l1, l2=l2, l3=l3)]
                    * (beta_hat[i] * beta_hat[j]) ** 2
                    * gamma_hat[matrix_position_to_vector_index_2way3([i, j],
                                                                    l1=l1, l2=l2, l3=l3)]
                    * gamma_hat[matrix_position_to_vector_index_2way3([i, k],
                                                                    l1=l1, l2=l2, l3=l3)]
                    * gamma_hat[matrix_position_to_vector_index_2way3([j, k],
                                                                    l1=l1, l2=l2, l3=l3)]
                    * delta_hat[table_position_to_vector_index3([i, j, k],
                                                                    l1=l1, l2=l2, l3=l3)]
                )
                discard_from_c_3way.append(
                    table_position_to_vector_index3([i, j, k],
                                                        l1=l1, l2=l2, l3=l3)
                )

        
        X_tilde = mains + two_ways
        Z_tilde = three_ways
        C = (
            intercept
            + X_main[:, np.arange(X_main.shape[1]) != k] @ beta_hat[np.arange(len(beta_hat)) != k]
            + X_2way[:, [idx for idx in range(X_2way.shape[1]) if idx not in discard_from_c_2way]] @ beta_2way[[idx for idx in range(len(beta_2way)) if idx not in discard_from_c_2way]]
            + X_3way[:, [idx for idx in range(X_3way.shape[1]) if idx not in discard_from_c_3way]] @ beta_3way[[idx for idx in range(len(beta_3way)) if idx not in discard_from_c_3way]]
        )  # C is the offset term

        Q_old = Q_bern(
            X=X, y=y, beta=beta_hat,
            gamma_vec=gamma_hat, delta_vec=delta_hat,
            lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
            w_beta=w, w_gamma=1, w_delta=1,
            l1=l1, l2=l2, l3=l3, 
            intercept=intercept
        )

        beta_old_value = beta_hat[k]
        # Use the minimizer
        beta_hat[k] = minimizer_Q_bern_beta(
            X=X_tilde, Z=Z_tilde, y=y, C=C,
            lambda_=lambda_beta, beta_old=beta_old_value,
            weight=w[k], scaled=True
        )  # Intercept is in C

        beta_2way = get_beta_vec_2way3(beta=beta_hat, l1=l1, l2=l2, l3=l3, 
                                        gamma=gamma_hat, only_beta=False)
        beta_3way = get_beta_vec_3way3(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3,
                                        delta=delta_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way, beta_3way]).reshape(-1, 1)
        eta = X @ beta_all + intercept

        Q_new = Q_bern(
            X=X, y=y, beta=beta_hat,
            gamma_vec=gamma_hat, delta_vec=delta_hat, 
            lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, 
            w_beta=w, w_gamma=1, w_delta=1, 
            l1=l1, l2=l2, l3=l3, 
            intercept=intercept
        )

        if Q_new - Q_old > abs(Q_old) * 1e-2:
            print("There might be numerical instability in update beta which was taken care of.")

        # Keep the new estimator only if it is better than the old one. Ensures numerical stability.
        if Q_new - Q_old >= 0:
            beta_hat[k] = beta_old_value



    print("Updated beta")
    return beta_hat







