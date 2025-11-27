from .Functions_for_updates import *
from .Helpers import *
def update_tau(X, y, beta_hat, gamma_hat, delta_hat, tau_hat, lambda_tau,
                l1=21, l2=14, l3=2, l4=3, intercept=0):
    beta_2way = get_beta_vec_2way4(
        beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4,
        gamma=gamma_hat, only_beta=False
    )

    beta_3way = get_beta_vec_3way4(
        beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4,
        delta=delta_hat, only_beta=False
    )

    beta_4way = get_beta_vec_4way4(
        beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4,
        tau=tau_hat, only_beta=True
    )

    X_4way = X[:, get_ranges4(l1, l2, l3, l4)[3]]

    if np.var(beta_4way) == 0:
        print("Four ways are 0.")
        return beta_4way * 0

    X_tilde = np.tile(beta_4way, (X_4way.shape[0], 1)) * X_4way

    X_c = X[:, (
        get_ranges4(l1, l2, l3, l4)[0]
        + get_ranges4(l1, l2, l3, l4)[1]
        + get_ranges4(l1, l2, l3, l4)[2]
    )]

    beta_c = np.concatenate((beta_hat, beta_2way, beta_3way))
    C = (X_c @ beta_c + intercept).reshape(-1,1)
    tau_hat_old = tau_hat.copy()
    Q_old = Q_bern(
        X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
        tau_vec=tau_hat, lambda_beta=0, lambda_gamma=0, lambda_delta=0,
        lambda_tau=lambda_tau, w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
        l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept
    )

    lasso_rez = irlasso_cb(
        X=X_tilde, Y=y, lambdas=[lambda_tau], w_lambda=None,
        centering=False, scaling=False, intercept=True,
        maxit=10, tol=0.0545, sd_tol=1e-6, verbose=False, C=C
    )
    tau_hat= np.array(lasso_rez['beta'][1:, 0, 0])


    Q_new = Q_bern(
        X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
        tau_vec=tau_hat, lambda_beta=0, lambda_gamma=0, lambda_delta=0,
        lambda_tau=lambda_tau, w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
        l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept
    )

    if Q_new - Q_old > abs(Q_old) * 1e-5:
        tau_hat = tau_hat_old
        print("There might be numerical instability in update tau "
              "which was taken care of by using old tau.")

    if Q_new - Q_old >= 0:
        tau_hat = tau_hat_old

    print("Updated tau")
    return tau_hat


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

def update_gamma(X, y, beta_hat, gamma_hat, delta_hat, tau_hat, lambda_gamma,
                l1=21, l2=14, l3=2, l4=3, w=1, intercept=0):
    range1 = np.arange(0, l1 )
    range2 = np.arange(l1 , l1 + l2 )
    range3 = np.arange(l1 + l2 , l1 + l2 + l3 )
    range4 = np.arange(l1 + l2 + l3 , l1 + l2 + l3 + l4 )

    ranges = get_ranges4(l1, l2, l3, l4)
    X_main = X[:, np.ravel(ranges[0])]
    X_2way = X[:, np.ravel(ranges[1])]
    X_3way = X[:, np.ravel(ranges[2])]
    X_4way = X[:, np.ravel(ranges[3])]

    w = np.ones(len(gamma_hat)) if np.isscalar(w) or len(np.atleast_1d(w)) <= 1 else w


    #print('i','j')
    for i in range1:
        ##print('a')
        for j in range2:
            ##print('b')
            pos_ij = get_position_vec_from_theta_matrix4([i, j], l1, l2, l3, l4)
            #print(pos_ij)
            if beta_hat[i] == 0 or beta_hat[j] == 0:
                #print("direct 0")
                gamma_hat[pos_ij] = 0
            else:
                ##print('c0')
                discard_from_c_2way = [pos_ij]
                discard_from_c_3way = []
                #print(delta_hat)
                #print(l1,l2,l3,l4)
                #print('c00')
                beta_2way = get_beta_vec_2way4(beta_hat, l1, l2, l3, l4, gamma_hat, only_beta=False)
                #print(beta_2way)
                beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, delta=delta_hat, l1=l1, l2=l2, l3=l3, l4=l4, only_beta=False)
                beta_4way = get_beta_vec_4way4(beta_3way, tau_hat, l1, l2, l3, l4, only_beta=False)
                #print('c1')
                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                eta = X @ beta_all + intercept
                
                #print("c2")
                #print(X_main[:, [0,4]])
                #print(X_main.shape)
                #print(X_2way[:,pos_ij])
                # Contribution of 2-way interactions without gamma
                two_ways = X_2way[:, pos_ij] * beta_hat[i  ] * beta_hat[j  ]
                #print("toways: ", two_ways)

                # Contribution of 3-way interactions without gamma
                three_ways = 0
                #print('c3')
                for k in np.concatenate([range3, range4]):
                    #print('c')
                    idx_ijk = psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)
                    three_ways += (
                        X_3way[:, idx_ijk]
                        * (beta_hat[i ] * beta_hat[j ] * beta_hat[k ]) ** 2
                        * gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)]
                        * gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)]
                        * delta_hat[idx_ijk]
                    )
                    discard_from_c_3way.append(idx_ijk)

                # Contribution of 4-way interactions without gamma
                four_ways = 0
                discard_from_c_4way = []
                for k in range3:
                    #print('d')
                    for l in range4:
                        #print('e')
                        idx_ijkl = phi_table_position_to_vector_index4([i, j, k, l], l1, l2, l3, l4)
                        term = (
                            X_4way[:, idx_ijkl]
                            * (beta_hat[i ] * beta_hat[j ]) ** 6
                            * (gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)]
                               * gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)]
                               * gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)]
                               * gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1, l2, l3, l4)]
                               * gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1, l2, l3, l4)]
                              ) ** 2
                            * delta_hat[psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)]
                            * delta_hat[psi_table_position_to_vector_index4([i, j, l], l1, l2, l3, l4)]
                            * delta_hat[psi_table_position_to_vector_index4([i, k, l], l1, l2, l3, l4)]
                            * delta_hat[psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)]
                            * tau_hat[idx_ijkl]
                            * (beta_hat[k ] * beta_hat[l ]) ** 6
                        )
                        four_ways += term
                        discard_from_c_4way.append(idx_ijkl)

                X_tilde = two_ways + three_ways
                Z_tilde = four_ways

                C = (
                    intercept
                    + X_main @ beta_hat
                    + X_2way[:, [k for k in range(X_2way.shape[1]) if k not in discard_from_c_2way]] @ beta_2way[[k for k in range(len(beta_2way)) if k not in discard_from_c_2way]]
                    + X_3way[:, [k for k in range(X_3way.shape[1]) if k not in discard_from_c_3way]] @ beta_3way[[k for k in range(len(beta_3way)) if k not in discard_from_c_3way]]
                    + X_4way[:, [k for k in range(X_4way.shape[1]) if k not in discard_from_c_4way]] @ beta_4way[[k for k in range(len(beta_4way)) if k not in discard_from_c_4way]]
                )

                Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,tau_vec=tau_hat,
                                lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, lambda_tau=0,
                                w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                                l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)

                gamma_old_value = gamma_hat[pos_ij].copy()
                new_val=minimizer_Q_bern_gamma(X=X_tilde, Z=Z_tilde, y=y, C=C,
                                                        lambda_ = lambda_gamma, beta_old = gamma_old_value,
                                                        weight=1, scaled=True)
                #print("new_val: ",new_val)
                gamma_hat[pos_ij] = new_val
                #print("gammhatpos ij: ", gamma_hat[pos_ij])
                #print(gamma_hat)

                # Recompute beta vectors after update
                beta_2way = get_beta_vec_2way4(beta_hat, l1, l2, l3, l4, gamma_hat, only_beta=False)
                beta_3way = get_beta_vec_3way4(beta_2way, delta_hat, l1, l2, l3, l4, only_beta=False)
                beta_4way = get_beta_vec_4way4(beta_3way, tau_hat, l1, l2, l3, l4, only_beta=False)

                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                eta = X @ beta_all + intercept

                Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,tau_vec=tau_hat,
                                lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta=0, lambda_tau=0,
                                w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                                l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)

                if Q_new - Q_old > abs(Q_old) * 1e-2:
                    print("There might be numerical instability in update gamma, which was taken care of.")

                if Q_new - Q_old >= 0:
                    gamma_hat[pos_ij] = gamma_old_value
                    
    #print('i','k')               
    for i in range1:
        for k in range3:

            if beta_hat[i] == 0 or beta_hat[k] == 0:
                gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)] = 0
                continue

            discard_from_c_2way = [get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)]
            discard_from_c_3way = []

            beta_2way = get_beta_vec_2way4(beta_hat, l1, l2, l3, l4, gamma_hat, only_beta=False)
            beta_3way = get_beta_vec_3way4(beta_2way, delta_hat, l1, l2, l3, l4, only_beta=False)
            beta_4way = get_beta_vec_4way4(beta_3way, tau_hat, l1, l2, l3, l4, only_beta=False)
            beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)

            # Contribution of 2-ways without gamma
            two_ways = X_2way[:, get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)] * beta_hat[i] * beta_hat[k]

            # Contribution of 3-ways without gamma
            three_ways = 0
            for j in range2:
                idx = psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)
                three_ways += (X_3way[:, idx] *
                               (beta_hat[i] * beta_hat[j] * beta_hat[k])**2 *
                               gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1, l2, l3, l4)] *
                               gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)] *delta_hat[idx])
                discard_from_c_3way.append(idx)

            for l in range4:
                idx = psi_table_position_to_vector_index4([i, k, l], l1, l2, l3, l4)
                three_ways += (X_3way[:, idx] *
                               (beta_hat[i] * beta_hat[k] * beta_hat[l])**2 *
                               gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)] *
                               gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1, l2, l3, l4)] *delta_hat[idx])
                discard_from_c_3way.append(idx)

            # Contribution of 4-ways without gamma
            four_ways = 0
            discard_from_c_4way = []
            for j in range2:
                for l in range4:
                    idx_phi = phi_table_position_to_vector_index4([i, j, k, l], l1, l2, l3, l4)
                    four_ways += (X_4way[:, idx_phi] *
                                  (beta_hat[i] * beta_hat[j])**6 *
                                  (gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1, l2, l3, l4)] *
                                   gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)] *
                                   gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)] *
                                   gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1, l2, l3, l4)] *
                                   gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1, l2, l3, l4)])**2 *
                                  delta_hat[psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)] *
                                  delta_hat[psi_table_position_to_vector_index4([i, j, l], l1, l2, l3, l4)] *
                                  delta_hat[psi_table_position_to_vector_index4([i, k, l], l1, l2, l3, l4)] *
                                  delta_hat[psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)] *
                                  tau_hat[idx_phi] *
                                  (beta_hat[k] * beta_hat[l])**6)
                    discard_from_c_4way.append(idx_phi)

            X_tilde = two_ways + three_ways
            Z_tilde = four_ways

            C = (intercept + X_main @ beta_hat +
                X_2way[:, [i for i in range(X_2way.shape[1]) if i not in discard_from_c_2way]] @
                beta_2way[[i for i in range(len(beta_2way)) if i not in discard_from_c_2way]] +
                X_3way[:, [i for i in range(X_3way.shape[1]) if i not in discard_from_c_3way]] @
                beta_3way[[i for i in range(len(beta_3way)) if i not in discard_from_c_3way]] +
                X_4way[:, [i for i in range(X_4way.shape[1]) if i not in discard_from_c_4way]] @
                beta_4way[[i for i in range(len(beta_4way)) if i not in discard_from_c_4way]])

            Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                        lambda_beta=0, lambda_gamma=lambda_gamma,
                        lambda_delta=0, lambda_tau=0,
                        w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                        l1=l1, l2=l2, l3=l3, l4=l4,
                        intercept=intercept)

            gamma_idx = get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)
            gamma_old_value = gamma_hat[gamma_idx]

            gamma_hat[gamma_idx] =minimizer_Q_bern_gamma(X=X_tilde, Z=Z_tilde, y=y, C=C,
                                                        lambda_ = lambda_gamma, beta_old = gamma_old_value,
                                                        weight=1, scaled=True)
            # Recompute beta vectors after gamma update
            beta_2way = get_beta_vec_2way4(beta_hat, l1, l2, l3, l4, gamma_hat, only_beta=False)
            beta_3way = get_beta_vec_3way4(beta_2way, delta_hat, l1, l2, l3, l4, only_beta=False)
            beta_4way = get_beta_vec_4way4(beta_3way, tau_hat, l1, l2, l3, l4, only_beta=False)

            Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                            lambda_beta=0, lambda_gamma=lambda_gamma,
                            lambda_delta=0, lambda_tau=0,
                            w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                            l1=l1, l2=l2, l3=l3, l4=l4,
                            intercept=intercept)

            if Q_new - Q_old > abs(Q_old) * 1e-2:
                print("There might be numerical instability in update gamma, which was taken care of.")

            if Q_new - Q_old >= 0:
                gamma_hat[gamma_idx] = gamma_old_value

   
   
    # i, l
    #print('i','l')
    for i in range1:
        for l in range4:
            #print(range4)
            if beta_hat[i] == 0 or beta_hat[l] == 0:
                gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)] = 0
                continue

            # Indices to discard for 2-way, 3-way, 4-way contributions
            discard_from_c_2way = [get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)]
            discard_from_c_3way = []
            discard_from_c_4way = []

            # Compute beta vectors for 2-way, 3-way, 4-way contributions
            beta_2way = get_beta_vec_2way4(beta_hat, l1, l2, l3, l4, gamma_hat, only_beta=False)
            beta_3way = get_beta_vec_3way4(beta_2way, delta_hat, l1, l2, l3, l4, only_beta=False)
            beta_4way = get_beta_vec_4way4(beta_3way, tau_hat, l1, l2, l3, l4, only_beta=False)
            beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)

            eta = X @ beta_all + intercept

            # Contribution of 2-ways without gamma
            two_ways = X_2way[:, get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)] \
                        * beta_hat[i] * beta_hat[l]

            # Contribution of 3-ways without gamma
            three_ways = 0
            for j in range2:  # Compute 3-ways ijl
                idx_3way = psi_table_position_to_vector_index4([i, j, l], l1, l2, l3, l4)
                three_ways += (X_3way[:, idx_3way] *
                            (beta_hat[i] * beta_hat[j] * beta_hat[l]) ** 2 *
                            gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1, l2, l3, l4)] *
                            gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1, l2, l3, l4)] *
                            delta_hat[idx_3way])
                discard_from_c_3way.append(idx_3way)

            for k in range3:  # Compute 3-ways ikl
                idx_3way = psi_table_position_to_vector_index4([i, k, l], l1, l2, l3, l4)
                three_ways += (X_3way[:, idx_3way] *
                            (beta_hat[i] * beta_hat[k] * beta_hat[l]) ** 2 *
                            gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1, l2, l3, l4)] *
                            gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)] *
                            delta_hat[idx_3way])
                discard_from_c_3way.append(idx_3way)

            # Contribution of 4-ways without gamma
            four_ways = 0
            for j in range2:
                for k in range3:
                    idx_4way = phi_table_position_to_vector_index4([i, j, k, l], l1, l2, l3, l4)
                    four_ways += (X_4way[:, idx_4way] *
                                   (beta_hat[i] * beta_hat[j]) ** 6 *
                                   (gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1, l2, l3, l4)] *
                                    gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)] *
                                    gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)] *
                                    gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1, l2, l3, l4)] *
                                    gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1, l2, l3, l4)]) ** 2 *
                                   delta_hat[psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)] *
                                   delta_hat[psi_table_position_to_vector_index4([i, j, l], l1, l2, l3, l4)] *
                                   delta_hat[psi_table_position_to_vector_index4([i, k, l], l1, l2, l3, l4)] *
                                   delta_hat[psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)] *
                                   tau_hat[idx_4way] *
                                   (beta_hat[k] * beta_hat[l]) ** 6)
                    discard_from_c_4way.append(idx_4way)

            X_tilde = two_ways + three_ways
            Z_tilde = four_ways

            C = (intercept + X_main @ beta_hat +
                X_2way[:, [i for i in range(X_2way.shape[1]) if i not in discard_from_c_2way]] @ beta_2way[[i for i in range(len(beta_2way)) if i not in discard_from_c_2way]] +
                X_3way[:, [i for i in range(X_3way.shape[1]) if i not in discard_from_c_3way]] @ beta_3way[[i for i in range(len(beta_3way)) if i not in discard_from_c_3way]] +
                X_4way[:, [i for i in range(X_4way.shape[1]) if i not in discard_from_c_4way]] @ beta_4way[[i for i in range(len(beta_4way)) if i not in discard_from_c_4way]])  # C is the offset term

            Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, tau_vec=tau_hat,
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, lambda_tau = 0, w_beta =1, w_gamma=1,
                            w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, scaled=True, intercept=intercept)
            gamma_old_value = gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)]

            # Use the minimizer
            gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)] = minimizer_Q_bern_gamma(X=X_tilde, Z=Z_tilde, y=y, C=C,
                                                        lambda_ = lambda_gamma, beta_old = gamma_old_value, weight=1, scaled=True)

            # Recompute beta vectors after gamma update
            beta_2way = get_beta_vec_2way4(beta_hat, l1, l2, l3, l4, gamma_hat, only_beta=False)
            beta_3way = get_beta_vec_3way4(beta_2way, delta_hat, l1, l2, l3, l4, only_beta=False)
            beta_4way = get_beta_vec_4way4(beta_3way, tau_hat, l1, l2, l3, l4, only_beta=False)
            #print('aaaa')
            beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
            #print('bbbb')
            eta = X @ beta_all + intercept
            Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, tau_vec=tau_hat,
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, lambda_tau = 0, w_beta =1, w_gamma=1,
                            w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)

            if Q_new - Q_old > abs(Q_old) * 1e-2:
                print("There might be numerical instability in update gamma.")

            # Keep the new estimator only if better than old one
            if Q_new - Q_old >= 0:
                gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)] = gamma_old_value
                
                
    # j, k
    #print('j,k')
    for j in range2:
        for k in range3:
            if beta_hat[j] == 0 or beta_hat[k] == 0:
                gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] = 0
            else:
                discard_from_c_2way = [get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)]
                discard_from_c_3way = []

                beta_2way = get_beta_vec_2way4(
                    beta=beta_hat,
                    l1=l1,
                    l2=l2,
                    l3=l3,
                    l4=l4,
                    gamma=gamma_hat,
                    only_beta=False
                )

                beta_3way = get_beta_vec_3way4(
                    beta_2way=beta_2way,
                    l1=l1,
                    l2=l2,
                    l3=l3,
                    l4=l4,
                    delta=delta_hat,
                    only_beta=False
                )

                beta_4way = get_beta_vec_4way4(
                    beta_3way=beta_3way,
                    l1=l1,
                    l2=l2,
                    l3=l3,
                    l4=l4,
                    tau=tau_hat,
                    only_beta=False
                )

                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                beta_all = beta_all.reshape(-1, 1)
                eta = X @ beta_all + intercept

                # Contribution of 2-ways without gamma
                two_ways = X_2way[:, get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * beta_hat[j] * beta_hat[k]

                # Contribution of 3-ways without gamma
                three_ways = 0
                for i in range1:
                    # Compute 3 ways ijk
                    three_ways += X_3way[:, psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        (beta_hat[i] * beta_hat[j] * beta_hat[k])**2 * \
                        gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)]
                    discard_from_c_3way.append(psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4))

                for l in range4:
                    # Compute 3 ways ikl
                    three_ways += X_3way[:, psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        (beta_hat[j] * beta_hat[k] * beta_hat[l])**2 * \
                        gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)]
                    discard_from_c_3way.append(psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4))

                # Contribution of 4-ways without gamma
                four_ways = 0
                discard_from_c_4way = []
                for i in range1:
                    for l in range4:
                        four_ways += X_4way[:, phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                            (beta_hat[i] * beta_hat[j])**6 * \
                            (gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                             gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                             gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                             gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                             gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)])**2 * \
                            delta_hat[psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                            delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                            delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                            delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                            tau_hat[phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                            (beta_hat[k] * beta_hat[l])**6
                        discard_from_c_4way.append(phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4))

                X_tilde = two_ways + three_ways
                Z_tilde = four_ways

                C = intercept + X_main @ beta_hat + X_2way[:, np.setdiff1d(range(X_2way.shape[1]), discard_from_c_2way)] @ beta_2way[np.setdiff1d(range(beta_2way.shape[0]), discard_from_c_2way)] + \
                    X_3way[:, np.setdiff1d(range(X_3way.shape[1]), discard_from_c_3way)] @ beta_3way[np.setdiff1d(range(beta_3way.shape[0]), discard_from_c_3way)] + \
                    X_4way[:, np.setdiff1d(range(X_4way.shape[1]), discard_from_c_4way)] @ beta_4way[np.setdiff1d(range(beta_4way.shape[0]), discard_from_c_4way)]

                Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, tau_vec=tau_hat,
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, lambda_tau = 0, w_beta =1, w_gamma=1,
                            w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, scaled=True, intercept=intercept)

                # Use the minimizer
                gamma_old_value = gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)]
                gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] = minimizer_Q_bern_gamma(X=X_tilde, Z=Z_tilde, y=y, C=C,
                                                        lambda_ = lambda_gamma, beta_old = gamma_old_value, weight=1, scaled=True)

                # Update betas
                beta_2way = get_beta_vec_2way4(beta_hat, l1, l2, l3, l4, gamma_hat, only_beta=False)
                beta_3way = get_beta_vec_3way4(beta_2way, delta_hat, l1, l2, l3, l4, only_beta=False)
                beta_4way = get_beta_vec_4way4(beta_3way, tau_hat, l1, l2, l3, l4, only_beta=False)
                
                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                eta = X @ beta_all + intercept
                Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, tau_vec=tau_hat,
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, lambda_tau = 0, w_beta =1, w_gamma=1,
                            w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, scaled=True, intercept=intercept)

                if Q_new - Q_old > abs(Q_old) * 1e-2:
                    print("There might be numerical instability in update gamma.")

                #  Keep the new estimator only if better than old one
                if Q_new - Q_old >= 0:
                    gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)] = gamma_old_value



    #j, l
    #print('j,l')
    for j in range2:
        for l in range4:
            if beta_hat[j] == 0 or beta_hat[l] == 0:
                gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] = 0
            else:
                discard_from_c_2way = [
                    get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)
                ]
                discard_from_c_3way = []

                beta_2way = get_beta_vec_2way4(
                    beta=beta_hat,
                    l1=l1,
                    l2=l2,
                    l3=l3,
                    l4=l4,
                    gamma=gamma_hat,
                    only_beta=False
                )
                beta_3way = get_beta_vec_3way4(
                    beta_2way=beta_2way,
                    l1=l1,
                    l2=l2,
                    l3=l3,
                    l4=l4,
                    delta=delta_hat,
                    only_beta=False
                )
                beta_4way = get_beta_vec_4way4(
                    beta_3way=beta_3way,
                    l1=l1,
                    l2=l2,
                    l3=l3,
                    l4=l4,
                    tau=tau_hat,
                    only_beta=False
                )

                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                eta = X @ beta_all + intercept

                # Contribution of 2-ways without gamma
                two_ways = X_2way[:, get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * beta_hat[j] * beta_hat[l]

                # Contribution of 3-ways without gamma
                three_ways = 0

                for i in range1:
                    three_ways += X_3way[:, psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * (beta_hat[i] * beta_hat[j] * beta_hat[l])**2 * \
                                  gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                  gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)]
                    discard_from_c_3way.append(psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4))

                for k in range3:
                    three_ways += X_3way[:, psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * (beta_hat[j] * beta_hat[k] * beta_hat[l])**2 * \
                                  gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                  gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)]
                    discard_from_c_3way.append(psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4))

                discard_from_c_4way = []
                four_ways = 0

                for i in range1:
                    for k in range3:
                        four_ways += X_4way[:, phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * (beta_hat[i] * beta_hat[j])**6 * \
                                     (gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                      gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                      gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                      gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                      gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)])**2 * \
                                     delta_hat[psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     tau_hat[phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     (beta_hat[k] * beta_hat[l])**6
                        discard_from_c_4way.append(phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4))

                X_tilde = two_ways + three_ways
                Z_tilde = four_ways

                C = intercept + X_main @ beta_hat + \
                    X_2way[:, [i for i in range(len(beta_2way)) if i not in discard_from_c_2way]] @ beta_2way[[i for i in range(len(beta_2way)) if i not in discard_from_c_2way]] + \
                    X_3way[:, [i for i in range(len(beta_3way)) if i not in discard_from_c_3way]] @ beta_3way[[i for i in range(len(beta_3way)) if i not in discard_from_c_3way]] + \
                    X_4way[:, [i for i in range(len(beta_4way)) if i not in discard_from_c_4way]] @ beta_4way[[i for i in range(len(beta_4way)) if i not in discard_from_c_4way]]

                Q_old =Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, tau_vec=tau_hat,
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, lambda_tau = 0, w_beta =1, w_gamma=1,
                            w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, scaled=True, intercept=intercept)

                gamma_old_value = gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)]
                gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] = minimizer_Q_bern_gamma(X=X_tilde, Z=Z_tilde, y=y, C=C, 
                                                                                                                            lambda_ = lambda_gamma, beta_old = gamma_old_value, weight=1, scaled=True)

                beta_2way = get_beta_vec_2way4(beta_hat, l1, l2, l3, l4, gamma_hat, only_beta=False)
                beta_3way = get_beta_vec_3way4(beta_2way, delta_hat, l1, l2, l3, l4, only_beta=False)
                beta_4way = get_beta_vec_4way4(beta_3way, tau_hat, l1, l2, l3, l4, only_beta=False)

                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                eta = X @ beta_all + intercept

                Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, tau_vec=tau_hat,
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, lambda_tau = 0, w_beta =1, w_gamma=1,
                            w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, scaled=True, intercept=intercept)

                if Q_new - Q_old > abs(Q_old) * 1e-2:
                    print("There might be numerical instability in update gamma.")

                if Q_new - Q_old >= 0:
                    gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] = gamma_old_value



    #print('k,l')
# k, l
    for k in range3:
        for l in range4:
            if beta_hat[k] == 0 or beta_hat[l] == 0:
                gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] = 0
            else:
                discard_from_c_2way = [
                    get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                ]

                discard_from_c_3way = []

                beta_2way = get_beta_vec_2way4(
                    beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4,
                    gamma=gamma_hat, only_beta=False
                )
                beta_3way = get_beta_vec_3way4(
                    beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4,
                    delta=delta_hat, only_beta=False
                )
                beta_4way = get_beta_vec_4way4(
                    beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4,
                    tau=tau_hat, only_beta=False
                )
                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                eta = X @ beta_all + intercept

                # Contribution of 2-ways without gamma
                two_ways = X_2way[:, get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                           * beta_hat[k] * beta_hat[l]

                # Contribution of 3-ways without gamma
                three_ways = 0
                for i in range1:
                    three_ways += X_3way[:, psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                    * (beta_hat[i] * beta_hat[k] * beta_hat[l])**2 \
                                    * gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                    * gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                    * delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)]
                    discard_from_c_3way.append(psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4))

                for j in range2:
                    three_ways += X_3way[:, psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                    * (beta_hat[j] * beta_hat[k] * beta_hat[l])**2 \
                                    * gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                    * gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                    * delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)]
                    discard_from_c_3way.append(psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4))

                discard_from_c_4way = []
                # Contribution of 4-ways without gamma
                four_ways = 0
                for i in range1:
                    for j in range2:
                        four_ways += X_4way[:, phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * (beta_hat[i] * beta_hat[j])**6 \
                                     * (gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)])**2 \
                                     * delta_hat[psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * tau_hat[phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] \
                                     * (beta_hat[k] * beta_hat[l])**6
                        discard_from_c_4way.append(phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4))

                X_tilde = two_ways + three_ways
                Z_tilde = four_ways

                C = intercept + X_main @ beta_hat \
                    + X_2way[:, [i for i in range(X_2way.shape[1]) if i not in discard_from_c_2way]] @ beta_2way[[i for i in range(len(beta_2way)) if i not in discard_from_c_2way]] \
                    + X_3way[:, [i for i in range(X_3way.shape[1]) if i not in discard_from_c_3way]] @ beta_3way[[i for i in range(len(beta_3way)) if i not in discard_from_c_3way]] \
                    + X_4way[:, [i for i in range(X_4way.shape[1]) if i not in discard_from_c_4way]] @ beta_4way[[i for i in range(len(beta_4way)) if i not in discard_from_c_4way]]

                Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, tau_vec=tau_hat,
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, lambda_tau = 0, w_beta =1, w_gamma=1,
                            w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, scaled=True, intercept=intercept)

                gamma_old_value = gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)]
                # Use the minimizer
                gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] = minimizer_Q_bern_gamma(X=X_tilde, Z=Z_tilde, y=y, C=C,
                                                        lambda_ = lambda_gamma, beta_old = gamma_old_value,
                                                        weight=1, scaled=True)

                # Recompute betas after gamma update
                beta_2way = get_beta_vec_2way4(beta_hat, l1, l2, l3, l4, gamma_hat, only_beta=False)
                beta_3way = get_beta_vec_3way4(beta_2way, delta_hat, l1, l2, l3, l4, only_beta=False)
                beta_4way = get_beta_vec_4way4(beta_3way, tau_hat, l1, l2, l3, l4, only_beta=False)
                beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                eta = X @ beta_all + intercept

                Q_new =Q_bern(X=X, y=y, beta=beta_hat, gamma_vec = gamma_hat, delta_vec = delta_hat, tau_vec=tau_hat,
                            lambda_beta=0, lambda_gamma=lambda_gamma, lambda_delta= 0, lambda_tau = 0, w_beta =1, w_gamma=1,
                            w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, scaled=True, intercept=intercept)

                if Q_new - Q_old > abs(Q_old) * 1e-2:
                    print("There might be numerical instability in update gamma which was taken care of.")
                if Q_new - Q_old >= 0:
                    gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] = gamma_old_value

    print("Updated gamma")
    return gamma_hat





########################### Update for 3 ways ###############################

################## Update the vector of hierarchical constants for the three-way interactions ##################
def update_delta(X, y, beta_hat, gamma_hat, delta_hat, tau_hat, lambda_delta,
                l1=21, l2=14, l3=2, l4=3, w=1, intercept=0):

    range1 = range(0, l1)
    range2 = range(l1, l1 + l2)
    range3 = range(l1 + l2, l1 + l2 + l3)
    range4 = range(l1 + l2 + l3, l1 + l2 + l3 + l4)

    ranges = get_ranges4(l1=l1, l2=l2, l3=l3, l4=l4)
    X_main = X[:, ranges[0]]
    X_2way = X[:, ranges[1]]
    X_3way = X[:, ranges[2]]
    X_4way = X[:, ranges[3]]

    w = np.ones(len(delta_hat)) if np.isscalar(w) or len(np.atleast_1d(w)) <= 1 else w

    
    for i in range1:
        for j in range2:
            for k in range3:
                if (gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1, l2, l3, l4)] == 0 or
                    gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)] == 0 or
                    gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)] == 0):

                    delta_hat[psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)] = 0
                else:
                    discard_from_c_3way = []

                    beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    gamma=gamma_hat, only_beta=False)
                    beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    delta=delta_hat, only_beta=False)
                    beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    tau=tau_hat, only_beta=False)
                    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                    eta = X @ beta_all + intercept

                    # Contribution of 3-ways without delta
                    three_ways = (X_3way[:, psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)] *
                                (beta_hat[i] * beta_hat[j] * beta_hat[k]) ** 2 *
                                gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)] *
                                gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)] *
                                gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1, l2, l3, l4)])

                    discard_from_c_3way.append(
                        psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)
                    )

                    assert len(discard_from_c_3way) == 1, "should discard only one three way in update delta"

                    # Contribution of 4-ways without delta
                    discard_from_c_4way = []
                    four_ways = 0
                    for l in range4:
                        four_ways += (X_4way[:, phi_table_position_to_vector_index4([i, j, k, l], l1, l2, l3, l4)] *
                                      (beta_hat[i] * beta_hat[j]) ** 6 *
                                      (gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1, l2, l3, l4)]) ** 2 *
                                      delta_hat[psi_table_position_to_vector_index4([i, j, l], l1, l2, l3, l4)] *
                                      delta_hat[psi_table_position_to_vector_index4([i, k, l], l1, l2, l3, l4)] *
                                      delta_hat[psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)] *
                                      tau_hat[phi_table_position_to_vector_index4([i, j, k, l], l1, l2, l3, l4)] *
                                      (beta_hat[k] * beta_hat[l]) ** 6)

                        discard_from_c_4way.append(
                            phi_table_position_to_vector_index4([i, j, k, l], l1, l2, l3, l4)
                        )

                    X_tilde = three_ways + four_ways
                    C = (X_main @ beta_hat +
                        X_2way @ beta_2way +
                        np.delete(X_3way, discard_from_c_3way, axis=1) @ np.delete(beta_3way, discard_from_c_3way) +
                        np.delete(X_4way, discard_from_c_4way, axis=1) @ np.delete(beta_4way, discard_from_c_4way) +
                        intercept)

                    Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                                    lambda_beta=0, lambda_gamma=0, lambda_delta=lambda_delta, lambda_tau=0,
                                    w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                                    l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)

                    delta_old_value = delta_hat[psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)]
                    delta_hat[psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)] = minimizer_Q_bern_delta(
                        X=X_tilde, y=y, C=C, lambda_=lambda_delta,
                        beta_old=delta_old_value, weight=1, scaled=True)

                    beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    gamma=gamma_hat, only_beta=False)
                    beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    delta=delta_hat, only_beta=False)
                    beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    tau=tau_hat, only_beta=False)
                    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                    eta = X @ beta_all + intercept

                    Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                                    lambda_beta=0, lambda_gamma=0, lambda_delta=lambda_delta, lambda_tau=0,
                                    w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                                    l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)

                    if Q_new - Q_old > abs(Q_old) * 1e-2:
                        print("There might be numerical instability in update delta which was taken care of.")

                    if Q_new - Q_old >= 0:
                        delta_hat[psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)] = delta_old_value



    # i, j, l loops
    for i in range1:
        for j in range2:
            for l in range4:

                if (gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)] == 0 or
                    gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] == 0 or
                    gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] == 0):

                    delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] = 0

                else:
                    discard_from_c_3way = []

                    beta_2way = get_beta_vec_2way4(
                        beta=beta_hat,
                        l1=l1, l2=l2, l3=l3, l4=l4,
                        gamma=gamma_hat,
                        only_beta=False
                    )

                    beta_3way = get_beta_vec_3way4(
                        beta_2way=beta_2way,
                        l1=l1, l2=l2, l3=l3, l4=l4,
                        delta=delta_hat,
                        only_beta=False
                    )

                    beta_4way = get_beta_vec_4way4(
                        beta_3way=beta_3way,
                        l1=l1, l2=l2, l3=l3, l4=l4,
                        tau=tau_hat,
                        only_beta=False
                    )

                    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way])
                    beta_all = beta_all.reshape(-1, 1)
                    eta = X @ beta_all + intercept

                    # Contribution of 3-ways without delta
                    three_ways = X_3way[:, psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                 (beta_hat[i] * beta_hat[j] * beta_hat[l])**2 * \
                                 gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                 gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                 gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)]

                    discard_from_c_3way.append(
                        psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)
                    )

                    assert len(discard_from_c_3way) == 1, "should discard only one three way in update delta"

                    discard_from_c_4way = []
                    # Contribution of 4-ways without delta
                    four_ways = 0
                    for k in range3:
                        four_ways += X_4way[:, phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     (beta_hat[i] * beta_hat[j])**6 * \
                                     (gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] *
                                      gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                                      gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] *
                                      gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                                      gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                                      gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)])**2 * \
                                     delta_hat[psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     tau_hat[phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                     (beta_hat[k] * beta_hat[l])**6

                        discard_from_c_4way.append(
                            phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                        )

                    X_tilde = three_ways + four_ways

                    # Offset term
                    C = X_main @ beta_hat + \
                        X_2way @ beta_2way + \
                        X_3way[:, [idx for idx in range(X_3way.shape[1]) if idx not in discard_from_c_3way]] @ beta_3way[[idx for idx in range(len(beta_3way)) if idx not in discard_from_c_3way]] + \
                        X_4way[:, [idx for idx in range(X_4way.shape[1]) if idx not in discard_from_c_4way]] @ beta_4way[[idx for idx in range(len(beta_4way)) if idx not in discard_from_c_4way]] + \
                        intercept

                    Q_old = Q_bern(
                        X=X, y=y,
                        beta=beta_hat,
                        gamma_vec=gamma_hat,
                        delta_vec=delta_hat,
                        tau_vec=tau_hat,
                        lambda_beta=0,
                        lambda_gamma=0,
                        lambda_delta=lambda_delta,
                        lambda_tau=0,
                        w_beta=1,
                        w_gamma=1,
                        w_delta=1,
                        w_tau=1,
                        l1=l1, l2=l2, l3=l3, l4=l4,
                        intercept=intercept
                    )

                    # Use the minimizer
                    delta_old_value = delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)]

                    delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] = minimizer_Q_bern_delta(
                        X=X_tilde,
                        y=y,
                        C=C,
                        lambda_=lambda_delta,
                        beta_old=delta_old_value,
                        weight=1,
                        scaled=True
                    )  # Intercept is in C

                    # Recompute betas
                    beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta=False)
                    beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta=False)
                    beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta=False)

                    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way])
                    beta_all = beta_all.reshape(-1, 1)
                    eta = X @ beta_all + intercept

                    Q_new = Q_bern(
                        X=X, y=y,
                        beta=beta_hat,
                        gamma_vec=gamma_hat,
                        delta_vec=delta_hat,
                        tau_vec=tau_hat,
                        lambda_beta=0,
                        lambda_gamma=0,
                        lambda_delta=lambda_delta,
                        lambda_tau=0,
                        w_beta=1,
                        w_gamma=1,
                        w_delta=1,
                        w_tau=1,
                        l1=l1, l2=l2, l3=l3, l4=l4,
                        intercept=intercept
                    )

                    if Q_new - Q_old > abs(Q_old) * 1e-2:
                        print("There might be numerical instability in update delta which was taken care of.")

                    # Keep the new estimator only if better
                    if Q_new - Q_old >= 0:
                        delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] = delta_old_value


    # i,k,l
    for i in range1:
        for k in range3:
            for l in range4:
                if (
                    gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] == 0 or
                    gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] == 0 or
                    gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] == 0
                ):
                    delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] = 0
                else:
                    discard_from_c_3way = []
                    beta_2way = get_beta_vec_2way4(
                        beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4,
                        gamma=gamma_hat, only_beta=False
                    )
                    beta_3way = get_beta_vec_3way4(
                        beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4,
                        delta=delta_hat, only_beta=False
                    )
                    beta_4way = get_beta_vec_4way4(
                        beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4,
                        tau=tau_hat, only_beta=False
                    )
                    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way])
                    beta_all = beta_all.reshape(-1, 1)
                    eta = X @ beta_all + intercept

                    # Contribution of 3-ways without delta
                    three_ways = 0
                    three_ways += X_3way[:, psi_table_position_to_vector_index4(
                        [i, k, l], l1=l1, l2=l2, l3=l3, l4=l4
                    )] * (beta_hat[i] * beta_hat[k] * beta_hat[l]) ** 2 * \
                        gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)]
                    discard_from_c_3way.append(
                        psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                    )
                    assert len(discard_from_c_3way) == 1, "should discard only one three way in update delta"

                    discard_from_c_4way = []
                    # Contribution of 4-ways without delta
                    four_ways = 0
                    for j in range2:
                        four_ways += X_4way[:, phi_table_position_to_vector_index4(
                            [i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4
                        )] * (beta_hat[i] * beta_hat[j]) ** 6 * (
                            gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] *
                            gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                            gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] *
                            gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                            gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                            gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)]
                        ) ** 2 * \
                        delta_hat[psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        tau_hat[phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        (beta_hat[k] * beta_hat[l]) ** 6
                        discard_from_c_4way.append(
                            phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                        )

                    X_tilde = three_ways + four_ways
                    C = (
                        X_main @ beta_hat +
                        X_2way @ beta_2way +
                        X_3way[:, [idx for idx in range(X_3way.shape[1]) if idx not in discard_from_c_3way]] @
                        beta_3way[[idx for idx in range(len(beta_3way)) if idx not in discard_from_c_3way]] +
                        X_4way[:, [idx for idx in range(X_4way.shape[1]) if idx not in discard_from_c_4way]] @
                        beta_4way[[idx for idx in range(len(beta_4way)) if idx not in discard_from_c_4way]] +
                        intercept
                    )  # C is the offset term

                    Q_old = Q_bern(
                        X=X, y=y, beta=beta_hat,
                        gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                        lambda_beta=0, lambda_gamma=0, lambda_delta=lambda_delta, lambda_tau=0,
                        w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                        l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept
                    )

                    # Use the minimizer
                    delta_old_value = delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)]
                    delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] = minimizer_Q_bern_delta(
                        X=X_tilde, y=y, C=C,
                        lambda_=lambda_delta, beta_old=delta_old_value,
                        weight=1, scaled=True
                    )  # Intercept is in C

                    beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta=False)
                    beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta=False)
                    beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta=False)
                    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way])
                    beta_all = beta_all.reshape(-1, 1)
                    eta = X @ beta_all + intercept

                    Q_new = Q_bern(
                        X=X, y=y, beta=beta_hat,
                        gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                        lambda_beta=0, lambda_gamma=0, lambda_delta=lambda_delta, lambda_tau=0,
                        w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                        l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept
                    )

                    if Q_new - Q_old > abs(Q_old) * 1e-2:
                        print("There might be numerical instability in update delta which was taken care of.")

                    # Keep the new estimator if only it is better
                    if Q_new - Q_old >= 0:
                        delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] = delta_old_value
                        


    # j,k,l loops
    for j in range2:
        for k in range3:
            for l in range4:
                if (gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)] == 0 or
                    gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1, l2, l3, l4)] == 0 or
                    gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1, l2, l3, l4)] == 0):

                    delta_hat[psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)] = 0
                else:
                    discard_from_c_3way = []
                    beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    gamma=gamma_hat, only_beta=False)  # with delta
                    beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    delta=delta_hat, only_beta=False)  # with gamma+delta
                    beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    tau=tau_hat, only_beta=False)  # with tau
                    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way])
                    beta_all = beta_all.reshape(-1, 1)
                    eta = X @ beta_all + intercept

                    # Contribution of 3-ways without delta
                    three_ways = 0
                    three_ways += (X_3way[:, psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)] *
                                    (beta_hat[j] * beta_hat[k] * beta_hat[l])**2 *
                                    gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1, l2, l3, l4)] *
                                    gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1, l2, l3, l4)] *
                                    gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)])
                    discard_from_c_3way.append(
                        psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)
                    )

                    assert len(discard_from_c_3way) == 1, "should discard only one three way in update delta"

                    # Contribution of 4-ways without delta
                    discard_from_c_4way = []
                    four_ways = 0
                    for i in range1:
                        four_ways += (X_4way[:, phi_table_position_to_vector_index4([i, j, k, l], l1, l2, l3, l4)] *
                                      (beta_hat[i] * beta_hat[j])**6 *
                                      (gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1, l2, l3, l4)] *
                                       gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1, l2, l3, l4)])**2 *
                                      delta_hat[psi_table_position_to_vector_index4([i, j, k], l1, l2, l3, l4)] *
                                      delta_hat[psi_table_position_to_vector_index4([i, j, l], l1, l2, l3, l4)] *
                                      delta_hat[psi_table_position_to_vector_index4([i, k, l], l1, l2, l3, l4)] *
                                      tau_hat[phi_table_position_to_vector_index4([i, j, k, l], l1, l2, l3, l4)] *
                                      (beta_hat[k] * beta_hat[l])**6)
                        discard_from_c_4way.append(
                            phi_table_position_to_vector_index4([i, j, k, l], l1, l2, l3, l4)
                        )

                    X_tilde = three_ways + four_ways

                    C = (X_main @ beta_hat +
                        X_2way @ beta_2way +
                        np.delete(X_3way, discard_from_c_3way, axis=1) @ np.delete(beta_3way, discard_from_c_3way) +
                        np.delete(X_4way, discard_from_c_4way, axis=1) @ np.delete(beta_4way, discard_from_c_4way) +
                        intercept)  # C is offset term

                    Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
                                    tau_vec=tau_hat, lambda_beta=0, lambda_gamma=0,
                                    lambda_delta=lambda_delta, lambda_tau=0,
                                    w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                                    l1=l1, l2=l2, l3=l3, l4=l4,
                                    intercept=intercept)

                    delta_old_value = delta_hat[psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)]

                    # Use minimizer
                    delta_hat[psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)] = \
                        minimizer_Q_bern_delta(X=X_tilde, y=y, C=C, lambda_=lambda_delta,
                                                beta_old=delta_old_value, weight=1, scaled=True)

                    # Recompute betas
                    beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    gamma=gamma_hat, only_beta=False)
                    beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    delta=delta_hat, only_beta=False)
                    beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4,
                                                    tau=tau_hat, only_beta=False)
                    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
                    eta = X @ beta_all + intercept

                    Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat,
                                    tau_vec=tau_hat, lambda_beta=0, lambda_gamma=0,
                                    lambda_delta=lambda_delta, lambda_tau=0,
                                    w_beta=1, w_gamma=1, w_delta=1, w_tau=1,
                                    l1=l1, l2=l2, l3=l3, l4=l4,
                                    intercept=intercept)

                    if Q_new - Q_old > abs(Q_old) * 1e-2:
                        print("There might be numerical instability in update delta which was taken care of.")

                    # Keep the new estimator only if it is better
                    if Q_new - Q_old >= 0:
                        delta_hat[psi_table_position_to_vector_index4([j, k, l], l1, l2, l3, l4)] = delta_old_value

    print("Updated delta")
    return delta_hat




############# Update the vector of main effects  ##################### 


def update_beta(X, y, beta_hat, gamma_hat, delta_hat, tau_hat, lambda_beta,
                l1=21, l2=14, l3=2, l4=3, w=1, intercept=0):
    
    range1 = np.arange(l1)
    range2 = np.arange(l1, l1+l2)
    range3 = np.arange(l1+l2, l1+l2+l3)
    range4 = np.arange(l1+l2+l3, l1+l2+l3+l4)
    
    X_main, X_2way, X_3way, X_4way = (
        X[:, get_ranges4(l1=l1, l2=l2, l3=l3, l4=l4)[i]] for i in range(4)
    )
    
    w = np.ones(len(beta_hat)) if np.isscalar(w) or len(np.atleast_1d(w)) <= 1 else w
    #print("w: ", w)
    beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta=False)
    beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta=False)
    beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta=False)
    beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
    eta = X @ beta_all + intercept
    #print("i")
    for i in range1:
        discard_from_c_main = [i]
        mains = X_main[:, i]
        #print('Mains: ', mains)
        discard_from_c_2way, two_ways = [], 0
        for j in np.concatenate([range2, range3, range4]):
            pos = get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)
            discard_from_c_2way.append(pos)
            two_ways += X_2way[:, pos] * beta_hat[j] * gamma_hat[pos]
        #print("2ways: ", two_ways)
        discard_from_c_3way, three_ways = [], 0
        for j in range2:
            for k in np.concatenate([range3, range4]):
                idx = psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)
                three_ways += X_3way[:, idx] * (beta_hat[j] * beta_hat[k])**2 * \
                              gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                              gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                              gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)] * delta_hat[idx]
                discard_from_c_3way.append(idx)
        for k in range3:
            for l in range4:
                idx = psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                three_ways += X_3way[:, idx] * (beta_hat[k] * beta_hat[l])**2 * \
                              gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                              gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                              gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] * delta_hat[idx]
                discard_from_c_3way.append(idx)
        
        discard_from_c_4way, four_ways = [], 0
        for j in range2:
            for k in range3:
                for l in range4:
                    idx = phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                    four_ways += X_4way[:, idx] * (beta_hat[j])**6 * \
                                 (gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                  gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                  gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                  gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                  gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                  gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)])**2 * \
                                 delta_hat[psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                 delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                 delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                 delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                                 tau_hat[idx] * (beta_hat[k] * beta_hat[l])**6
                    discard_from_c_4way.append(idx)
        
        X_tilde, Z_tilde, T_tilde = mains + two_ways, three_ways, four_ways
        C = intercept + X_main[:, np.arange(len(beta_hat)) != i] @ beta_hat[np.arange(len(beta_hat)) != i] + \
            X_2way[:, [x for x in range(X_2way.shape[1]) if x not in discard_from_c_2way]] @ beta_2way[[x for x in range(len(beta_2way)) if x not in discard_from_c_2way]] + \
            X_3way[:, [x for x in range(X_3way.shape[1]) if x not in discard_from_c_3way]] @ beta_3way[[x for x in range(len(beta_3way)) if x not in discard_from_c_3way]] + \
            X_4way[:, [x for x in range(X_4way.shape[1]) if x not in discard_from_c_4way]] @ beta_4way[[x for x in range(len(beta_4way)) if x not in discard_from_c_4way]]
        #print(np.sum(abs(C-intercept- X_main[:, np.arange(len(beta_hat)) != i] @ beta_hat[np.arange(len(beta_hat)) != i])), "should be 0")
        
        Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, lambda_tau=0,
                        w_beta=w, w_gamma=1, w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)
        
        beta_old_value = beta_hat[i]
        #print(X_tilde)
        #print(Z_tilde)
        #print(T_tilde)
        beta_hat[i] = minimizer_Q_bern_beta(X=X_tilde, Z=Z_tilde, t=T_tilde, y=y, C=C,
                                            lambda_=lambda_beta, beta_old=beta_old_value, weight=w[i], scaled=True)
        
        beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta=False)
        beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta=False)
        beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
        eta = X @ beta_all + intercept
        
        Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, lambda_tau=0,
                        w_beta=w, w_gamma=1, w_delta=1, w_tau=1, l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)
        #print('aaa')
        #print(Q_new-Q_old)
        if Q_new - Q_old > abs(Q_old) * 1e-2:
            print(Q_new, Q_old)
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
        for ik in list(range1) + list(range3) + list(range4):
            small_idx = min(ik, j)
            big_idx = max(ik, j)

            pos = get_position_vec_from_theta_matrix4(
                [small_idx, big_idx], l1=l1, l2=l2, l3=l3, l4=l4
            )
            #print(ik, pos)
            discard_from_c_2way.append(pos)
            two_ways = two_ways + X_2way[:, pos] * beta_hat[ik] * gamma_hat[pos]

        discard_from_c_3way = []
        # Contribution of 3-ways without beta
        three_ways = 0
        for i in range1:
            for k in list(range3) + list(range4):
                pos = psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)
                three_ways = three_ways + X_3way[:, pos] * (beta_hat[i] * beta_hat[k])**2 * \
                    gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    delta_hat[pos]
                discard_from_c_3way.append(pos)

        for k in range3:
            for l in range4:
                pos = psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                three_ways = three_ways + X_3way[:, pos] * (beta_hat[l] * beta_hat[k])**2 * \
                    gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    delta_hat[pos]
                discard_from_c_3way.append(pos)

        discard_from_c_4way = []
        # Contribution of 4-ways without beta
        four_ways = 0
        for i in range1:
            for k in range3:
                for l in range4:
                    pos = phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                    four_ways = four_ways + X_4way[:, pos] * (beta_hat[i])**6 * \
                        (gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)])**2 * \
                        delta_hat[psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        tau_hat[pos] * (beta_hat[k] * beta_hat[l])**6
                    discard_from_c_4way.append(pos)

        X_tilde = mains + two_ways
        Z_tilde = three_ways
        T_tilde = four_ways
        C = intercept + X_main[:, np.arange(X_main.shape[1]) != j] @ beta_hat[np.arange(len(beta_hat)) != j] + \
            X_2way[:, [i for i in range(X_2way.shape[1]) if i not in discard_from_c_2way]] @ beta_2way[[i for i in range(len(beta_2way)) if i not in discard_from_c_2way]] + \
            X_3way[:, [i for i in range(X_3way.shape[1]) if i not in discard_from_c_3way]] @ beta_3way[[i for i in range(len(beta_3way)) if i not in discard_from_c_3way]] + \
            X_4way[:, [i for i in range(X_4way.shape[1]) if i not in discard_from_c_4way]] @ beta_4way[[i for i in range(len(beta_4way)) if i not in discard_from_c_4way]]

        Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, lambda_tau=0,
                        w_beta=w, w_gamma=1, w_delta=1, w_tau=1,
                        l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)

        beta_old_value = beta_hat[j]
        # Use the minimizer
        beta_hat[j] = minimizer_Q_bern_beta(
            X=X_tilde, Z=Z_tilde, t=T_tilde, y=y, C=C,
            lambda_=lambda_beta, beta_old=beta_old_value, weight=w[j], scaled=True
        )

        beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta=False)
        beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta=False)
        beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
        eta = X @ beta_all + intercept

        Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, lambda_tau=0,
                        w_beta=w, w_gamma=1, w_delta=1, w_tau=1,
                        l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)

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
        for il in list(range1) + list(range2) + list(range4):
            small_idx = min(il, k)
            big_idx = max(il, k)

            discard_from_c_2way.append(
                get_position_vec_from_theta_matrix4(
                    [small_idx, big_idx], l1=l1, l2=l2, l3=l3, l4=l4
                )
            )

            two_ways += (
                X_2way[:, get_position_vec_from_theta_matrix4([small_idx, big_idx],
                                                            l1=l1, l2=l2, l3=l3, l4=l4)]
                * beta_hat[il]
                * gamma_hat[get_position_vec_from_theta_matrix4([small_idx, big_idx],
                                                                l1=l1, l2=l2, l3=l3, l4=l4)]
            )

        discard_from_c_3way = []
        # Contribution of 3-ways without beta
        three_ways = 0
        for i in list(range1) + list(range2):
            for l in range4:
                three_ways += (
                    X_3way[:, psi_table_position_to_vector_index4([i, k, l],
                                                                l1=l1, l2=l2, l3=l3, l4=l4)]
                    * (beta_hat[i] * beta_hat[l]) ** 2
                    * gamma_hat[get_position_vec_from_theta_matrix4([i, k],
                                                                    l1=l1, l2=l2, l3=l3, l4=l4)]
                    * gamma_hat[get_position_vec_from_theta_matrix4([k, l],
                                                                    l1=l1, l2=l2, l3=l3, l4=l4)]
                    * gamma_hat[get_position_vec_from_theta_matrix4([i, l],
                                                                    l1=l1, l2=l2, l3=l3, l4=l4)]
                    * delta_hat[psi_table_position_to_vector_index4([i, k, l],
                                                                    l1=l1, l2=l2, l3=l3, l4=l4)]
                )
                discard_from_c_3way.append(
                    psi_table_position_to_vector_index4([i, k, l],
                                                        l1=l1, l2=l2, l3=l3, l4=l4)
                )

        for i in range1:
            for j in range2:
                three_ways += (
                    X_3way[:, psi_table_position_to_vector_index4([i, j, k],
                                                                l1=l1, l2=l2, l3=l3, l4=l4)]
                    * (beta_hat[i] * beta_hat[j]) ** 2
                    * gamma_hat[get_position_vec_from_theta_matrix4([i, j],
                                                                    l1=l1, l2=l2, l3=l3, l4=l4)]
                    * gamma_hat[get_position_vec_from_theta_matrix4([i, k],
                                                                    l1=l1, l2=l2, l3=l3, l4=l4)]
                    * gamma_hat[get_position_vec_from_theta_matrix4([j, k],
                                                                    l1=l1, l2=l2, l3=l3, l4=l4)]
                    * delta_hat[psi_table_position_to_vector_index4([i, j, k],
                                                                    l1=l1, l2=l2, l3=l3, l4=l4)]
                )
                discard_from_c_3way.append(
                    psi_table_position_to_vector_index4([i, j, k],
                                                        l1=l1, l2=l2, l3=l3, l4=l4)
                )

        discard_from_c_4way = []
        # Contribution of 4-ways without beta
        four_ways = 0
        for i in range1:
            for j in range2:
                for l in range4:
                    four_ways += (
                        X_4way[:, phi_table_position_to_vector_index4([i, j, k, l],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                        * (beta_hat[i]) ** 6
                        * (gamma_hat[get_position_vec_from_theta_matrix4([i, k],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                           * gamma_hat[get_position_vec_from_theta_matrix4([i, l],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                           * gamma_hat[get_position_vec_from_theta_matrix4([j, k],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                           * gamma_hat[get_position_vec_from_theta_matrix4([j, l],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                           * gamma_hat[get_position_vec_from_theta_matrix4([k, l],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                           * gamma_hat[get_position_vec_from_theta_matrix4([i, j],
                                                                          l1=l1, l2=l2, l3=l3, l4=l4)]) ** 2
                        * delta_hat[psi_table_position_to_vector_index4([i, j, k],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                        * delta_hat[psi_table_position_to_vector_index4([i, j, l],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                        * delta_hat[psi_table_position_to_vector_index4([i, k, l],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                        * delta_hat[psi_table_position_to_vector_index4([j, k, l],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                        * tau_hat[phi_table_position_to_vector_index4([i, j, k, l],
                                                                        l1=l1, l2=l2, l3=l3, l4=l4)]
                        * (beta_hat[j] * beta_hat[l]) ** 6
                    )
                    discard_from_c_4way.append(
                        phi_table_position_to_vector_index4([i, j, k, l],
                                                            l1=l1, l2=l2, l3=l3, l4=l4)
                    )

        X_tilde = mains + two_ways
        Z_tilde = three_ways
        T_tilde = four_ways
        C = (
            intercept
            + X_main[:, np.arange(X_main.shape[1]) != k] @ beta_hat[np.arange(len(beta_hat)) != k]
            + X_2way[:, [idx for idx in range(X_2way.shape[1]) if idx not in discard_from_c_2way]] @ beta_2way[[idx for idx in range(len(beta_2way)) if idx not in discard_from_c_2way]]
            + X_3way[:, [idx for idx in range(X_3way.shape[1]) if idx not in discard_from_c_3way]] @ beta_3way[[idx for idx in range(len(beta_3way)) if idx not in discard_from_c_3way]]
            + X_4way[:, [idx for idx in range(X_4way.shape[1]) if idx not in discard_from_c_4way]] @ beta_4way[[idx for idx in range(len(beta_4way)) if idx not in discard_from_c_4way]]
        )  # C is the offset term

        Q_old = Q_bern(
            X=X, y=y, beta=beta_hat,
            gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
            lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, lambda_tau=0,
            w_beta=w, w_gamma=1, w_delta=1, w_tau=1,
            l1=l1, l2=l2, l3=l3, l4=l4,
            intercept=intercept
        )

        beta_old_value = beta_hat[k]
        # Use the minimizer
        beta_hat[k] = minimizer_Q_bern_beta(
            X=X_tilde, Z=Z_tilde, t=T_tilde, y=y, C=C,
            lambda_=lambda_beta, beta_old=beta_old_value,
            weight=w[k], scaled=True
        )  # Intercept is in C

        beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4,
                                        gamma=gamma_hat, only_beta=False)
        beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4,
                                        delta=delta_hat, only_beta=False)
        beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4,
                                        tau=tau_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
        eta = X @ beta_all + intercept

        Q_new = Q_bern(
            X=X, y=y, beta=beta_hat,
            gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
            lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, lambda_tau=0,
            w_beta=w, w_gamma=1, w_delta=1, w_tau=1,
            l1=l1, l2=l2, l3=l3, l4=l4,
            intercept=intercept
        )

        if Q_new - Q_old > abs(Q_old) * 1e-2:
            print("There might be numerical instability in update beta which was taken care of.")

        # Keep the new estimator only if it is better than the old one. Ensures numerical stability.
        if Q_new - Q_old >= 0:
            beta_hat[k] = beta_old_value



    #l
    #print("l")
    for l in range4:
        discard_from_c_main = [l]
        mains = np.array(X_main[:, l])

        discard_from_c_2way = []
        # Contribution of 2-ways without beta
        two_ways = 0
        for i in list(range1) + list(range2) + list(range3):
            pos = get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)
            discard_from_c_2way.append(pos)
            two_ways += X_2way[:, pos] * beta_hat[i] * gamma_hat[pos]

        discard_from_c_3way = []
        # Contribution of 3-ways without beta
        three_ways = 0
        for i in range1:
            for k in range3:
                pos3 = psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                three_ways += X_3way[:, pos3] * (beta_hat[i] * beta_hat[k])**2 * \
                    gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    delta_hat[pos3]
                discard_from_c_3way.append(pos3)

        for i in range1:
            for j in range2:
                pos3 = psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)
                three_ways += X_3way[:, pos3] * (beta_hat[i] * beta_hat[j])**2 * \
                    gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                    delta_hat[pos3]
                discard_from_c_3way.append(pos3)

        discard_from_c_4way = []
        # Contribution of 4-ways without beta
        four_ways = 0
        for i in range1:
            for j in range2:
                for k in range3:
                    pos4 = phi_table_position_to_vector_index4([i, j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)
                    four_ways += X_4way[:, pos4] * (beta_hat[i])**6 * \
                        (gamma_hat[get_position_vec_from_theta_matrix4([i, k], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([i, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([j, k], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([j, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([k, l], l1=l1, l2=l2, l3=l3, l4=l4)] *
                         gamma_hat[get_position_vec_from_theta_matrix4([i, j], l1=l1, l2=l2, l3=l3, l4=l4)])**2 * \
                        delta_hat[psi_table_position_to_vector_index4([i, j, k], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([i, j, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([i, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        delta_hat[psi_table_position_to_vector_index4([j, k, l], l1=l1, l2=l2, l3=l3, l4=l4)] * \
                        tau_hat[pos4] * (beta_hat[j] * beta_hat[k])**6
                    discard_from_c_4way.append(pos4)

        X_tilde = mains + two_ways
        Z_tilde = three_ways
        T_tilde = four_ways

        C = intercept + X_main[:, [i for i in range(X_main.shape[1]) if i != l]] @ beta_hat[[i for i in range(len(beta_hat)) if i != l]] + \
            X_2way[:, [i for i in range(X_2way.shape[1]) if i not in discard_from_c_2way]] @ beta_2way[[i for i in range(len(beta_2way)) if i not in discard_from_c_2way]] + \
            X_3way[:, [i for i in range(X_3way.shape[1]) if i not in discard_from_c_3way]] @ beta_3way[[i for i in range(len(beta_3way)) if i not in discard_from_c_3way]] + \
            X_4way[:, [i for i in range(X_4way.shape[1]) if i not in discard_from_c_4way]] @ beta_4way[[i for i in range(len(beta_4way)) if i not in discard_from_c_4way]]

        Q_old = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, lambda_tau=0,
                        w_beta=w, w_gamma=1, w_delta=1, w_tau=1,
                        l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)

        beta_old_value = beta_hat[l]
        # Use the minimizer
        beta_hat[l] = minimizer_Q_bern_beta(X=X_tilde, Z=Z_tilde, t=T_tilde, y=y, C=C,
                                            lambda_=lambda_beta, beta_old=beta_old_value,
                                            weight=w[l], scaled=True)

        beta_2way = get_beta_vec_2way4(beta=beta_hat, l1=l1, l2=l2, l3=l3, l4=l4, gamma=gamma_hat, only_beta=False)
        beta_3way = get_beta_vec_3way4(beta_2way=beta_2way, l1=l1, l2=l2, l3=l3, l4=l4, delta=delta_hat, only_beta=False)
        beta_4way = get_beta_vec_4way4(beta_3way=beta_3way, l1=l1, l2=l2, l3=l3, l4=l4, tau=tau_hat, only_beta=False)
        beta_all = np.concatenate([beta_hat, beta_2way, beta_3way, beta_4way]).reshape(-1, 1)
        eta = X @ beta_all + intercept

        Q_new = Q_bern(X=X, y=y, beta=beta_hat, gamma_vec=gamma_hat, delta_vec=delta_hat, tau_vec=tau_hat,
                        lambda_beta=lambda_beta, lambda_gamma=0, lambda_delta=0, lambda_tau=0,
                        w_beta=w, w_gamma=1, w_delta=1, w_tau=1,
                        l1=l1, l2=l2, l3=l3, l4=l4, intercept=intercept)

        if Q_new - Q_old > abs(Q_old) * 1e-2:
            print("There might be numerical instability in update beta which was taken care of.")

        if Q_new - Q_old >= 0:
            beta_hat[l] = beta_old_value

    print("Updated beta")
    return beta_hat







