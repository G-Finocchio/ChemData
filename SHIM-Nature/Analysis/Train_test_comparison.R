library(ggplot2)
library(randomForest)
directory_here<-getwd()
directory_libs <- "../libs"
setwd(directory_libs)
# Get a list of all R script files in the directory
files <- as.list( list.files(directory_libs, pattern = "\\.R$", full.names = TRUE) ) 
# Source each file
lapply(files, source)
setwd(directory_here)

l1=21
l2=14
l3=2
l4=3

# Ranges for main effects, two-way interactions, three-way interactions, and four-way interactions
range_main<-unlist(get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[1])
range_theta<-unlist(get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[2])
range_psi<-unlist(get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[3])
range_phi<- unlist(get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[4])

# Load data
data<-read_load_BHA_data()
x.train_main<-data$x.train[, range_main]
x.train<-data$x.train
y.train<-data$y.train
x.test_main<-data$x.test[, range_main]
x.test<-data$x.test
y.test<-data$y.test


# RF with best parameters
best_mtry<-24
fit.rf=randomForest(y.train~.,data=data.frame(x.train_main), mtry = best_mtry, ntree = 500)
pred.rf=predict(fit.rf,data.frame(x.test))
r2(y.test, pred.rf)

ObservedYield <- as.vector(y.test)
PredictedYield <- as.vector(pred.rf)

# Plot prediction vs observed
plot_predicted_yield(PredictedYield=PredictedYield, ObservedYield=ObservedYield, title="RF", fn_save="../Results/Rf_pred_test.pdf")


# cB Lasso
best_lmd=3e-4
res_lasso<-irlasso.cb(X=x.train, Y=y.train, lambda=3e-4, w.lambda=NULL, beta0=NULL,
                          centering=FALSE, scaling=FALSE, intercept=T,
                          maxit=20, tol=0.05, sd.tol=1e-6,
                          verbose=TRUE)

coefs_lasso<-array(res_lasso$beta[-1,1,1])
sum(coefs_lasso!=0)
intercept<-res_lasso$beta[1,1,1]
pred.lasso<-kappa1(x.test%*%coefs_lasso+intercept) # prediction 
r2(y.test,pred.lasso)

ObservedYield <- as.vector(y.test)
PredictedYield <- as.vector(pred.lasso)
print(pred.lasso[2:10])
print(pred.rf[2:10])

# Plot prediction vs observed
plot_predicted_yield(PredictedYield=PredictedYield, ObservedYield=ObservedYield, title="CB Lasso", fn_save="../Results/Lasso_pred_test.pdf")


################ Do SHIM with best coefs on train-test #############################

# Remove the comment symbol to run SHIM from scratch on train set!

# GLM on main effects for initialization
res_lasso_train<-irlasso.cb(X=x.train_main, Y=y.train, lambda=0, w.lambda=NULL, beta0=NULL,
                          centering=FALSE, scaling=FALSE, intercept=T,
                          maxit=20, tol=0.05, sd.tol=1e-6,
                          verbose=TRUE)

# Use lasso for init for SHIM

coefs_lasso<-array(res_lasso_train$beta[-1,1,1])
interc_init<-res_lasso_train$beta[1,1,1]

beta_main_lasso<-coefs_lasso[unlist(range_main)]
beta_2way_lasso<-rep(0,length(unlist(range_theta)))
beta_3way_lasso<-rep(0, length(unlist(range_psi)))
beta_4way_lasso<-rep(0, length(unlist(range_phi)))
w_beta=pmin(10^7, abs(1/beta_main_lasso))
beta_hat<-beta_main_lasso
gamma_hat<- beta_2way_lasso
delta_hat<-beta_3way_lasso
tau_hat<-beta_4way_lasso



#####################################   DO GRID SEARCH   ###################


range_lmd_beta<-c(1e-5)
range_lmd_gamma<-c(1e-4)
lambda_delta<-1e-3
lambda_tau<-3e-3
interc_init<-0 # Search again for intercept as first step of SHIM
for (lambda_beta in range_lmd_beta){
  for (lambda_gamma in range_lmd_gamma)
  { 
    paste("lmdbeta: ", lambda_beta, "lmdgamma: ", lambda_gamma)
    my_shim<-SHIM_4way(X=x.train, y=y.train, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, tau_init=tau_hat, l1=l1, l2=l2, l3=l3, l4=l4)
    
    start.time <- Sys.time()
    
    fitted<-my_shim$fit(X=x.train, y=y.train, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, lambda_tau=lambda_tau, w_beta = w_beta,  tol=2e-2, compute_Q = Q_bern, intercept = interc_init, use_intercept = TRUE)
    
    beta_all_shim<-fitted$beta_all   
    intercept_all<-fitted$intercept
    coefs_all<-c(intercept_all, beta_all_shim)
    saveRDS(coefs_all, paste0("../coefficients/Coefficients_traintest_lmd_", lambda_beta,"_", lambda_gamma,'_', lambda_delta, ".rds") )
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,4)
    print(time.taken)
    
  }
}

coefs_all<-readRDS(paste0("../coefficients/Coefficients_traintest_lmd_", 1e-5,"_", 1e-4,'_', 1e-3, ".rds"))
intercept_all<-coefs_all[1]
beta_all<-coefs_all[-1]

pred.shim=kappa1(x.test%*%beta_all + intercept_all)
r2(y.test, pred.shim)
ObservedYield <- as.vector(y.test)
PredictedYield <- as.vector(pred.shim)

# Plot prediction vs observed
plot_predicted_yield(PredictedYield=PredictedYield, ObservedYield=ObservedYield, title=" CB SH Lasso", fn_save="../Results/SHLasso_pred_test.pdf")




