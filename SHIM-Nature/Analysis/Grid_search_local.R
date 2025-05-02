directory_here<-getwd()
# Specify the directory containing the R scripts
directory_libs <- "../libs"
setwd(directory_libs)
# Get a list of all R script files in the directory
files <- as.list( list.files(directory_libs, pattern = "\\.R$", full.names = TRUE) ) 

getwd()
# Source each file
lapply(files, source)
class(files)
setwd(directory_here)

l1=21
l2=14
l3=2
l4=3

range_main<-unlist(get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[1])
range_theta<-unlist(get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[2])
range_psi<-unlist(get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[3])
range_phi<-unlist(get_ranges4(l1=l1,l2=l2,l3=l3,l4=l4)[4])


# Load data
data<-read_load_BHA_data()
xxxx.all<-data$xxxx.all
y.all<-data$y.all
x.all<-xxxx.all[,range_main]



############    GLM on main effects for initialization   ########################################


res_lasso_all<-irlasso.cb(X=x.all, Y=y.all, lambda=0, w.lambda=NULL, beta0=NULL,
                          centering=FALSE, scaling=FALSE, intercept=T,
                          maxit=20, tol=0.05, sd.tol=1e-6,
                          verbose=TRUE)

# Use lasso for init for SHIM

coefs_lasso<-array(res_lasso_all$beta[-1,1,1])
interc_init<-res_lasso_all$beta[1,1,1]

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


range_lmd_beta<-c(5e-5) # One can chose any list of values of lambda to try. The current values are the ones for the final estimation. 
range_lmd_gamma<-c( 7e-4) # e.g. c(3e-4, 7e-4,1e-3,3e-3)
lambda_delta<-1e-3
lambda_tau<-3e-3
interc_init<-0
for (lambda_beta in range_lmd_beta){
  for (lambda_gamma in range_lmd_gamma)
  { 
    paste("lmdbeta: ", lambda_beta, "lmdgamma: ", lambda_gamma)
    my_shim<-SHIM_4way(X=xxxx.all, y=y.all, beta_init = beta_hat, gamma_init = gamma_hat, delta_init = delta_hat, tau_init=tau_hat, l1=l1, l2=l2, l3=l3, l4=l4)
    
    start.time <- Sys.time()
    
    
    fitted<-my_shim$fit(X=xxxx.all, y=y.all, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, lambda_tau=lambda_tau, w_beta = w_beta, 
                        w_gamma = 1, w_delta = 1, w_tau=1, tol=5e-3, compute_Q = Q_bern, intercept = interc_init, use_intercept = TRUE)
    
    
    beta_all_shim<-fitted$beta_all   
    intercept_all<-fitted$intercept
    coefs_all<-c(intercept_all, beta_all_shim)
    saveRDS(coefs_all, paste("../coefficients/Coefficients_lmd_", lambda_beta,"_", lambda_gamma,'_', lambda_delta) )
    
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,4)
    time.taken
    print(time.taken)
    
  }
}

cfs<-readRDS(paste("../coefficients/Coefficients_lmd_", 5e-05,"_", 7e-04,'_', 1e-3))
max(abs(cfs-coefs_all))
