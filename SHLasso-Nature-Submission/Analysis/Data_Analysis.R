library(ggplot2)

directory_here<-getwd()
# Specify the directory containing the R scripts
directory_libs <- "../libs"
setwd(directory_libs)
# Get a list of all R script files in the directory
files <- as.list( list.files(directory_libs, pattern = "\\.R$", full.names = TRUE) ) 

getwd()
# Source each file from libs
lapply(files, source)
class(files)
setwd(directory_here)

# Specify the number of levels for each factor
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
x.all<-xxxx.all[,range_main] # get the matrix of main effects



##########    GLM on main effects for initialization   ##########  


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



##########  Do CB SHIM with given initialization  ##########  (note that different initialization can be used)


lambda_beta<-5e-5
lambda_gamma<-7e-4
lambda_delta<-1e-3
lambda_tau<-3e-3
interc_init<-0

my_shim<-SHIM_4way(X=xxxx.all, y=y.all, beta_init = beta_hat, 
                   gamma_init = gamma_hat, delta_init = delta_hat, tau_init=tau_hat, l1=l1, l2=l2, l3=l3, l4=l4)
    
start.time <- Sys.time()
    
    
fitted<-my_shim$fit(X=xxxx.all, y=y.all, lambda_beta = lambda_beta, lambda_gamma = lambda_gamma, lambda_delta = lambda_delta, lambda_tau=lambda_tau, 
                    w_beta = w_beta, tol=5e-3, compute_Q = Q_bern, intercept = interc_init, use_intercept = TRUE)
    
    
beta_all_shim<-fitted$beta_all   
intercept_shim<-fitted$intercept
coefs_all<-c(intercept_shim, beta_all_shim)
fn_coefs<-paste("../coefficients/Final_coefficients.rds") # file where the coefficients are stored
saveRDS(coefs_all, fn_coefs )
end.time <- Sys.time()
time.taken <- round(end.time - start.time,4)
time.taken
print(time.taken)
    
# Print the number of coefficients
p<-sum(coefs_all!=0)
cat("Number of coefficients kept in the model is", paste0(p, "."), "\n")

y.pred<-kappa1(xxxx.all%*%beta_all_shim+intercept_shim) # prediction 

# Prepare to plot
ObservedYield <- as.vector(y.all)
PredictedYield <- as.vector(y.pred)
cat("The R2 score on entire data using only", p, "coefficients is", paste0( round(r2(y.all, y.pred),3),".") ,"\n")

# Create the table with final coefficients
Create_table_results(fn_coefs = fn_coefs)


