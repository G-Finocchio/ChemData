# Functions can be accessed via the command source("tools.r")



# The next function implements Partial Least Squares on
# uni/multivariate response Y and design X using KERNELPLS
# The function allows for centering/scaling (FALSE by default)
# The function allows for intercept (FALSE by default)

kernelpls2 <- function(X, Y, ncomp, 
                       centering=FALSE, scaling=FALSE, intercept=FALSE,
                       maxit=50, tol=1e-8,
                       verbose=FALSE){
  
  if (verbose) print("Performing KERNELPLS...")
  
  # Get variables
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Center variables (column-wise)
  mu.X <- rep(0, ncol(X))
  mu.Y <- rep(0, ncol(Y))
  if (centering) {
    mu.X <- apply(X, 2, mean)
    mu.Y <- apply(Y, 2, mean)
    X <- apply(X, 2, function(v){v-mean(v)})
    Y <- apply(Y, 2, function(v){v-mean(v)})
  }
  
  # Scale variables (column-wise)
  sd.X <- rep(1, ncol(X))
  sd.Y <- rep(1, ncol(Y))
  if (scaling) {
    sd.X <- apply(X, 2, function(v){max(sd(v), tol)})
    sd.Y <- apply(Y, 2, function(v){max(sd(v), tol)})
    X <- apply(X, 2, function(v){v/max(sd(v), tol)})
    Y <- apply(Y, 2, function(v){v/max(sd(v), tol)})
  }
  
  # Intercept
  if (intercept) { 
    X <- cbind(1, X) # intercept is added as first column
  }    
  
  # Get parameters
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(Y)[2]
  
  ## Initialization
  R <- matrix(0, ncol=ncomp, nrow=p)        # projection
  P <- matrix(0, ncol=ncomp, nrow=p)        # X loadings
  tQ <- matrix(0, ncol=m, nrow=ncomp)       # Y loadings; transposed
  B <- array(0, c(p, m, ncomp))             # coefficients
  W <- P                                    # weights
  U <- matrix(0, ncol=ncomp, nrow=n)        # scores
  T <- matrix(0, ncol=ncomp, nrow=n)        # scores
  tsqs <- rep.int(1, ncomp)                 # t't
  
  ## 1.
  XtY <- crossprod(X, Y)
  
  for (a in 1:ncomp) {
    
    ## 2.
    if (m == 1) {
      
      w.a <- XtY / sqrt(c(crossprod(XtY)))
      
    } else {
      
      if (m < p) {
        
        q <- eigen(crossprod(XtY), symmetric = TRUE)$vectors[,1]
        w.a <- XtY %*% q
        w.a <- w.a / sqrt(c(crossprod(w.a)))
        
      } else {
        
        w.a <- eigen(XtY %*% t(XtY), symmetric = TRUE)$vectors[,1]
        
      }
    }
    
    ## 3.
    r.a <- w.a
    
    if (a > 5) {
      
      ## This is faster when a > 5:
      r.a <- r.a - colSums(crossprod(w.a, P[,1:(a-1), drop=FALSE])%*%t(R[,1:(a-1), drop=FALSE]))
      
    } else if (a > 1) {
      
      for (j in 1:(a - 1)) {
        r.a <- r.a - c(P[,j] %*% w.a) * R[,j]
      }
      
    }
    
    ## 4.
    t.a <- X %*% r.a
    tsq <- c(crossprod(t.a))
    p.a <- crossprod(X, t.a) / tsq
    q.a <- crossprod(XtY, r.a) / tsq
    
    ## 5.
    XtY <- XtY - (tsq * p.a) %*% t(q.a)
    
    ## 6. 7. 8.
    R[,a]  <- r.a
    P[,a]  <- p.a
    tQ[a,] <- q.a
    B[,,a] <- R[,1:a, drop=FALSE] %*% tQ[1:a,, drop=FALSE]
    tsqs[a] <- tsq
    
    ## Extra step to calculate Y scores:
    u.a <- Y %*% q.a / c(crossprod(q.a))
    
    ## Make u orth to previous X scores:
    if (a > 1) {
      u.a <- u.a - T %*% (crossprod(T, u.a) / tsqs)
    }
    U[,a] <- u.a
    T[,a] <- t.a
    W[,a] <- w.a
    
  }
  
  # Compute coeffs for original variables
  B.tilde <- array(0, c(p, m, ncomp))
  
  if (intercept) {
    
    for (nc in 1:ncomp) {
      for (mc in 1:m) {
        B.tilde[1,mc,nc] <- sd.Y[mc]*B[1,mc,nc] + mu.Y[mc] - sd.Y[mc]*(mu.X/sd.X)%*%B[2:p,mc,nc]  # intercept
        B.tilde[2:p,mc,nc] <- sd.Y[mc]*B[2:p,mc,nc]/sd.X                                          # core
      }
    }
    
  } else {
    
    for (nc in 1:ncomp) {
      for (mc in 1:m) {
        B.tilde[1:p,mc,nc] <- sd.Y[mc]*B[1:p,mc,nc]/sd.X  # no intercept, core
      }
    }
    
  }
  
  return(list(B=B,      # coeffs for centered/scaled variables
              T=T,      # scores for centered/scaled variables
              P=P,      # X loadings for centered/scaled variables
              W=W,      # weights for centered/scaled variables
              U=U,      # scores for centered/scaled variables
              Q=t(tQ),  # Y loadings for centered/scaled variables
              R=R,      # projection for centered/scaled variables
              muX=mu.X, muY=mu.Y,  # vectors of column-wise mean
              sdX=sd.X, sdY=sd.Y,  # vectors of column-wise sd
              Btilde=B.tilde       # coeffs for unscaled variables
              )
         )
  
}



# The next function implements Partial Least Squares in any case
# For the parameter method, "kernel" and "nipals" call the respective algorithms
# The function allows for centering/scaling (kernel only)
# The function allows for intercept (kernel only)

plslm <- function(X, Y, ncomp,
                  centering=FALSE, scaling=FALSE, intercept=FALSE,
                  maxit=50, tol=1e-8, 
                  method="kernel",
                  verbose=FALSE){
  
  if (method=="kernel") {
    if (verbose) print("Performing PLS with kernelpls...")
    return(kernelpls2(X=X, Y=Y, ncomp=ncomp, 
                      centering=centering, scaling=scaling, intercept=intercept, 
                      maxit=maxit, tol=tol,
                      verbose=verbose))
  }
  else {
    if (verbose) print("Invalid method selection, performing kernelpls as default...")
    return(kernelpls2(X=X, Y=Y, ncomp=ncomp, 
                      centering=centering, scaling=scaling, intercept=intercept, 
                      maxit=maxit, tol=tol,
                      verbose=verbose))
  }
  
}



# The next function implements PLSGLM with Continuous Bernoulli response
# The parameter beta0 is the initialization, choosing NULL initializes at MLE 
# The function allows for centering/scaling (TRUE by default)
# The function allows for intercept (TRUE by default)

plsglm.cb <- function(X, Y, ncomp, beta0=NULL,
                      centering=TRUE, scaling=TRUE, intercept=TRUE,
                      maxit=50, tol=1e-6,
                      verbose=FALSE){
  
  if (verbose) print("Performing PLSGLM-NEW")
  
  # CB Likelihood
  kappa1 <- function(x) {1-1/x-1/(1-exp(x))}
  kappa2 <- function(x) {1/x^2+1/(2-2*cosh(x))}
  
  # Get variables
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Center variables (column-wise)
  mu.X <- rep(0, ncol(X))
  mu.Y <- rep(0, ncol(Y))
  if (centering) {
    mu.X <- apply(X, 2, mean)
    #mu.Y <- apply(Y, 2, mean)
    X <- apply(X, 2, function(v){v-mean(v)})
    #Y <- apply(Y, 2, function(v){v-mean(v)})
  }
  
  # Scale variables (column-wise)
  sd.X <- rep(1, ncol(X))
  sd.Y <- rep(1, ncol(Y))
  if (scaling) {
    sd.X <- apply(X, 2, function(v){max(sd(v), tol)})
    #sd.Y <- apply(Y, 2, function(v){max(sd(v), tol)})
    X <- apply(X, 2, function(v){v/max(sd(v), tol)})
    #Y <- apply(Y, 2, function(v){v/max(sd(v), tol)})
  }
  
  # Intercept
  if (intercept) { 
    X <- cbind(1, X)
  }
  
  # Get parameters
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(Y)[2]
  
  # Initialization
  if (is.null(beta0)){ beta0 <- rep(0, p) }
  
  beta.old <- array(beta0, dim=c(p, 1, ncomp))
  Z.old <- array(0, dim=c(n, 1, ncomp))
  W.old <- array(0, dim=c(n, n, ncomp))
  
  beta <- array(beta0, dim=c(p, 1, ncomp))
  Z <- array(0, dim=c(n, 1, ncomp))
  W <- array(0, dim=c(n, n, ncomp))
  
  R <- list()
  
  nc <- ncomp
  it.stop <- rep(0, ncomp)
  
  # Run until convergence or stop
  counter <- 0
  
  repeat{
    if (verbose) print(paste("Performing PLSGLM iter", counter, sep=" "))
    
    # Keep old variables
    beta.old <- beta
    Z.old <- Z
    W.old <- W
    
    # Compute W and Z
    eta <- array(0, dim=c(n, 1, ncomp))
    k.p <- array(0, dim=c(n, 1, ncomp))
    k.pp <- array(0, dim=c(n, 1, ncomp))
    
    if (length(nc)==0) {
      if (verbose) print("No comps left...")
      break
    }
    
    if (verbose) print(paste("Comps left", paste(nc, collapse=" ")))
    for (m in nc) {
      
      eta[,,m] <- X%*%beta[,,m]
      
      k.p[,,m] <- kappa1(eta[,,m])
      k.p[,,m][which(k.p[,,m]<0.001)] <- 0.001
      k.p[,,m][which(k.p[,,m]>0.999)] <- 0.999
      
      k.pp[,,m] <- kappa2(eta[,,m])
      k.pp[,,m][which(k.pp[,,m]<0.001)] <- 0.001
      k.pp[,,m][which(k.pp[,,m]>0.999)] <- 0.999
      
      W[,,m] <- diag(as.vector(k.pp[,,m]))
      
      Z[,,m] <- eta[,,m] + diag(1/as.vector(k.pp[,,m]))%*%(Y-as.vector(k.p[,,m]))
      
      # Weighted X.W and Z.W on largest component
      X.W <- as.matrix(sqrt(W[,,m])%*%X)
      Z.W <- as.matrix(sqrt(W[,,m])%*%Z[,,m])
      
      # Compute coefficients for Z.W ~ X.W
      # No center/scale X.W, Z.W
      # No intercept X.W
      fit.pls <- plslm(X=X.W, Y=Z.W, ncomp=m,
                       centering=FALSE, scaling=FALSE, intercept=FALSE)
      R[[m]] <- fit.pls$R
      beta[,,m] <- fit.pls$Btilde[,,m]
      
    }
    
    #print(paste("Norm of beta", paste(apply(beta^2, 3, sum), collapse=" "), sep=" "))
    
    epsilon <- sqrt(apply((beta-beta.old)^2, 3, sum)/apply((beta.old)^2, 3, sum) )
    #print(paste("Divergence", paste(epsilon, collapse=" "), sep=" "))
    print(paste("Min Divergence", min(epsilon[nc]), sep=" "))
    
    log.like <- apply(beta, 3, function(v) sum(kappa1(X%*%v)*Y+(1-kappa1(X%*%v))*(1-Y)))
    #print(paste("Loglike", paste(log.like, collapse=" "), sep=" "))
    
    log.like.ratio <- log.like - apply(beta.old, 3, function(v) sum(kappa1(X%*%v)*Y+(1-kappa1(X%*%v))*(1-Y)))
    #print(paste("Loglike ratio", paste(log.like.ratio, collapse=" "), sep=" "))
    
    if (sum(is.nan(epsilon[nc]))>0) {
      nan.stop <- which(is.nan(epsilon))
      if (verbose) print(paste("Divergence NaN comps", paste(nan.stop, collapse=" ")))
      for (m in nc) {
        beta[,,m] <- beta.old[,,m]
        Z[,,m] <- Z.old[,,m]
        W[,,m] <- W.old[,,m]
      }
      nc <- setdiff(nc, nan.stop)
    }
    
    if (min(epsilon[nc])<tol) { 
      nc.stop <- which(epsilon<tol)
      it.stop[nc.stop] <- counter
      if (verbose) print(paste("Divergence stop comps", paste(nc.stop, collapse=" ")))
      nc <- setdiff(nc, nc.stop)
    }
    
    if (counter==maxit) { 
      if (verbose) print("Maximum iterarion, no convergence...")
      it.stop[which(it.stop==0)] <- counter
      break
    }
    else {
      counter <- counter+1
    }
    
  }
  
  # Compute coeffs for original variables
  beta.tilde <- beta*0
  
  for (m in 1:ncomp) {
    
    if (intercept) {
      
      beta.tilde[,,m][1] <- sd.Y[1]*beta[,,m][1] + mu.Y[1] - sd.Y[1]*(mu.X/sd.X)%*%beta[,,m][2:p] # intercept
      beta.tilde[,,m][2:p] <- sd.Y[1]*beta[,,m][2:p]/sd.X                                         # core
      
    } else {
      
      for (mc in 1:m) {
        beta.tilde[,,m][1:p] <- sd.Y[1]*beta[,,m][1:p]/sd.X                                       # no intercept, core
      }
    }
  }
  
  gc()
  
  return(list(BETA=beta.tilde, 
              beta=beta, 
              Z=Z, W=W, 
              R=R,
              it=it.stop))
}

