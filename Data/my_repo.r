# Functions can be accessed via the command source("my_repo.r")


# The next function computes the performance of a classification probability

my.prob.score <- function(Y, prob, time.train, time.test, cut.level=0.5, tp.level=0){
  
  # split times
  fit.train <- prob[time.train]
  fit.test <- prob[time.test]
  
  # threshold estimator
  fit.train[which(fit.train<=cut.level)] <- 0
  fit.train[which(fit.train>cut.level)] <- 1
  
  fit.test[which(fit.test<=cut.level)] <- 0
  fit.test[which(fit.test>cut.level)] <- 1
  
  # compute TPR/FNR train
  pos.train <- which(Y[time.train]==1)
  true.pos.train <- NULL
  for (event.train in which(fit.train==1) ) {
    event.score <- min(abs(event.train - pos.train))
    if (event.score <= tp.level){
      true.pos.train <- c(true.pos.train, event.train)
    }
  }
  true.pos.train <- sort(unique(true.pos.train))
  true.pos.rate.train <- round(length(true.pos.train)/length(pos.train), digits=6)
  if (is.na(true.pos.rate.train)) true.pos.rate.train <- 1
  false.neg.rate.train <- round(1 - true.pos.rate.train, digits=6)
  
  # compute TNR/FPR train
  neg.train <- which(Y[time.train]==0)
  true.neg.train <- which(fit.train==0 & fit.train==Y[time.train])
  true.neg.rate.train <- round(length(true.neg.train)/length(neg.train), digits=6)
  if (is.na(true.neg.rate.train)) true.neg.rate.train <- 1
  false.pos.rate.train <- round(1 - true.neg.rate.train, digits=6)
  
  # compute ACC train
  acc.train <- round((length(true.pos.train)+length(true.neg.train))/length(time.train), digits=6)
  
  # compute TPR/FNR test
  pos.test <- which(Y[time.test]==1)
  true.pos.test <- NULL
  for (event.test in which(fit.test==1) ) {
    event.score <- min(abs(event.test - pos.test))
    if (event.score <= tp.level){
      true.pos.test <- c(true.pos.test, event.test)
    }
  }
  true.pos.test <- sort(unique(true.pos.test))
  true.pos.rate.test <- round(length(true.pos.test)/length(pos.test), digits=6)
  if (is.na(true.pos.rate.test)) true.pos.rate.test <- 1
  false.neg.rate.test <- round(1 - true.pos.rate.test, digits=6)
  
  # compute TNR/FPR test
  neg.test <- which(Y[time.test]==0)
  true.neg.test <- which(fit.test==0 & fit.test==Y[time.test])
  true.neg.rate.test <- round(length(true.neg.test)/length(neg.test), digits=6)
  if (is.na(true.neg.rate.test)) true.neg.rate.test <- 1
  false.pos.rate.test <- round(1 - true.neg.rate.test, digits=6)
  
  # compute ACC test
  acc.test <- round((length(true.pos.test)+length(true.neg.test))/length(time.test), digits=6)
  
  score <- list(TPRtrain=true.pos.rate.train,
                FNRtrain=false.neg.rate.train,
                TNRtrain=true.neg.rate.train,
                FPRtrain=false.pos.rate.train,
                ACCtrain=acc.train,
                TPRtest=true.pos.rate.test,
                FNRtest=false.neg.rate.test,
                TNRtest=true.neg.rate.test,
                FPRtest=false.pos.rate.test,
                ACCtest=acc.test)
  
  return(score)
}



# The next function visualizes computes the performance of score list

my.score.visual <- function(list.score, k.max, filename.scores, height=12, width=24) {
  
  TPRtrain <- NULL
  TPRtest <- NULL
  TNRtrain <- NULL
  TNRtest <- NULL
  ACCtrain <- NULL
  ACCtest <- NULL
  l2.score <- NULL
  m.score <- NULL
  
  for (seed_id in 1:k.max) {
    if (is.numeric(list.score[[seed_id]]$TPRtrain) & 
        is.numeric(list.score[[seed_id]]$TPRtest) &
        is.numeric(list.score[[seed_id]]$TNRtrain) &
        is.numeric(list.score[[seed_id]]$TNRtest) &
        is.numeric(list.score[[seed_id]]$ACCtrain) &
        is.numeric(list.score[[seed_id]]$ACCtest) &
        is.numeric(list.score[[seed_id]]$L2) &
        is.numeric(list.score[[seed_id]]$ncomp)) {
      
      TPRtrain <- c(TPRtrain, list.score[[seed_id]]$TPRtrain)
      TPRtest <- c(TPRtest, list.score[[seed_id]]$TPRtest)
      TNRtrain <- c(TNRtrain, list.score[[seed_id]]$TNRtrain)
      TNRtest <- c(TNRtest, list.score[[seed_id]]$TNRtest)
      ACCtrain <- c(ACCtrain, list.score[[seed_id]]$ACCtrain)
      ACCtest <- c(ACCtest, list.score[[seed_id]]$ACCtest)
      l2.score <- c(l2.score, list.score[[seed_id]]$L2)
      m.score <- c(m.score, list.score[[seed_id]]$ncomp)
      
    }
  }
  
  pdf(file=paste(filename.scores, ".pdf", sep=""), height=height, width=width)
  par(mfrow=c(2,4))
  hist(TPRtrain, k.max, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  abline(v=mean(TPRtrain), col=2)
  abline(v=c(mean(TPRtrain)-sd(TPRtrain), mean(TPRtrain)+sd(TPRtrain)), col=4)
  legend("topleft", col=c(1,2,4), lwd=2, lty=1, cex=1.5, c(paste("TPR score"), 
                                                           paste("mean", format(mean(TPRtrain), nsmall=6)), 
                                                           paste("sdev", format(sd(TPRtrain), nsmall=6)) ) ) 
  
  hist(TNRtrain, k.max, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  abline(v=mean(TNRtrain), col=2)
  abline(v=c(mean(TNRtrain)-sd(TNRtrain), mean(TNRtrain)+sd(TNRtrain)), col=4)
  legend("topleft", col=c(1,2,4), lwd=2, lty=1, cex=1.5, c(paste("TNR score"), 
                                                           paste("mean", format(mean(TNRtrain), nsmall=6)), 
                                                           paste("sdev", format(sd(TNRtrain), nsmall=6)) ) ) 
  
  hist(ACCtrain, k.max, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  abline(v=mean(ACCtrain), col=2)
  abline(v=c(mean(ACCtrain)-sd(ACCtrain), mean(ACCtrain)+sd(ACCtrain)), col=4)
  legend("topleft", col=c(1,2,4), lwd=2, lty=1, cex=1.5, c(paste("ACC score"), 
                                                           paste("mean", format(mean(ACCtrain), nsmall=6)), 
                                                           paste("sdev", format(sd(ACCtrain), nsmall=6)) ) ) 
  
  hist(m.score, k.max, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  abline(v=mean(m.score), col=2)
  abline(v=c(mean(m.score)-sd(m.score), mean(m.score)+sd(m.score)), col=4)
  legend("topright", col=c(1,2,4), lwd=2, lty=1, cex=1.5, c(paste("ncomp"),
                                                            paste("mean", format(mean(m.score), nsmall=6)),
                                                            paste("sdev", format(sd(m.score), nsmall=6)) ) )
  
  hist(TPRtest, k.max, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  abline(v=mean(TPRtest), col=2)
  abline(v=c(mean(TPRtest)-sd(TPRtest), mean(TPRtest)+sd(TPRtest)), col=4)
  legend("topleft", col=c(1,2,4), lwd=2, lty=1, cex=1.5, c(paste("TPR score"), 
                                                           paste("mean", format(mean(TPRtest), nsmall=6)), 
                                                           paste("sdev", format(sd(TPRtest), nsmall=6)) ) ) 
  
  hist(TNRtest, k.max, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  abline(v=mean(TNRtest), col=2)
  abline(v=c(mean(TNRtest)-sd(TNRtest), mean(TNRtest)+sd(TNRtest)), col=4)
  legend("topleft", col=c(1,2,4), lwd=2, lty=1, cex=1.5, c(paste("TNR score"), 
                                                           paste("mean", format(mean(TNRtest), nsmall=6)), 
                                                           paste("sdev", format(sd(TNRtest), nsmall=6)) ) ) 
  
  hist(ACCtest, k.max, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  abline(v=mean(ACCtest), col=2)
  abline(v=c(mean(ACCtest)-sd(ACCtest), mean(ACCtest)+sd(ACCtest)), col=4)
  legend("topleft", col=c(1,2,4), lwd=2, lty=1, cex=1.5, c(paste("ACC score"), 
                                                           paste("mean", format(mean(ACCtest), nsmall=6)), 
                                                           paste("sdev", format(sd(ACCtest), nsmall=6)) ) ) 
  
  hist(l2.score, k.max, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  abline(v=mean(l2.score), col=2)
  abline(v=c(mean(l2.score)-sd(l2.score), mean(l2.score)+sd(l2.score)), col=4)
  legend("topright", col=c(1,2,4), lwd=2, lty=1, cex=1.5, c(paste("L2 score"), 
                                                            paste("mean", format(mean(l2.score), nsmall=6)), 
                                                            paste("sdev", format(sd(l2.score), nsmall=6)) ) ) 
  
  dev.off()
  
}



# The next function loads variables from local repo

my.read.data <- function(dataset_name="gra", dataset_id="1", path_private="C:/Users/gianl/Desktop/Gianluca/"){
  
  if (dataset_name=="icn") {
    
    
    path_to_data <- paste(path_private, "ionchannel/", dataset_id, "/traj.dat", sep="")
    path_to_perm <- paste(path_private, "ionchannel/", dataset_id, "/perm.dat", sep="")
    
    data.read <- as.matrix(read.table(path_to_data,header=F))
    res.read <- mean(diff(data.read[,1]))
    sample.read <- nrow(data.read)
    
    col.read <- list("2"=9745, "3"=7741, "5"=9745, "6"=7741)
    X.read <- data.read[-1, 2:col.read[[dataset_id]]]
    ION.read <- data.read[-1,(col.read[[dataset_id]]+1):ncol(data.read)]
    
    ind.read <- as.matrix(read.table(path_to_perm,header=F))[,1]/res.read
    y.events.read <- rep(0, sample.read-1)
    y.events.read[ind.read] <- 1
    
  }
  
  if (dataset_name=="gra") {
    
    path_to_data <- paste(path_private, "gramicidin_ext/gramicidin", dataset_id, ".dat", sep="")
    path_to_perm <- paste(path_private, "gramicidin_ext/k_out", dataset_id, ".dat", sep="")
    
    data.read <- as.matrix(read.table(path_to_data,header=F))
    X.read <- data.read[,2:1657]
    ION.read <- data.read[,1658:ncol(data.read)]
    
    total.read <- 2000000
    sample.read <- nrow(X.read)-1
    ll.read <- readLines(path_to_perm)
    ind.read <- NULL
    for (line in ll.read) {
      if (grepl("permeations at times (ps): ", line, fixed=TRUE)==TRUE) {
        words <- scan(text=line, what="")
        index <- which(words=="(ps):")
        ind.read <- c(ind.read, as.integer(words[10:length(words)]) )
      }
    }
    ind.read <- sample.read*sort(ind.read)/total.read
    y.events.read <- rep(0, sample.read)
    y.events.read[ind.read] <- 1
    
  }
  
  return(list(X=X.read, ION=ION.read, Y=y.events.read))
  
}



# The next function simulates toy latent data

my.simul.data <- function(seed_id, dataset_name="sim",
                          n=500, p=100, m=5, sigma.p=2, sigma.e=0.25,
                          ifplot=TRUE){
  
  dataset_id <- as.character(seed_id)
  print(paste("Simulating dataset ", dataset_name, dataset_id, sep=""))
  print(paste("Using parameters: ",
              "n = ", n, "; ",
              "p = ", p, "; ",
              "m = ", m, "; ",
              "sigma.p = ", sigma.p, "; ",
              "sigma.e = ", sigma.e, "; ", sep=""))
  
  dataname.sim <- paste(dataset_name, dataset_id, sep="")
  filename.sim <- paste("_n", n, "_p", p, "_m", m, 
                        "_sp", sigma.p, "_se", sigma.e, sep="")
  
  # simulate latent design Z
  print("Simulating latent design Z")
  set.seed(seed_id)
  Z.sim <- matrix(0, nrow=n, ncol=m)
  if (dataset_name=="sim") {
    for (col in 1:m) {
      Z.sim[,col] <- rnorm(n, mean=0, sd=1)
    }
  } else if (dataset_name=="simdyn") {
    for (row in 2:n) {
      Z.sim[row,] <- Z.sim[row-1,] + rnorm(m, mean=0, sd=1)
    }
    for (col in 1:m) {
      Z.sim[,col] <- Z.sim[,col] - mean(Z.sim[,col])
      Z.sim[,col] <- Z.sim[,col]/max(abs(Z.sim[,col]))
    }
  }
  if (ifplot) {
    print("Plotting simulated latent design")
    pdf(file=paste(dataname.sim, "_Z.sim", filename.sim, ".pdf", sep=""))
    plot(Z.sim[,1], main="Simulated latent design's trajectory") 
    dev.off()
  }
  
  # simulate loading P
  print("Simulating loadings P")
  set.seed(seed_id+1)
  P.sim <- matrix(0, nrow=p, ncol=m)
  for (row in 1:p) {
    P.sim[row,] <- rnorm(m, mean=0, sd=sigma.p)
  }
  print(paste("Conditioning number of P is", kappa(t(P.sim)%*%P.sim)))
  
  # simulate design X = ZP' + E without centering
  print("Simulating design X = ZP' + E")
  set.seed(seed_id+2)
  X.sim <- Z.sim%*%t(P.sim)
  for (col in 1:p) {
    X.sim[,col] <- X.sim[,col] + rnorm(n, mean=0, sd=sigma.e)
    #X.sim[,col] <- X.sim[,col] - mean(X.sim[,col])
  }
  if (ifplot) {
    print("Plotting simulated design")
    pdf(file=paste(dataname.sim, "_X.sim", filename.sim, ".pdf", sep=""))
    plot(X.sim[,1], main="Simulated design's trajectory") 
    dev.off()
  }
  print(paste("Conditioning number of X is", kappa(t(X.sim)%*%X.sim)))
  
  # simulate latent coeffs Alpha
  print("Simulate latent coefficients Alpha")
  set.seed(seed_id+3)
  alpha.sim <- runif(m,5,10)
  if (ifplot) { 
    print("Plotting simulated latent coefficients")
    pdf(file=paste(dataname.sim, "_alpha.sim", filename.sim, ".pdf", sep=""))
    plot(alpha.sim, main="Simulated latent coefficients") 
    dev.off()
  }
  
  # compute latent probabilities Prob
  print("Compute latent probabilities Prob from linear link Z*Alpha")
  eta.sim <- Z.sim%*%alpha.sim
  prob.sim <- 1/(1+exp(-eta.sim))
  if (ifplot) { 
    print("Plotting simulated latent probabilities")
    pdf(file=paste(dataname.sim, "_prob.sim", filename.sim, ".pdf", sep=""))
    plot(prob.sim, main="Simulated latent probabilities") 
    dev.off()
  }
  
  # simulate binary response Y
  print("Simulate binary response Y from latent probabilities")
  set.seed(seed_id+4)
  y.sim <- rep(0,n)
  for (k in 1:n) {
    y.sim[k] <- rbinom(1, size=1, prob=prob.sim[k])
  }
  if (ifplot) { 
    print("Plotting simulated binary response")
    pdf(file=paste(dataname.sim, "_y.sim", filename.sim, ".pdf", sep=""))
    plot(y.sim, main="Simulated binary response") 
    dev.off()
  }
  
  # compute coeffs Beta = P*inv(P'P)*Alpha
  print("Compute coefficients Beta = P*inv(P'P)*Alpha")
  beta.sim <- P.sim%*%solve(t(P.sim)%*%P.sim)%*%alpha.sim
  if (ifplot) { 
    print("Plotting simulated coefficients")
    pdf(file=paste(dataname.sim, "_beta.sim", filename.sim, ".pdf", sep=""))
    plot(beta.sim, main="Simulated coefficients") 
    dev.off()
  }
  
  return(list(X=X.sim, 
              Y=y.sim, 
              beta=beta.sim, 
              alpha=alpha.sim, 
              prob=prob.sim,
              P=P.sim,
              Z=Z.sim))
  
}


# The next function simulates toy latent data
# It allows to use an input latent variable Z
# If Z is NULL, the function my.simul.data() is called instead
my.simul.data2 <- function(seed_id, dataset_name="sim",
                           n=500, p=100, m=5, sigma.p=2, sigma.e=0.25, Z=NULL,
                           ifplot=TRUE) {
  
  if (is.null(Z)) {
    
    return(my.simul.data(seed_id, dataset_name, n, p, m, sigma.p, sigma.e, ifplot))
    
  } else {
    
    dataset_id <- as.character(seed_id)
    dataname.sim <- paste(dataset_name, dataset_id, sep="")
    filename.sim <- paste("_n", n, "_p", p, "_m", m, 
                          "_sp", sigma.p, "_se", sigma.e, sep="")
    
    # if Z is not NULL we select the corresponding rows and columns
    # we center column-wise and then proceed with the scheme in my.simul.data()
    Z.sim <- as.matrix(Z[1:n,1:m])
    Z.sim <- scale(Z.sim, center=TRUE, scale=FALSE)
    if (ifplot) {
      print("Plotting simulated latent design")
      pdf(file=paste(dataname.sim, "_Z.sim", filename.sim, ".pdf", sep=""), height=12, width=24)
      plot(Z.sim[,1], main="Simulated latent design's trajectory")
      dev.off()
    }
    
    # simulate loading P
    print("Simulating loadings P")
    set.seed(seed_id+1)
    P.sim <- matrix(0, nrow=p, ncol=m)
    for (row in 1:p) {
      P.sim[row,] <- rnorm(m, mean=0, sd=sigma.p)
    }
    print(paste("Conditioning number of P is", kappa(t(P.sim)%*%P.sim)))
    
    # simulate design X = ZP' + E without centering
    print("Simulating design X = ZP' + E")
    set.seed(seed_id+2)
    X.sim <- Z.sim%*%t(P.sim)
    for (col in 1:p) {
      X.sim[,col] <- X.sim[,col] + rnorm(n, mean=0, sd=sigma.e)
      #X.sim[,col] <- X.sim[,col] - mean(X.sim[,col])
    }
    if (ifplot) {
      print("Plotting simulated design")
      pdf(file=paste(dataname.sim, "_X.sim", filename.sim, ".pdf", sep=""), height=12, width=24)
      plot(X.sim[,1], main="Simulated design's trajectory")
      dev.off()
    }
    print(paste("Conditioning number of X is", kappa(t(X.sim)%*%X.sim)))
    
    # simulate latent coeffs Alpha
    print("Simulate latent coefficients Alpha")
    set.seed(seed_id+3)
    alpha.sim <- runif(m,5,10)
    if (ifplot) { 
      print("Plotting simulated latent coefficients")
      pdf(file=paste(dataname.sim, "_alpha.sim", filename.sim, ".pdf", sep=""))
      plot(alpha.sim, main="Simulated latent coefficients")
      dev.off()
    }
    
    # compute latent probabilities Prob
    print("Compute latent probabilities Prob from linear link Z*Alpha")
    eta.sim <- Z.sim%*%alpha.sim
    prob.sim <- 1/(1+exp(-eta.sim))
    if (ifplot) { 
      print("Plotting simulated latent probabilities")
      pdf(file=paste(dataname.sim, "_prob.sim", filename.sim, ".pdf", sep=""), height=12, width=24)
      plot(prob.sim, main="Simulated latent probabilities")
      dev.off()
    }
    
    # simulate binary response Y
    print("Simulate binary response Y from latent probabilities")
    set.seed(seed_id+4)
    y.sim <- rep(0,n)
    for (k in 1:n) {
      y.sim[k] <- rbinom(1, size=1, prob=prob.sim[k])
    }
    if (ifplot) { 
      print("Plotting simulated binary response")
      pdf(file=paste(dataname.sim, "_y.sim", filename.sim, ".pdf", sep=""), height=12, width=24)
      plot(y.sim, main="Simulated binary response") 
      dev.off()
    }
    
    # compute coeffs Beta = P*inv(P'P)*Alpha
    print("Compute coefficients Beta = P*inv(P'P)*Alpha")
    beta.sim <- P.sim%*%solve(t(P.sim)%*%P.sim)%*%alpha.sim
    if (ifplot) { 
      print("Plotting simulated coefficients")
      pdf(file=paste(dataname.sim, "_beta.sim", filename.sim, ".pdf", sep=""))
      plot(beta.sim, main="Simulated coefficients") 
      dev.off()
    }
    
    return(list(X=X.sim, 
                Y=y.sim, 
                beta=beta.sim, 
                alpha=alpha.sim, 
                prob=prob.sim,
                P=P.sim,
                Z=Z.sim))
    
    
    
  }
}



# The next function simulates a toy general latent data
# In this version all quantities are deterministic

my.simul.data3 <- function(seed_id, dataset_name="sim",
                           n=500, p=100, m=5, sd.ex=0.25, alpha=NULL,
                           ifplot=TRUE){
  
  dataset_id <- as.character(seed_id)
  print(paste("Simulating dataset ", dataset_name, dataset_id, sep=""))
  print(paste("Using parameters: ",
              "n = ", n, "; ",
              "p = ", p, "; ",
              "m = ", m, "; ", 
              sep=""))
  
  dataname.sim <- paste(dataset_name, dataset_id, sep="")
  filename.sim <- paste("_n", n, 
                        "_p", p, 
                        "_m", m,
                        sep="")
  
  # simulate latent design Z = M.Z + Z.c(Sigma.Z)^1/2 
  print("Simulating latent design Z")
  set.seed(seed_id)
  
  m.Z <- 1:m                          # latent mean vector
  M.Z <- matrix(0, nrow=n, ncol=m)    # latent mean matrix
  Sigma.Z <- diag(1:m)                # latent covariance matrix
  Z.c <- matrix(0, nrow=n, ncol=m)    # latent normalized matrix
  for (col in 1:m) {
    M.Z[,col] <- m.Z[col]
    Z.c[,col] <- rnorm(n, mean=0, sd=1)
  }
  Z.sim <- M.Z + Z.c%*%sqrt(Sigma.Z)  # latent design 
  
  if (ifplot) {
    print("Plotting simulated latent design")
    pdf(file=paste(dataname.sim, "_Z", filename.sim, ".pdf", sep=""))
    plot(Z.c[,1], main="Simulated Z.c")
    plot(Z.sim[,1], main="Simulated Z") 
    dev.off()
  }
  
  # simulate loading P
  print("Simulating loadings P")
  set.seed(seed_id+1)
  
  P.sim <- matrix(0, nrow=p, ncol=m)
  for (col in 1:m) {
    P.sim[,col] <- rnorm(p, mean=0, sd=1)
  }
  print(paste("Conditioning number of P is", kappa(t(P.sim)%*%P.sim)))
  
  # simulate design X = M.X + X.c(Sigma.X)^1/2 + E.X with X.c = Z.c P'
  print("Simulating design X")
  set.seed(seed_id+2)
  
  m.X <- 1:p                          # mean vector
  M.X <- matrix(0, nrow=n, ncol=p)    # mean matrix
  Sigma.X <- diag(1:p)                # covariance matrix
  E.X <- matrix(0, nrow=n, ncol=p)    # residual matrix
  X.c <- Z.c%*%t(P.sim)               # normalized matrix
  for (col in 1:p) {
    M.X[,col] <- m.X[col]
    E.X[,col] <- rnorm(p, mean=0, sd=sd.ex)
  }
  X.sim <- M.X + X.c%*%sqrt(Sigma.X) + E.X # design 
  
  if (ifplot) {
    print("Plotting simulated design")
    pdf(file=paste(dataname.sim, "_X", filename.sim, ".pdf", sep=""))
    plot(X.c[,1], main="Simulated X.c")
    plot(X.sim[,1], main="Simulated X")
    dev.off()
  }
  print(paste("Conditioning number of X is", kappa(t(X.sim)%*%X.sim)))
  
  # simulate latent coeffs Alpha 
  print("Simulate latent coefficients")
  set.seed(seed_id+3)
  
  if (is.null(alpha)) alpha.c <- c(-5, rep(5, m))                             # coeffs Z.c
  
  alpha.sim <- rep(0, m+1)                                                    # from Z.c to Z
  alpha.sim[1] <- alpha.c[1]-m.Z%*%sqrt(solve(Sigma.Z))%*%(alpha.c[2:(m+1)])  # intercept Z
  alpha.sim[2:(m+1)] <- sqrt(solve(Sigma.Z))%*%(alpha.c[2:(m+1)])             # coeffs Z
  
  
  if (ifplot) { 
    print("Plotting simulated latent coefficients")
    pdf(file=paste(dataname.sim, "_alpha", filename.sim, ".pdf", sep=""))
    plot(alpha.c, main="Simulated coefficients for Z.c") 
    plot(alpha.sim, main="Simulated coefficients for Z") 
    dev.off()
  }
  
  # compute coeffs Beta 
  print("Compute coefficients")
  
  beta.c <- rep(0, p+1)                                                       # from Z.c to X.c
  beta.c[1] <- alpha.c[1]                                                     # intercept X.c
  beta.c[2:(p+1)] <- P.sim%*%solve(t(P.sim)%*%P.sim)%*%(alpha.c[2:(m+1)])     # coeffs X.c
  
  beta.sim <- rep(0, p+1)                                                     # from X.c to X
  beta.sim[1] <- beta.c[1]-m.X%*%sqrt(solve(Sigma.X))%*%(beta.c[2:(p+1)])     # intercept X
  beta.sim[2:(p+1)] <- sqrt(solve(Sigma.X))%*%(beta.c[2:(p+1)])               # coeffs X
  
  if (ifplot) { 
    print("Plotting simulated coefficients")
    pdf(file=paste(dataname.sim, "_beta", filename.sim, ".pdf", sep=""))
    plot(beta.c, main="Simulated coefficients for X.c") 
    plot(beta.sim, main="Simulated coefficients for X") 
    dev.off()
  }
  
  # compute latent probabilities Prob
  print("Compute latent probabilities Prob from linear link")
  
  eta.c <- cbind(1, Z.c)%*%alpha.c
  prob.c <- 1/(1+exp(-eta.c))
  
  eta.sim <- cbind(1, Z.sim)%*%alpha.sim
  prob.sim <- 1/(1+exp(-eta.sim))
  
  if (ifplot) { 
    print("Plotting simulated latent probabilities")
    pdf(file=paste(dataname.sim, "_prob", filename.sim, ".pdf", sep=""))
    plot(prob.c, main="Simulated latent probabilities from Z.c")
    plot(prob.sim, main="Simulated latent probabilities from Z") 
    dev.off()
  }
  
  # simulate binary response Y
  print("Simulate binary response Y from latent probabilities")
  set.seed(seed_id+4)
  
  Y.c <- rep(0,n)
  for (k in 1:n) {
    Y.c[k] <- rbinom(1, size=1, prob=prob.c[k])
  }
  Y.sim <- Y.c
  
  if (ifplot) { 
    print("Plotting simulated binary response")
    pdf(file=paste(dataname.sim, "_Y", filename.sim, ".pdf", sep=""))
    plot(Y.sim, main="Simulated binary response from Z.c") 
    plot(Y.sim, main="Simulated binary response from Z") 
    dev.off()
  }
  
  return(list(MZ=M.Z, Zc=Z.c, SigmaZ=Sigma.Z, Z=Z.sim,
              MX=M.X, Xc=X.c, SigmaX=Sigma.X, X=X.sim, EX=E.X,
              P=P.sim,
              Ac=alpha.c, A=alpha.sim,
              Bc=beta.c, B=beta.sim,
              Prc=prob.c, Pr=prob.sim, 
              Yc=Y.c, Y=Y.sim
              ))
  
}



# The next function is an auxiliary check of unidimensionality

my.is.unidim <- function(v){
  V <- as.matrix(v)
  return(dim(V)[2]==1)
}



# The next function subtracts a vector from the column of a matrix

my.vect.from.col <- function(x, v){
  X <- as.matrix(x)
  for (col in 1:length(v)) {
    X[,col] <- X[,col] - v[col]
  }
  return(X)
}



# The next function does pre-processing based on variance

my.preproc <- function(X, to_keep) {
  
  X <- as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  p3 <- p/3
  
  v.prep <- rep(0,p)
  
  for (col in 1:p) {
    v.prep[col] <- sd(X[,col])
  }
  
  v.level <- sort(v.prep, decreasing=TRUE)[to_keep]
  
  X.prep <- NULL
  
  for (col in 1:p) {
    if(v.prep[col]>=v.level) { X.prep <- cbind(X.prep, X[,col]) }
  }
  
  return(X.prep)
  
}



# The next function lags the rows of an input matrix

my.lag <- function(x, s.lag=0){
  
  X <- as.matrix(x)
  n <- dim(X)[1]
  
  time.lag <- (s.lag+1):n
  X.lag <- X[time.lag,]
  
  s <- 1
  while(s.lag>=s) {
    
    X.lag <- cbind(X.lag, X[time.lag-s,])
    s <- s+1
    
  }
  
  return(X.lag)
}



# The next function undoes the symmetry in ionchannel data
my.undo.symmetry <- function(x, 
                             center=c(4.1763,4.1405,0), 
                             r1=diag(c(1,1,1)), 
                             r2=diag(c(-1,-1,1)),  
                             r3=t(cbind(c(0,-1,0), c(1,0,0), c(0,0,1))), 
                             r4=t(cbind(c(0,1,0), c(-1,0,0), c(0,0,1))) )  {
  # Starting point
  x.new <- x
  
  # Useful quantities
  LL <- ncol(x.new)
  LL3 <- LL/3
  LL4 <- LL/4
  LL12 <- LL/12
  
  # Center the system around center
  x.new[,(1:LL3)*3-2] <- x.new[,(1:LL3)*3-2] - center[1]
  x.new[,(1:LL3)*3-1] <- x.new[,(1:LL3)*3-1] - center[2]
  x.new[,(1:LL3)*3] <- x.new[,(1:LL3)*3] - center[3]
  
  # Transform the first fold with r1
  fold <- 1
  R1 <- diag(LL)
  
  for (k in 1:LL12) {
    chunk <- ((fold-1)*LL4+3*(k-1)+1):(((fold-1)*LL4+3*(k-1)+3))
    R1[chunk, chunk] <- r1
  }
  x.new[,((fold-1)*LL4+1):(fold*LL4)] <- x.new[,((fold-1)*LL4+1):(fold*LL4)]%*%R1[((fold-1)*LL4+1):(fold*LL4),((fold-1)*LL4+1):(fold*LL4)]
  
  # Transform the second fold with r2
  fold <- 2
  R2 <- diag(LL)
  
  for (k in 1:LL12) {
    chunk <- ((fold-1)*LL4+3*(k-1)+1):(((fold-1)*LL4+3*(k-1)+3))
    R2[chunk, chunk] <- r2
  }
  x.new[,((fold-1)*LL4+1):(fold*LL4)] <- x.new[,((fold-1)*LL4+1):(fold*LL4)]%*%R2[((fold-1)*LL4+1):(fold*LL4),((fold-1)*LL4+1):(fold*LL4)]
  
  # Transform the third fold with r3
  fold <- 3
  R3 <- diag(LL)
  
  for (k in 1:LL12) {
    chunk <- ((fold-1)*LL4+3*(k-1)+1):(((fold-1)*LL4+3*(k-1)+3))
    R3[chunk, chunk] <- r3
  }
  x.new[,((fold-1)*LL4+1):(fold*LL4)] <- x.new[,((fold-1)*LL4+1):(fold*LL4)]%*%R3[((fold-1)*LL4+1):(fold*LL4),((fold-1)*LL4+1):(fold*LL4)]
  
  # Transform the fourth fold with r4
  fold <- 4
  R4 <- diag(LL)
  
  for (k in 1:LL12) {
    chunk <- ((fold-1)*LL4+3*(k-1)+1):(((fold-1)*LL4+3*(k-1)+3))
    R4[chunk, chunk] <- r4
  }
  x.new[,((fold-1)*LL4+1):(fold*LL4)] <- x.new[,((fold-1)*LL4+1):(fold*LL4)]%*%R4[((fold-1)*LL4+1):(fold*LL4),((fold-1)*LL4+1):(fold*LL4)]
  
  # Centering the system back
  x.new[,(1:LL3)*3-2] <- x.new[,(1:LL3)*3-2] + center[1]
  x.new[,(1:LL3)*3-1] <- x.new[,(1:LL3)*3-1] + center[2]
  x.new[,(1:LL3)*3] <- x.new[,(1:LL3)*3] + center[3]
  
  # visual
  plot(x[1], main="Before transformation")
  points((1:LL3*3-1), x[1,(1:LL3*3-1)], col=2)
  points((1:LL3*3), x[1,(1:LL3*3)], col=3)
  abline(v=c(LL4, 2*LL4, 3*LL4, 4*LL4), col=4)
  
  plot(x.new[1,], main="After transformation")
  points((1:LL3*3-1), x.new[1,(1:LL3*3-1)], col=2)
  points((1:LL3*3), x.new[1,(1:LL3*3)], col=3)
  abline(v=c(LL4, 2*LL4, 3*LL4, 4*LL4), col=4)
  
  return(x.new)
}




# The next function filters out the noise from a ion matrix.
# The x,y,z - coordinates are put to zero when the ion is not "trapped".
# A new multivariate permeation variable is computed accordingly.

my.ion.filter <- function(ION, Delta=25, theta=0.1){
  
  ION <- as.matrix(ION)
  
  # parameters
  p <- ncol(ION)
  p3 <- p/3
  n <- nrow(ION) - nrow(ION)%%Delta
  K <- n/Delta
  
  time.total <- 1:n
  
  # new variables
  ion <- ION[time.total,]
  ion.new <- 0*ion
  y.new <- matrix(0, nrow=n, ncol=p3)
  y.cross <- matrix(0, nrow=n, ncol=p3)
  y.trap <- matrix(0, nrow=n, ncol=p3)
  z.new <- rep(0,n)
  
  for (atom in 1:p3) {
    
    # center atom coordinates
    atom.x <- ion[,atom*3-2]
    atom.x <- atom.x - mean(atom.x)
    atom.y <- ion[,atom*3-1]
    atom.y <- atom.y - mean(atom.y)
    atom.z <- ion[,atom*3]
    atom.z <- atom.z - mean(atom.z)
    
    # compute I.xy and ion.new
    I.xy <- NULL
    for (k in 1:K) {
      
      local <- ((k-1)*Delta+1):(k*Delta)
      
      m.x <- sum(atom.x[local])/Delta
      v.x <- sum(atom.x[local]^2)/(Delta-1)
      m.y <- sum(atom.y[local])/Delta
      v.y <- sum(atom.y[local]^2)/(Delta-1)
      
      if (v.x<theta & abs(m.x)<theta & v.y<theta & abs(m.y)<theta) {
        
        I.xy <- c(I.xy, local)
        
        ion.new[local,atom*3-2] <- atom.x[local]
        ion.new[local,atom*3-1] <- atom.y[local]
        ion.new[local,atom*3] <- atom.z[local]
        
      }
      
    }
    
    # compute responses binary responses
    if (!is.null(I.xy)) {
      
      B.xy <- split(I.xy, cumsum(c(1,diff(I.xy) != 1)))
      
      for (k in 1:length(B.xy)) {
        
        I.xy.k <- B.xy[[k]]
        
        local.min <- min(I.xy.k)
        local.max <- max(I.xy.k)
        
        # y.trap has ones at first and last times when x,y-coords are trapped
        y.trap[local.min, atom] <- 1
        y.trap[local.max, atom] <- 1
        
        m.z.min <- mean(atom.z[local.min:(local.min+5)])
        m.z.max <- mean(atom.z[(local.max-5):local.max])
        
        if (sign(m.z.min)*sign(m.z.max) < 0) {
          
          # y.new has ones at exit times if z-coord has crossing
          y.new[max(I.xy.k), atom] <- 1
          
        }
        
        # y.cross has ones at crossing times of z-coord when x,y-coords are trapped
        J.xy.k <- I.xy.k[which(abs(atom.z[I.xy.k])>theta)]
        J.xy.k.cross <- J.xy.k[1+which(diff(sign(atom.z[I.xy.k]))!=0)]
        y.cross[J.xy.k.cross, atom] <- 1
        
      }
    }
    
    # compute z.new
    if (!is.null(I.xy)) {
      
      B.xy <- split(I.xy, cumsum(c(1,diff(I.xy) != 1)))
      
      for (k in 1:length(B.xy)) {
        
        I.xy.k <- B.xy[[k]]
        local.min <- min(I.xy.k)
        local.max <- max(I.xy.k)
        
        m.z.min <- sum(atom.z[local.min:(local.min+5)])/5
        m.z.max <- sum(atom.z[(local.max-5):local.max])/5
        
        if (sign(m.z.min)*sign(m.z.max) < 0) {
          
          z.new[I.xy.k] <- atom.z[I.xy.k]
          
        }
      }
    }
    
  }
  
  return(list(ion=ion.new, y=y.new, yc=y.cross, yt=y.trap, z=z.new))
}





# The next function takes a vector X and centers according to a parameter
# "none" performs no centering

my.center1 <- function(X, centering="none"){
  
  if (centering=="none"){
    X.center <- X
  } 
  else if (centering=="rigid"){
    X.center <- X-mean(X)
  } 
  else if (centering=="column"){
    X.center <- X-mean(X)
  }
  else {
    X.center <- X
  }
  return(X.center)
}


# The next function takes a matrix X and centers according to a parameter
# "rigid" centers the x,y,z coordinates by their global mean
# "column" centers X column-wise
# "none" performs no centering

my.center2 <- function(X, centering="none"){
  
  if (my.is.unidim(X)){
    X.center <- my.center1(X, centering)
  }
  else {
    n <- nrow(X)
    p <- ncol(X) # assumed to be a multiple of 3
    p3 <- p/3
    X.center <- X
    
    if (centering=="none"){
      
      X.center <- X
      
    } 
    else if (centering=="rigid"){
      
      x.mean <- 0
      y.mean <- 0
      z.mean <- 0
      
      for (col in 1:p3){
        x.mean <- x.mean + mean(X[,col*3-2])/p3
        y.mean <- y.mean + mean(X[,col*3-1])/p3
        z.mean <- z.mean + mean(X[,col*3-0])/p3
      }
      X.center[,(1:p3)*3-2] <- X.center[,(1:p3)*3-2] - x.mean
      X.center[,(1:p3)*3-1] <- X.center[,(1:p3)*3-1] - y.mean
      X.center[,(1:p3)*3-0] <- X.center[,(1:p3)*3-0] - z.mean
      
    } 
    else if (centering=="column"){
      
      for (col in 1:p3){
        X.center[,col*3-2] <- X.center[,col*3-2] - mean(X.center[,col*3-2])
        X.center[,col*3-1] <- X.center[,col*3-1] - mean(X.center[,col*3-1])
        X.center[,col*3-0] <- X.center[,col*3-0] - mean(X.center[,col*3-0])
      }
      
    } 
    else {
      
      X.center <- X
      
    }
  }
  return(X.center)
}



# The next function implements Partial Least Squares on
# univariate response Y and design X using NIPALS

my.pls1 <- function(X, Y, ncomp, maxit=50, tol=1e-8){
  
  Xh <- as.matrix(X)
  Yh <- as.matrix(Y)
  
  T <- NULL
  W <- NULL
  Q <- NULL
  U <- NULL
  P <- NULL
  D <- NULL
  C <- NULL
  W <- NULL
  
  for (h in 1:ncomp) {
    
    nr <- 0   # number of iterations
    uh <- Yh  # initialize residual of response Y
    end_algo <- FALSE   # stop condition
    
    # Inner ster of PLS2 NIPALS
    while (!end_algo) {
      
      nr <- nr + 1  # increase iterations
      
      wh <- t(Xh)%*%uh                      # this matches t(X)%*%Y
      wh <- wh/as.vector(sqrt(t(wh)%*%wh))  
      
      th <- Xh%*%wh   # projection of X into T
      
      ch <- t(Yh)%*%th                      # this matches t(Y)%*%T
      ch <- ch/as.vector(sqrt(t(ch)%*%ch))  
      
      uhnew <- Yh%*%ch  # projection of Y into U
      
      uhnew <- Yh%*%ch
      deltau <- uhnew-uh
      unorm <- as.numeric(sqrt(t(deltau)%*%deltau))
      
      if (unorm<tol | is.na(unorm)){
        end_algo <- TRUE
      }
      else if (nr > maxit) {
        end_algo <- TRUE
      }
      
      uh <- uhnew
      
    }
    
    ph <- t(Xh)%*%th/as.vector(t(th)%*%th)  # p scores
    qh <- t(Yh)%*%uh/as.vector(t(uh)%*%uh)  # Q scores
    dh <- t(uh)%*%uh/as.vector(t(th)%*%th)  # D scores
    
    Xh <- Xh - th%*%t(ph)                   # X deflation
    Yh <- Yh - (th%*%t(ch)*as.vector(dh))   # Y deflation
    
    T <- cbind(T, th)
    Q <- cbind(Q, qh)
    U <- cbind(U, uh)
    P <- cbind(P, ph)
    D <- cbind(D, dh)
    C <- cbind(C, ch)
    W <- cbind(W, wh)
    
  }
  
  R <- W%*%solve(t(P)%*%W)
  B <- R%*%t(C)
  
  return(list(P=P, T=T, Q=Q, U=U, D=D, W=W, C=C, R=R, B=B))
}



# The next function implements Partial Least Squares on
# multivariate response Y and design X using NIPALS

my.pls2 <- function(X, Y, ncomp, maxit=50, tol=1e-8){
  
  if (my.is.unidim(Y)){
    
    return(my.pls1(X, Y, ncomp, maxit, tol))
    
  }
  else {
    
    Xh <- as.matrix(X) 
    Yh <- as.matrix(Y)
    
    T <- NULL
    W <- NULL
    Q <- NULL
    U <- NULL
    P <- NULL
    D <- NULL
    C <- NULL
    W <- NULL
    
    for (h in 1:ncomp) {
      
      nr <- 0
      uh <- Yh[,1]
      end_algo <- FALSE
      
      # Inner ster of PLS2 NIPALS
      while (!end_algo) {
        
        nr <- nr + 1
        
        wh <- t(Xh)%*%uh
        wh <- wh/as.vector(sqrt(t(wh)%*%wh))
        
        th <- Xh%*%wh
        
        ch <- t(Yh)%*%th
        ch <- ch/as.vector(sqrt(t(ch)%*%ch))
        
        uhnew <- Yh%*%ch
        deltau <- uhnew-uh
        unorm <- as.numeric(sqrt(t(deltau)%*%deltau))
        
        if (unorm<tol | is.na(unorm)){
          end_algo <- TRUE
        }
        else if (nr > maxit) {
          end_algo <- TRUE
        }
        
        uh <- uhnew
        
      }
      
      ph <- t(Xh)%*%th/as.vector(t(th)%*%th)
      qh <- t(Yh)%*%uh/as.vector(t(uh)%*%uh)
      dh <- t(uh)%*%uh/as.vector(t(th)%*%th)
      
      Xh <- Xh - th%*%t(ph)
      Yh <- Yh - (th%*%t(ch)*as.vector(dh))
      
      T <- cbind(T, th)
      Q <- cbind(Q, qh)
      U <- cbind(U, uh)
      P <- cbind(P, ph)
      D <- cbind(D, dh)
      C <- cbind(C, ch)
      W <- cbind(W, wh)
      
    }
    
    R <- W%*%solve(t(P)%*%W)
    B <- R%*%t(C)
    
    return(list(P=P, T=T, Q=Q, U=U, D=D, W=W, C=C, R=R, B=B))
    
  }
  
}




# The next function implements Partial Least Squares with
# univariate response and Cross Validation for component selection
# The function splits the observations in blocks and use all
# but one to predict

my.pls1.cv <- function(X, Y, ncomp_max, nblock=10, maxit=50, tol=1e-8){
  
  n <- length(Y)
  p <- ncol(X)
  p3 <- p/3
  
  time.total <- 1:n
  delta <- n/nblock
  mean.score <- rep(0,ncomp_max)
  
  for (h in 1:ncomp_max) {
    
    for (k in 1:nblock) {
      
      time.local.test <- ((k-1)*delta+1):(k*delta)
      time.local.train <- setdiff(time.total, time.local.test)
      
      B.local.train <- my.pls1(X[time.local.train,], Y[time.local.train], h, maxit, tol)$B
      
      eta <- X%*%B.local.train
      score <- sum((Y[time.local.test]-eta[time.local.test])^2)/length(time.local.test)
      mean.score[h] <- mean.score[h] + score/nblock
      
    }
    print(paste("PLS comp", h, "has score", mean.score[h], sep=" "))
    
  }
  
  h.best <- which(mean.score==min(mean.score))
  return(my.pls1(X, Y, h.best, maxit, tol))
  
}



# The next function implements Partial Least Squares with
# multivariate response and Cross Validation for component selection
# The function splits the observations in blocks and use all
# but one to predict

my.pls2.cv <- function(X, Y, ncomp_max, nblock=10, maxit=50, tol=1e-8){
  
  if (my.is.unidim(Y)) {
    
    return(my.pls1.cv(X, Y, ncomp_max, nblock, maxit, tol))
    
  }
  else {
    
    n <- nrow(Y)
    p <- ncol(X)
    p3 <- p/3
    
    time.total <- 1:n
    delta <- n/nblock
    mean.score <- rep(0,ncomp_max)
    
    for (h in 1:ncomp_max) {
      
      for (k in 1:nblock) {
        
        time.local.test <- ((k-1)*delta+1):(k*delta)
        time.local.train <- setdiff(time.total, time.local.test)
        
        B.local.train <- my.pls2(X[time.local.train,], Y[time.local.train,], h, maxit, tol)$B
        
        eta <- X%*%B.local.train
        score <- (Y[time.local.test,]-eta[time.local.test])^2/length(time.local.test)
        score <- sum(apply(score,2,sum))/ncol(Y)
        mean.score[h] <- mean.score[h] + score/nblock
        
      }
      print(paste("PLS comp", h, "has score", mean.score[h], sep=" "))
      
    }
    
    h.best <- which(mean.score==min(mean.score))
    return(my.pls2(X, Y, h.best, maxit, tol))
    
  }
  
}




# The next function implements Partial Least Squares on
# uni/multivariate response Y and design X using KERNELPLS
# The function allows for centering/scaling (FALSE by default)
# The function allows for intercept (FALSE by default)

my.kernelpls2 <- function(X, Y, ncomp, 
                          centering=FALSE, scaling=FALSE, intercept=FALSE,
                          maxit=50, tol=1e-8){
  print("Performing KERNELPLS...")
  
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

my.plsr <- function(X, Y, ncomp, 
                    centering=FALSE, scaling=FALSE, intercept=FALSE,
                    maxit=50, tol=1e-8, method="kernel"){
  
  if (method=="nipals") {
    print("Performing PLS with nipals...")
    return(my.pls2(X, Y, ncomp, maxit, tol))
  }
  else if (method=="kernel") {
    print("Performing PLS with kernelpls...")
    return(my.kernelpls2(X=X, Y=Y, ncomp=ncomp, centering=centering, scaling=scaling, intercept=intercept, maxit=maxit, tol=tol))
  }
  else {
    print("Invalid method selection, performing kernelpls as default...")
    return(my.kernelpls2(X=X, Y=Y, ncomp=ncomp, centering=centering, scaling=scaling, intercept=intercept, maxit=maxit, tol=tol))
  }
  
}



# The next function implements IRLS with univariate Bernoulli response
# The parameter beta0 is the initialization, choosing NULL initializes at zero (or MLE)
# The function allows for centering/scaling (TRUE by default)
# The function allows for intercept (TRUE by default)

my.irls <- function(X, Y, beta0=NULL,
                    centering=TRUE, scaling=FALSE, intercept=TRUE,
                    lambda = 0,
                    maxit=50, tol=1e-6,
                    gof.check=FALSE){
  
  print("Performing IRLS")
  
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
    X <- cbind(1, X) # intercept is added as first column
  }
  
  # Get parameters
  n  <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(Y)[2]
  
  # Initialization
  beta <- beta0
  if (is.null(beta)){ beta <- rep(0,p) }
  Z <- NULL
  W <- NULL
  
  kappa1<-function(x) 1/(1+exp(-x))
  kappa2<-function(x) exp(-x)/(1+exp(-x))^2
  
  counter <- 0
  repeat{
    print(paste("Performing IRLS iter", counter, sep=" "))
    
    # Compute coefficients for X, Y
    Z.old <- Z
    W.old <- W
    beta[which(is.na(beta))] <- 0
    beta.old <- beta
    
    eta <- X%*%beta
    
    k.p <- kappa1(eta)
    k.p[which(k.p<0.01)] <- 0.01
    k.p[which(k.p>0.99)] <- 0.99
    
    k.pp <- kappa2(eta)
    k.pp[which(k.pp<0.01)] <- 0.01
    k.pp[which(k.pp>0.99)] <- 0.99
    
    W <- diag(as.vector(k.pp))
    
    Z <- eta + diag(1/as.vector(k.pp))%*%(Y-as.vector(k.p))
    
    # Weighted X.W and Z.W
    X.W <- as.matrix(sqrt(W)%*%X)
    Z.W <- as.matrix(sqrt(W)%*%Z)
    
    # Compute coefficients for Z.W ~ X.W
    # No center/scale X.W, Z.W
    # No intercept X.W
    beta <- solve(t(X.W)%*%X.W + diag(lambda,p))%*%t(X.W)%*%Z.W
    
    print(paste("Condition number is", kappa(t(X.W)%*%X.W + diag(lambda,p))))
    
    epsilon <- sqrt(sum((beta-beta.old)^2)/sum(beta.old^2))
    print(paste("Divergence", epsilon, sep=" "))
    
    like <- prod(kappa1(X%*%beta)^Y)*prod((1-kappa1(X%*%beta))^(1-Y))
    like.ratio <- like/(prod(kappa1(X%*%beta.old)^Y)*prod((1-kappa1(X%*%beta.old))^(1-Y)))
    GOF <- log(like.ratio)
    print(paste("Likelihood", like, sep=" "))
    print(paste("Likelihood ratio", like.ratio, sep=" "))
    print(paste("Goodness of fit", GOF, sep=" "))
    
    if (epsilon<tol) { 
      print("Divergence stop...")
      break 
    }
    
    if ((counter>0) & gof.check & (GOF<log(1+tol))) {
      print("GOF too small stop...")
      beta <- beta.old
      Z <- Z.old
      W <- W.old
      counter <- counter-1
      break
    }
    
    if(epsilon<tol) { 
      break 
    }
    
    if(counter==maxit) { 
      print("Maximum iterarion, no convergence...")
      break
    }
    else {
      counter <- counter+1
    }
    
  }
  
  # Compute coeffs for original variables
  beta.tilde <- beta*0
  
  if (intercept) {
    
    beta.tilde[1] <- sd.Y[1]*beta[1] + mu.Y[1] - sd.Y[1]*(mu.X/sd.X)%*%beta[2:p] # intercept
    beta.tilde[2:p] <- sd.Y[1]*beta[2:p]/sd.X                                    # core
    
  } else {
    
    beta.tilde[1:p] <- sd.Y[1]*beta[1:p]/sd.X                                    # no intercept
    
  }
  
  return(list(BETA=beta.tilde,
              Z=Z, W=W,
              it=counter))
}



# The next function implements RIRLS for Bernoulli response

my.rirls <- function(X, Y, lambdaR=0, Sigma=NULL, maxit=50, tol=1e-8){
  print("Performing RIRLS.")
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  n  <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(Y)[2]
  
  if (is.null(Sigma)){ Sigma <- diag(1,p) }
  
  beta <- rep(0,p)
  Z <- NULL
  W <- NULL
  
  counter <- 0
  repeat{
    print(paste("Performing RIRLS iter", counter, sep=" "))
    
    beta.old <- beta
    
    eta <- X%*%beta
    
    k.p <- 1/(1+exp(-eta))
    k.p <- k.p[,1]
    
    k.pp <- exp(-eta)/(1+exp(-eta))^2
    k.pp <- k.pp[,1]
    
    W <- diag(k.pp)
    
    if (min(k.pp) > tol) {
      
      Z <- eta + diag(1/k.pp)%*%(Y-k.p)
      
      # # adaptive choice of lambdaR
      # bic.lambda <- NULL
      # lambda.seq <- seq(0, 1000, by=10)
      # for (lambda in lambda.seq) {
      #   bic.lambda <- c(bic.lambda, -2*(t(X%*%beta)%*%Y -sum(log(1+exp(X%*%beta)))) +log(n)*(sum(diag(X%*%solve(t(X)%*%W%*%X+lambda*Sigma^2)%*%t(X)%*%W))))
      # }
      # lambdaR <- lambda.seq[which(bic.lambda==min(bic.lambda))]
      
      beta <- solve(t(X)%*%W%*%X + lambdaR*Sigma^2)%*%t(X)%*%W%*%Z
      
    }
    else { print("Estimated probabilities are too small.") }
    
    epsilon <- sqrt(sum((beta-beta.old)^2)/sum(beta.old^2))
    print(paste("Divergence", epsilon, sep=" "))
    
    if(epsilon<=tol) { break }
    
    if(counter==maxit) { print("Maximum iterarion, no convergence."); break }
    
    counter <- counter+1
    
  }
  
  return(list(beta=beta, Z=Z, W=W, lambda=lambdaR))
  
}




# The next function implements WPLS with
# univariate response

my.wpls1 <- function(X, Sigma, Y, W, ncomp, maxit=50, tol=1e-8){
  print("Performing WPLS1.")
  
  X <- as.matrix(X)
  Sigma <- as.matrix(Sigma)
  Y <- as.matrix(Y)
  W <- as.matrix(W)
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # initialize auxiliary
  t.k <- matrix(1, nrow=n, ncol=1)
  Sig.inv <- diag(1/diag(Sigma))
  E.0 <- X%*%Sig.inv # to output
  E.k <- X%*%Sig.inv
  f.k <- Y
  omega.k <- matrix(0, nrow=p, ncol=1)
  Psi <- diag(1, nrow=p)
  
  # initialize required
  Psi.tilde <- NULL
  q.m <- NULL
  B <- array(0, c(p, 1, ncomp+1)) 
  
  for (k in 0:ncomp){
    print(paste("Performing WPLS1 comp", k, "of", ncomp, sep=" "))
    
    norm <- (t(t.k)%*%W%*%t.k)[1,1]
    
    q.k <- ((t(t.k)%*%W%*%f.k)/norm)[1,1]
    q.m <- rbind(q.m, q.k)
    
    p.k <- (t(E.k)%*%W%*%t.k)/norm
    
    f.k <- f.k - q.k*t.k
    
    E.k <- E.k - t.k%*%t(p.k)
    
    Psi <- Psi%*%(diag(1,nrow=p)-omega.k%*%t(p.k))
    
    omega.k <- t(E.k)%*%W%*%f.k
    
    t.k <- E.k%*%omega.k
    
    Psi.tilde <- cbind(Psi.tilde, Psi%*%omega.k)
    
    B[,,k+1] <- Sig.inv%*%(Psi.tilde[,-1])%*%(q.m[-1])
    
  }
  
  Psi.tilde <- Psi.tilde[,-1] # the first column is zero, the remaining part is a (p,ncomp)-matrix
  
  q.m <- q.m[-1] # the first element is zero, the remaining part is a ncomp-vector
  
  return(list(Xscale=E.0, B=B[,,2:(ncomp+1)]))
  
}



# The next function implements PLSGLM with univariate response
# The parameter beta0 is the initialization, choosing NULL initializes at zero (or MLE)
# The function allows for centering/scaling (TRUE by default)
# The function allows for intercept (TRUE by default)

my.plsglm <- function(X, Y, ncomp, beta0=NULL,
                      centering=TRUE, scaling=TRUE, intercept=TRUE,
                      maxit=50, tol=1e-6,
                      gof.check=FALSE){
  
  print("Performing PLSGLM")
  
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
    X <- cbind(1, X) # intercept is added as first column
  }
  
  # Get parameters
  n  <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(Y)[2]
  
  # Initialization
  beta <- beta0
  if (is.null(beta)){ beta <- rep(0,p) }
  Z <- NULL
  W <- NULL
  
  kappa1<-function(x) 1/(1+exp(-x))
  kappa2<-function(x) exp(-x)/(1+exp(-x))^2
  
  history <- list(gof=NULL, likeratio=NULL, like=NULL, div=NULL)
  
  counter <- 0
  repeat{
    print(paste("Performing PLSGLM iter", counter, sep=" "))
    
    # Stores relevant quantities for GOF check
    if (counter>0) {
      fit.pls.old <- fit.pls
      Z.old <- Z
      W.old <- W
    }
    
    # Compute coefficients for X, Y
    beta[which(is.na(beta))] <- 0
    beta.old <- beta
    
    eta <- X%*%beta
    
    k.p <- kappa1(eta)
    k.p <- k.p[,1]
    
    k.pp <- kappa2(eta)
    k.pp <- k.pp[,1]
    k.pp[which(k.pp<1e-16)] <- 1e-16
    
    W <- diag(k.pp)
    
    if (sum(is.na(W))>0) {
      print("Estimated probabilities are NA...")
      
    } else {
      
      Z <- eta + diag(1/k.pp)%*%(Y-k.p)
      
      # Weighted X.W and Z.W
      X.W <- as.matrix(sqrt(W)%*%X)
      Z.W <- as.matrix(sqrt(W)%*%Z)
      
      # Compute coefficients for Z.W ~ X.W
      # No center/scale X.W, Z.W
      # No intercept X.W
      fit.pls <- my.plsr(X=X.W, Y=Z.W, ncomp=ncomp,
                         centering=FALSE, scaling=FALSE, intercept=FALSE,
                         maxit=50, tol=1e-8, method="kernel")
      R <- fit.pls$R
      beta.out <- fit.pls$Btilde
      beta <- beta.out[,,ncomp]
      
    }
    
    print(paste("Norm of beta", sum(beta^2), sep=" "))
    
    epsilon <- sqrt(sum((beta-beta.old)^2)/sum(beta.old^2))
    print(paste("Divergence", epsilon, sep=" "))
    
    like <- prod(kappa1(X%*%beta)^Y)*prod((1-kappa1(X%*%beta))^(1-Y))
    like.ratio <- like/(prod(kappa1(X%*%beta.old)^Y)*prod((1-kappa1(X%*%beta.old))^(1-Y)))
    GOF <- log(like.ratio)
    print(paste("Likelihood", like, sep=" "))
    print(paste("Likelihood ratio", like.ratio, sep=" "))
    print(paste("Goodness of fit", GOF, sep=" "))
    
    history$gof <- c(history$gof, GOF)
    history$likeratio <- c(history$likeratio, like.ratio)
    history$like <- c(history$like, like)
    history$div <- c(history$div, epsilon)
    
    if (epsilon<tol) { 
      print("Divergence stop...")
      break 
    }
    
    if ((counter>0) & gof.check & (GOF<log(1+tol))) {
      print("GOF too small stop...")
      beta.out <- fit.pls.old$Btilde
      beta <- beta.out[,,ncomp]
      Z <- Z.old
      W <- W.old
      R <- fit.pls.old$R
      counter <- counter-1
      break
    }
    
    if (counter==maxit) { 
      print("Maximum iterarion, no convergence...")
      break
    }
    else {
      counter <- counter+1
    }
    
  }
  
  # Compute coeffs for original variables
  beta.tilde <- beta.out*0
  
  for (m in 1:ncomp) {
    
    if (intercept) {
      
      beta.tilde[,,m][1] <- sd.Y[1]*beta.out[,,m][1] + mu.Y[1] - sd.Y[1]*(mu.X/sd.X)%*%beta.out[,,m][2:p] # intercept
      beta.tilde[,,m][2:p] <- sd.Y[1]*beta.out[,,m][2:p]/sd.X                                             # core
      
    } else {
      
      for (mc in 1:m) {
        beta.tilde[,,m][1:p] <- sd.Y[1]*beta.out[,,m][1:p]/sd.X                                           # no intercept, core
      }
    }
  }
  
  # Compute weight matrix
  W.out <- array(0, dim=c(n, n, ncomp))
  
  for (m in 1:ncomp) {
    W.out[,,m] <- diag(kappa2(X%*%beta.out[,,m])[,1])
  }
  
  return(list(BETA=beta.tilde, 
              beta=beta.out, 
              Z=Z, W=W, 
              R=R, 
              w=W.out,
              it=counter,
              history=history))
}




# The next function implements a variation of PLSGLM

my.plsglm.new <- function(X, Y, ncomp, beta0=NULL,
                          centering=TRUE, scaling=TRUE, intercept=TRUE,
                          maxit=50, tol=1e-6,
                          kappa1=function(x) {1-1/x-1/(1-exp(x))},
                          kappa2=function(x) {1/x^2+1/(2-2*cosh(x))},
                          verbose=FALSE){
  
  print("Performing PLSGLM-NEW")
  
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
    X <- cbind(1, X) # intercept is added as first column
  }
  
  # Get parameters
  n  <- dim(X)[1]
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
      k.p[,,m][which(k.p[,,m]<0.01)] <- 0.01
      k.p[,,m][which(k.p[,,m]>0.99)] <- 0.99
      
      k.pp[,,m] <- kappa2(eta[,,m])
      k.pp[,,m][which(k.pp[,,m]<0.01)] <- 0.01
      k.pp[,,m][which(k.pp[,,m]>0.99)] <- 0.99
      
      W[,,m] <- diag(as.vector(k.pp[,,m]))
      
      Z[,,m] <- eta[,,m] + diag(1/as.vector(k.pp[,,m]))%*%(Y-as.vector(k.p[,,m]))
      
      # Weighted X.W and Z.W on largest component
      X.W <- as.matrix(sqrt(W[,,m])%*%X)
      Z.W <- as.matrix(sqrt(W[,,m])%*%Z[,,m])
      
      # Compute coefficients for Z.W ~ X.W
      # No center/scale X.W, Z.W
      # No intercept X.W
      fit.pls <- my.plsr(X=X.W, Y=Z.W, ncomp=m,
                         centering=FALSE, scaling=FALSE, intercept=FALSE,
                         maxit=50, tol=1e-8, method="kernel")
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
  
  return(list(BETA=beta.tilde, 
              beta=beta, 
              Z=Z, W=W, 
              R=R,
              it=it.stop))
}






# The next function implements PLSGLM with Continuous Bernoulli response
# The parameter beta0 is the initialization, choosing NULL initializes at MLE 
# The function allows for centering/scaling (TRUE by default)
# The function allows for intercept (TRUE by default)

my.plsglm.cb <- function(X, Y, ncomp, beta0=NULL,
                         centering=TRUE, scaling=TRUE, intercept=TRUE,
                         maxit=50, tol=1e-6){
  
  print("Performing PLSGLM-CB")
  
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
    X <- cbind(1, X) # intercept is added as first column
  }
  
  # Get parameters
  n  <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(Y)[2]
  
  # Initialization
  beta <- beta0
  if (is.null(beta)){ 
    if (intercept) {
      beta <- lm(Y~X-1)$coefficients 
    } else {
      beta <- lm(Y~X-1)$coefficients
    }
  }
  Z <- NULL
  W <- NULL
  
  kappa1<-function(x) 1-1/x-1/(1-exp(x))
  kappa2<-function(x) 1/x^2+1/(2-2*cosh(x))
  
  counter <- 0
  repeat{
    print(paste("Performing PLSGLM-CB iter", counter, sep=" "))
    
    # Compute coefficients for X, Y
    beta[which(is.na(beta))] <- 0
    beta.old <- beta
    
    eta <- X%*%beta
    
    k.p <- kappa1(eta)
    k.p <- k.p[,1]
    
    k.pp <- kappa2(eta)
    k.pp <- k.pp[,1]
    k.pp[which(k.pp<1e-16)] <- 1e-16
    
    W <- diag(k.pp)
    
    if (sum(is.na(W))>0) {
      print("Estimated probabilities are NA...")
      
    } else {
      
      Z <- eta + diag(1/k.pp)%*%(Y-k.p)
      
      # Weighted X.W and Z.W
      X.W <- as.matrix(sqrt(W)%*%X)
      Z.W <- as.matrix(sqrt(W)%*%Z)
      
      # Compute coefficients for Z.W ~ X.W
      # No center/scale X.W, Z.W
      # No intercept X.W
      fit.pls <- my.plsr(X=X.W, Y=Z.W, ncomp=ncomp,
                         centering=FALSE, scaling=FALSE, intercept=FALSE,
                         maxit=50, tol=1e-8, method="kernel")
      R <- fit.pls$R
      beta.out <- fit.pls$Btilde
      beta <- beta.out[,,ncomp]
      
    }
    
    epsilon <- sqrt(sum((beta-beta.old)^2)/sum(beta.old^2))
    print(paste("Divergence", epsilon, sep=" "))
    
    if(epsilon<tol) { 
      break 
    }
    
    if(counter==maxit) { 
      print("Maximum iterarion, no convergence...")
      break
    }
    else {
      counter <- counter+1
    }
    
  }
  
  # Compute coeffs for original variables
  beta.tilde <- beta.out*0
  
  for (m in 1:ncomp) {
    
    if (intercept) {
      
      beta.tilde[,,m][1] <- sd.Y[1]*beta.out[,,m][1] + mu.Y[1] - sd.Y[1]*(mu.X/sd.X)%*%beta.out[,,m][2:p] # intercept
      beta.tilde[,,m][2:p] <- sd.Y[1]*beta.out[,,m][2:p]/sd.X                                             # core
      
    } else {
      
      for (mc in 1:m) {
        beta.tilde[,,m][1:p] <- sd.Y[1]*beta.out[,,m][1:p]/sd.X                                           # no intercept, core
      }
      
    }
    
  }
  
  # Compute weight matrix
  W.out <- array(0, dim=c(n, n, ncomp))
  
  for (m in 1:ncomp) {
    W.out[,,m] <- diag(kappa2(X%*%beta.out[,,m])[,1])
  }
  
  return(list(BETA=beta.tilde, 
              beta=beta.out, 
              Z=Z, W=W, 
              R=R, 
              w=W.out,
              it=counter))
  
}



# The next function implements PCA with PRCOMP
# Original code at https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/prcomp.R

my.pca <- function(x, 
                   centering=TRUE, scaling=TRUE, 
                   ncomp=ncol(x)) {
  
  print("Performing PCA")
  
  # Get variables
  X <- as.matrix(x)
  
  # Center variables (column-wise)
  mu.X <- rep(0, ncol(X))
  if (centering) {
    mu.X <- apply(X, 2, mean)
    X <- apply(X, 2, function(v){v-mean(v)})
  }
  
  # Scale variables (column-wise)
  sd.X <- rep(1, ncol(X))
  if (scaling) {
    sd.X <- apply(X, 2, function(v){max(sd(v), tol)})
    X <- apply(X, 2, function(v){v/max(sd(v), tol)})
  }
  
  # Get parameters
  n  <- dim(X)[1]
  p <- dim(X)[2]
  
  # Fit PCA
  fit.pca <- prcomp(x=X, center=FALSE, scale.=FALSE, rank.=ncomp, retx=FALSE)
  
  R <- fit.pca$rotation     # rotation matrix for center/scale X
  sdev <- fit.pca$sdev
  
  return(list(R=R, 
              muX=mu.X, sdX=sd.X,
              sdev=sdev))
}





# The next function implements TCCA as in [de Groot (2019)]

my.tcca <- function(x, ncomp, s.lag=0){
  
  X <- as.matrix(x)
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  time.0 <- 1:(n-s.lag)
  time.1 <- (s.lag+1):n
  
  # means
  print("Computing means...")
  mu.0 <- apply(X[time.0,], 2, mean)
  mu.1 <- apply(X[time.1,], 2, mean)
  
  # centered data
  print("Centering...")
  X.0 <- my.vect.from.col(X[time.0,], mu.0)
  X.1 <- my.vect.from.col(X[time.1,], mu.1)
  
  # instantaneous covariance matrices
  print("Computing covariances...")
  C.00 <- crossprod(X.0, X.0)/length(time.0)
  C.01 <- crossprod(X.0, X.1)/length(time.0)
  C.11 <- crossprod(X.1, X.1)/length(time.1)
  
  # Koopman matrix 
  print("Computing K...")
  K <- solve(sqrtm(C.00))%*%C.01%*%solve(sqrtm(C.01))
  
  # singular value decomposition of K
  print("Computing svd of K...")
  usv <- svd(K, nu=ncomp, nv=ncomp)
  S <- usv$d
  U <- usv$u
  V <- usv$v
  
  # left and right singular functions
  print("Computing functions...")
  f.0 <- t(U)%*%solve(sqrtm(C.00))%*%my.vect.from.col(X, mu.0)
  f.1 <- t(V)%*%solve(sqrtm(C.00))%*%my.vect.from.col(X, mu.1)
  
  return(list(X0=X.0, X1=X.1, mu0=mu.0, mu1=mu.1, C00=C.00, C01=C.01, C11=C.11, K=K, S=S, U=U, V=V, f0=f.0, f1=f.1))
  
}




# The next function implements D-GPR-PLS for univariate response
# Source in Liu et al. (2019) - Dynamic Nonlinear PLS Modeling Using GPR

my.dgprpls <- function(X, Y, ncomp, maxit=50, tol=1e-8, method="kernel") {
  return(TRUE)
}



