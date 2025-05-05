############################ Function to create the table with the results #####################


# Imports
library(dplyr)
library(writexl)

Create_table_results<-function(fn_coefs){


# Read data
dx.train <- as.data.frame(read.table("Data/XtrainNScaled.txt",header=TRUE))
dx.test <- as.data.frame(read.table("Data/XtestNScaled.txt",header=TRUE))
dy.train <- as.data.frame(read.table("Data/YTrain.txt",header=FALSE))
dy.test <- as.data.frame(read.table("Data/YTest.txt",header=FALSE))

# Response variable
y.all <- as.matrix(rbind(dy.train,dy.test))
n.all <- nrow(y.all)
tt <- (1:n.all)/n.all

# Add three artificial descriptors for additives
xc <- rbind(dx.train,dx.test)
set.seed(2134)
c1 <- factor(xc[,1])
c2 <- factor(xc[,3])
c3 <- factor(xc[,4])
levels(c1) <- runif(22)
levels(c2) <- rnorm(22)
levels(c3) <- runif(22)
xc <- cbind(cbind(as.numeric(c1),as.numeric(c2),as.numeric(c3)),xc)
colnames(xc)[1:3] <- c("add_new1","add_new2","add_new3")

xs <- xc[,c(4,23,50,60)]
colnames(xs) <- c("additive","aryl_halide","base","ligand")

a <- rep(NA,ncol(xs))
for (i in 1:ncol(xs)) a[i] <- length(unique(xs[,i]))
xcf <- matrix(NA,nrow(xs),sum(a))
b <- cumsum(a)
colnames(xcf) <- rep(colnames(xs),times=a)
colnum <- order(unique(xs[,1]))
for (i in 2:ncol(xs)) colnum <- c(colnum,order(unique(xs[,i])))
colnames(xcf) <- paste(colnames(xcf),colnum)
for (i in 1:nrow(xs)) {
  for (j in 1:length(a)) {
    res <- rep(0, a[j])
    where <- match( xs[i,j], unique(xs[,j]) )
    res[ where ] <- 1 
    xcf[i,(max(b[j-1],0)+1):b[j]] <- res
  }
}
for (i in 1:length(b)) {
  ind <- match(xcf[,b[i]],1)==1
  xcf[ind,(max(b[i-1],0)+1):b[i]] <- -1
}
xcf <- xcf[,-b]
x.all <- xcf

# Identify label ijkl for yield
rownames(y.all) <- as.character(1:n.all)
for (ijkl in 1:n.all) {
  add.i <- colnames(x.all)[which(x.all[ijkl,]!=0)][1]
  add.I <- sum(grepl("additive",colnames(x.all)[which(x.all[ijkl,]!=0)]))
  if (add.I>1) add.i <- "additive 22"                                          # this is "additive 22"
  ary.j <- colnames(x.all)[which(x.all[ijkl,]!=0)][1+add.I]
  ary.J <- sum(grepl("aryl_halide",colnames(x.all)[which(x.all[ijkl,]!=0)]))
  if (ary.J>1) ary.j <- "aryl_halide 15"                                       # this is "aryl_halide 15"
  bas.k <- colnames(x.all)[which(x.all[ijkl,]!=0)][1+add.I+ary.J]
  bas.K <- sum(grepl("base",colnames(x.all)[which(x.all[ijkl,]!=0)]))
  if (bas.K>1) bas.k <- "base 1"                                              # this is "base 1"
  lig.l <- colnames(x.all)[which(x.all[ijkl,]!=0)][1+add.I+ary.J+bas.K]
  lig.L <- sum(grepl("ligand",colnames(x.all)[which(x.all[ijkl,]!=0)]))
  if (lig.L>1) lig.l <- "ligand 3"                                            # this is "ligand 3"
  rownames(y.all)[ijkl] <- paste(add.i, ary.j, bas.k, lig.l, sep=":")
}

# Mixed terms with 2-levels combinations
xx <- rep(1,nrow(xcf))
bb <- cumsum(a-1)
for (j in 1:3) {
  for (i in (max(bb[j-1],0)+1):bb[j]) {
    xxp <- xcf[,i]*xcf[,-c(1:bb[j])]
    colnames(xxp) <- paste(colnames(xcf)[i],colnames(xcf[,-c(1:bb[j])]),sep=":")
    xx <- cbind(xx,xxp)
  }
}
xx <- cbind(xcf,xx[,-1])
xx.all <- xx

# Mixed terms with 3-levels combinations
xcf1 <- xcf
colnames(xcf1) <- c(rep("additive",21),rep("aryl_halide",14),rep("base",2),rep("ligand",3))
xx1 <- xx[,-c(1:40)]

xxx <- rep(1,nrow(xcf))
ind <- rep(TRUE,ncol(xx1))
for (j in 1:2) {
  ind <- as.logical((!grepl(colnames(xcf1)[bb[j]],colnames(xx1)))*(ind))
  for (i in (max(bb[j-1],0)+1):bb[j]) {
    xxxp <- xcf[,i]*xx1[,ind]
    colnames(xxxp) <- paste(colnames(xcf)[i],colnames(xx1)[ind],sep=":")
    xxx <- cbind(xxx,xxxp)
  }
}
xxx <- cbind(xx,xxx[,-1])
xxx.all <- xxx

# Mixed terms with 4-levels combinations
xxx1 <- xxx[,-c(1:515)]

xxxx <- rep(NA,nrow(xcf))
for (i in 1:21) {
  xxxxp <- xcf[,i]*xxx1[,1597:1680]
  colnames(xxxxp) <- paste(colnames(xcf)[i],colnames(xxx1)[1597:1680],sep=":")
  xxxx <- cbind(xxxx,xxxxp)
}
xxxx <- cbind(xxx,xxxx[,-1])
xxxx.all <- as.matrix(xxxx)
rownames(xxxx.all) <- rownames(y.all)



# Combinations
ind.1w <- 1:ncol(x.all)
ind.2w <- (ncol(x.all)+1):ncol(xx.all)
ind.3w <- (ncol(xx.all)+1):ncol(xxx.all)
ind.4w <- (ncol(xxx.all)+1):ncol(xxxx.all)

################################
## Yield Tables
################################

# Initialize table
df.yields.colnames <- c("addi C3 NMR shift", "aryl hali C1 NMR", "base N1 elec charge", "liga C10 NMR shift",
                        "y.obs", "mu.hat", "g(mu.hat)", "intercept",
                        "addi", "aryl", "base", "liga",
                        "addi:aryl", "addi:base", "addi:liga", "aryl:base", "aryl:liga", "base:liga",
                        "addi:aryl:base", "addi:aryl:liga", "addi:base:liga", "aryl:base:liga", "addi:aryl:base:liga")
df.yields <- data.frame(matrix(ncol=length(df.yields.colnames),nrow=0))
colnames(df.yields) <- df.yields.colnames

# Get coefficients of best model
beta.hat.shimNS.all <-  readRDS(fn_coefs)
eta.hat.shimNS.all <- as.numeric(xxxx.all%*%beta.hat.shimNS.all[-1] + beta.hat.shimNS.all[1])
mu.hat.shimNS.all <- kappa1(eta.hat.shimNS.all)
coef.shim <- rbind(c("intercept",beta.hat.shimNS.all[1]),data.frame(colnames(xxxx.all),beta.hat.shimNS.all[-1]))
coef.shim[,2] <- as.numeric(coef.shim[,2])
colnames(coef.shim) <- c("V1", "V2")

# Get top yields
n.yield.max <- length(y.all[,1])
IJKL <- order(y.all[,1],decreasing=TRUE)[1:n.yield.max]
y.IJKL <- as.data.frame(y.all[IJKL,1], row.names=rownames(y.all)[IJKL])  # top yields

for (ijkl in rownames(y.IJKL)) {
  
  # For the given yield get: relevant features, non-standard features, yield value, yield position in data, yield estimated mean
  tags.ijkl <- strsplit(ijkl, split=":")[[1]]
  spec.ijkl <- tags.ijkl[(grepl("additive 22",tags.ijkl))|
                           (grepl("aryl_halide 15",tags.ijkl))|
                           (grepl("base 1",tags.ijkl))|
                           (grepl("ligand 3",tags.ijkl))]        
  y.ijkl <- format(y.IJKL[ijkl,], digits=8, nsmall=8)
  ind.ijkl <- which(rownames(y.all)==ijkl)
  mu.ijkl <- mu.hat.shimNS.all[ind.ijkl]
  
  # Get coefficients with active features (we cross-multiply to get correct signs for positive and negative effects)
  active.ijkl <- which(xxxx.all[ind.ijkl,]!=0)
  coef.shim.ijkl <- rbind(coef.shim[1,], coef.shim[-1,][active.ijkl,])
  coef.shim.ijkl$V2 <- c(coef.shim.ijkl$V2[1], as.numeric(xxxx.all[ind.ijkl,active.ijkl])*coef.shim.ijkl$V2[-1])
  
  # Deal with one-level combinations
  for (spec in spec.ijkl) {
    spec.nonum <- strsplit(spec, split=" ")[[1]][1]
    cut.ijkl <- (!grepl(":",coef.shim.ijkl[,1]))&(grepl(spec.nonum,coef.shim.ijkl[,1]))
    coef.shim.ijkl[cut.ijkl,]
    new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
    new.val <- sum(coef.shim.ijkl[new.ind,2])
    new.tag <- spec
    coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
    coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
  }
  
  # Deal with two-level combinations
  for (tags in setdiff(tags.ijkl, spec.ijkl)) {
    for (spec in spec.ijkl) {
      spec.nonum <- strsplit(spec, split=" ")[[1]][1]
      cut.ijkl <- (grepl("1",lengths(regmatches(coef.shim.ijkl[,1], 
                                                gregexpr(":", coef.shim.ijkl[,1])))))&(grepl(tags,coef.shim.ijkl[,1]))&(grepl(spec.nonum,coef.shim.ijkl[,1]))
      new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
      new.val <- sum(coef.shim.ijkl[new.ind,2])
      new.tag <- paste(sort(c(tags, spec)),collapse=":")
      coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
      coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
    }
  }
  
  for (spec1 in spec.ijkl) {
    for (spec2 in setdiff(spec.ijkl,c(spec1))) {
      spec1.nonum <- strsplit(spec1, split=" ")[[1]][1]
      spec2.nonum <- strsplit(spec2, split=" ")[[1]][1]
      cut.ijkl <- (grepl("1",lengths(regmatches(coef.shim.ijkl[,1], 
                                                gregexpr(":", coef.shim.ijkl[,1])))))&(grepl(spec1.nonum,coef.shim.ijkl[,1]))&(grepl(spec2.nonum,coef.shim.ijkl[,1]))
      new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
      new.val <- sum(coef.shim.ijkl[new.ind,2])
      new.tag <- paste(sort(c(spec1, spec2)),collapse=":")
      coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
      coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
    }
  }
  
  # Deal with three-level combinations
  for (tags1 in setdiff(tags.ijkl, spec.ijkl)) {
    for (tags2 in setdiff(setdiff(tags.ijkl, spec.ijkl),c(tags1))) {
      for (spec in spec.ijkl) {
        spec.nonum <- strsplit(spec, split=" ")[[1]][1]
        cut.ijkl <- (grepl("2",lengths(regmatches(coef.shim.ijkl[,1], gregexpr(":", coef.shim.ijkl[,1])))))&
          (grepl(tags1,coef.shim.ijkl[,1]))&(grepl(tags2,coef.shim.ijkl[,1]))&(grepl(spec.nonum,coef.shim.ijkl[,1]))
        new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
        new.val <- sum(coef.shim.ijkl[new.ind,2])
        new.tag <- paste(sort(c(tags1, tags2, spec)),collapse=":")
        coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
        coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
      }
    }
  }
  
  for (tags in setdiff(tags.ijkl, spec.ijkl)) {
    for (spec1 in spec.ijkl) {
      for (spec2 in setdiff(spec.ijkl,c(spec1))) {
        spec1.nonum <- strsplit(spec1, split=" ")[[1]][1]
        spec2.nonum <- strsplit(spec2, split=" ")[[1]][1]
        cut.ijkl <- (grepl("2",lengths(regmatches(coef.shim.ijkl[,1], gregexpr(":", coef.shim.ijkl[,1])))))&
          (grepl(tags,coef.shim.ijkl[,1]))&(grepl(spec1.nonum,coef.shim.ijkl[,1]))&(grepl(spec2.nonum,coef.shim.ijkl[,1]))
        new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
        new.val <- sum(coef.shim.ijkl[new.ind,2])
        new.tag <- paste(sort(c(tags, spec1, spec2)),collapse=":")
        coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
        coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
      }
    }
  }
  
  for (spec1 in spec.ijkl) {
    for (spec2 in setdiff(spec.ijkl,c(spec1))) {
      for (spec3 in setdiff(spec.ijkl,c(spec1,spec2))) {
        spec1.nonum <- strsplit(spec1, split=" ")[[1]][1]
        spec2.nonum <- strsplit(spec2, split=" ")[[1]][1]
        spec3.nonum <- strsplit(spec3, split=" ")[[1]][1]
        cut.ijkl <- (grepl("2",lengths(regmatches(coef.shim.ijkl[,1], gregexpr(":", coef.shim.ijkl[,1])))))&
          (grepl(spec1.nonum,coef.shim.ijkl[,1]))&(grepl(spec2.nonum,coef.shim.ijkl[,1]))&(grepl(spec3.nonum,coef.shim.ijkl[,1]))
        new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
        new.val <- sum(coef.shim.ijkl[new.ind,2])
        new.tag <- paste(sort(c(spec1, spec2, spec3)),collapse=":")
        coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
        coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
      }
    }
  }
  
  # Deal with four-level combinations
  if (length(spec.ijkl)==1) {
    tags1 <- setdiff(tags.ijkl, spec.ijkl)[1]
    tags2 <- setdiff(tags.ijkl, spec.ijkl)[2]
    tags3 <- setdiff(tags.ijkl, spec.ijkl)[3]
    spec <- spec.ijkl[1]
    spec.nonum <- strsplit(spec, split=" ")[[1]][1]
    cut.ijkl <- (grepl("3",lengths(regmatches(coef.shim.ijkl[,1], gregexpr(":", coef.shim.ijkl[,1])))))&
      (grepl(tags1,coef.shim.ijkl[,1]))&(grepl(tags2,coef.shim.ijkl[,1]))&(grepl(tags3,coef.shim.ijkl[,1]))&(grepl(spec.nonum,coef.shim.ijkl[,1]))
    new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
    new.val <- sum(coef.shim.ijkl[new.ind,2])
    new.tag <- paste(sort(c(tags1, tags2, tags3, spec)),collapse=":")
    coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
    coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
  }
  
  if (length(spec.ijkl)==2) {
    tags1 <- setdiff(tags.ijkl, spec.ijkl)[1]
    tags2 <- setdiff(tags.ijkl, spec.ijkl)[2]
    spec1 <- spec.ijkl[1]
    spec2 <- spec.ijkl[2]
    spec1.nonum <- strsplit(spec1, split=" ")[[1]][1]
    spec2.nonum <- strsplit(spec2, split=" ")[[1]][1]
    cut.ijkl <- (grepl("3",lengths(regmatches(coef.shim.ijkl[,1], gregexpr(":", coef.shim.ijkl[,1])))))&
      (grepl(tags1,coef.shim.ijkl[,1]))&(grepl(tags2,coef.shim.ijkl[,1]))&(grepl(spec1.nonum,coef.shim.ijkl[,1]))&(grepl(spec2.nonum,coef.shim.ijkl[,1]))
    new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
    new.val <- sum(coef.shim.ijkl[new.ind,2])
    new.tag <- paste(sort(c(tags1, tags2, spec1, spec2)),collapse=":")
    coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
    coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
  }
  
  if (length(spec.ijkl)==3) {
    tags <- setdiff(tags.ijkl, spec.ijkl)[1]
    spec1 <- spec.ijkl[1]
    spec2 <- spec.ijkl[2]
    spec3 <- spec.ijkl[3]
    spec1.nonum <- strsplit(spec1, split=" ")[[1]][1]
    spec2.nonum <- strsplit(spec2, split=" ")[[1]][1]
    spec3.nonum <- strsplit(spec3, split=" ")[[1]][1]
    cut.ijkl <- (grepl("3",lengths(regmatches(coef.shim.ijkl[,1], gregexpr(":", coef.shim.ijkl[,1])))))&
      (grepl(tags,coef.shim.ijkl[,1]))&(grepl(spec1.nonum,coef.shim.ijkl[,1]))&(grepl(spec2.nonum,coef.shim.ijkl[,1]))&(grepl(spec3.nonum,coef.shim.ijkl[,1]))
    new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
    new.val <- sum(coef.shim.ijkl[new.ind,2])
    new.tag <- paste(sort(c(tags, spec1, spec2, spec3)),collapse=":")
    coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
    coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
  }
  
  if (length(spec.ijkl)==4) {
    spec1 <- spec.ijkl[1]
    spec2 <- spec.ijkl[2]
    spec3 <- spec.ijkl[3]
    spec4 <- spec.ijkl[4]
    spec1.nonum <- strsplit(spec1, split=" ")[[1]][1]
    spec2.nonum <- strsplit(spec2, split=" ")[[1]][1]
    spec3.nonum <- strsplit(spec3, split=" ")[[1]][1]
    spec4.nonum <- strsplit(spec4, split=" ")[[1]][1]
    cut.ijkl <- (grepl("3",lengths(regmatches(coef.shim.ijkl[,1], gregexpr(":", coef.shim.ijkl[,1])))))&
      (grepl(spec1.nonum,coef.shim.ijkl[,1]))&(grepl(spec2.nonum,coef.shim.ijkl[,1]))&(grepl(spec3.nonum,coef.shim.ijkl[,1]))&(grepl(spec4.nonum,coef.shim.ijkl[,1]))
    new.ind <- (1:length(coef.shim.ijkl[,1]))[cut.ijkl]
    new.val <- sum(coef.shim.ijkl[new.ind,2])
    new.tag <- paste(sort(c(spec1, spec2, spec3, spec4)),collapse=":")
    coef.shim.ijkl <- coef.shim.ijkl[-new.ind,]
    coef.shim.ijkl <- rbind(coef.shim.ijkl, list(new.tag,new.val))
  }
  
  # Recover original levels
  add.ijkl <- as.numeric(strsplit(sort(tags.ijkl)[1],split=" ")[[1]][2])
  ary.ijkl <- as.numeric(strsplit(sort(tags.ijkl)[2],split=" ")[[1]][2])
  bas.ijkl <- as.numeric(strsplit(sort(tags.ijkl)[3],split=" ")[[1]][2])
  lig.ijkl <- as.numeric(strsplit(sort(tags.ijkl)[4],split=" ")[[1]][2])
  
  add.true <- paste("additive ", apply(xs, 2, unique)$additive[which(apply(xs, 2, function(x) order(unique(x)))$additive==add.ijkl)], sep="")
  ary.true <- paste("aryl_halide ", apply(xs, 2, unique)$aryl_halide[which(apply(xs, 2, function(x) order(unique(x)))$aryl_halide==ary.ijkl)], sep="")
  bas.true <- paste("base ", apply(xs, 2, unique)$base[which(apply(xs, 2, function(x) order(unique(x)))$base==bas.ijkl)], sep="")
  lig.true <- paste("ligand ", apply(xs, 2, unique)$ligand[which(apply(xs, 2, function(x) order(unique(x)))$ligand==lig.ijkl)], sep="")
  
  coef.shim.ijkl <- coef.shim.ijkl[order(coef.shim.ijkl$V1),]
  coef.shim.ijkl$V2 <- format(as.numeric(coef.shim.ijkl$V2),nsmall=7,digits=7)
  coef.shim.ijkl$V1 <- sort(c("intercept", add.true, ary.true, bas.true, lig.true, paste(add.true, ary.true, bas.true, lig.true, sep=":"),
                              paste(add.true, ary.true, sep=":"), paste(add.true, bas.true, sep=":"), paste(add.true, lig.true, sep=":"),
                              paste(ary.true, bas.true, sep=":"), paste(ary.true, lig.true, sep=":"), paste(bas.true, lig.true, sep=":"), 
                              paste(add.true, ary.true, bas.true, sep=":"), paste(add.true, ary.true, lig.true, sep=":"),
                              paste(add.true, bas.true, lig.true, sep=":"), paste(ary.true, bas.true, lig.true, sep=":") ))
  coef.shim.ijkl <- rbind( c("g.link(mu.hat)", format(eta.hat.shimNS.all[ind.ijkl],nsmall=7,digits=7)), coef.shim.ijkl)
  coef.shim.ijkl <- rbind( c("mu.hat", format(mu.hat.shimNS.all[ind.ijkl],nsmall=7,digits=7)), coef.shim.ijkl)
  coef.shim.ijkl <- rbind( c("y.obs", format(y.all[ind.ijkl,],nsmall=7,digits=7)), coef.shim.ijkl)
  
  # Update global table
  df.yields <- rbind(df.yields, 
                     c(strsplit(add.true,split=" ")[[1]][2], strsplit(ary.true,split=" ")[[1]][2], strsplit(bas.true,split=" ")[[1]][2], strsplit(lig.true,split=" ")[[1]][2],
                       coef.shim.ijkl[1,2], coef.shim.ijkl[2,2], coef.shim.ijkl[3,2], coef.shim.ijkl[18,2],
                       coef.shim.ijkl[4,2], coef.shim.ijkl[12,2], coef.shim.ijkl[16,2], coef.shim.ijkl[19,2],
                       coef.shim.ijkl[5,2], coef.shim.ijkl[9,2], coef.shim.ijkl[11,2], 
                       coef.shim.ijkl[13,2], coef.shim.ijkl[15,2], coef.shim.ijkl[17,2],
                       coef.shim.ijkl[6,2], coef.shim.ijkl[8,2], coef.shim.ijkl[10,2], coef.shim.ijkl[14,2], coef.shim.ijkl[7,2]) )
  colnames(df.yields) <- df.yields.colnames
  
}
gc()

# Save global table
write_xlsx(df.yields, path=paste("../Results/Table_of_coefficients.xlsx"))
print("Table is stored in folder Results")
}

