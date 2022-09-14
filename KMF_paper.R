######################################################################
## 1. COMPARISON OF LINEAR AND CONTINUOUS BERNOULLI MODELS
######################################################################
source("tools.r")
library(ggplot2)

#######################
## Data processing
#######################

# Read data
x <- data.matrix(read.table("Data/Xtrain.txt",header=TRUE))
x.test <- data.frame(read.table("Data/Xtest.txt",header=TRUE))
y <- as.vector(read.table("Data/YTrain.txt",header=FALSE)[,1])
y.test <- as.vector(read.table("Data/YTest.txt",header=FALSE)[,1])

# Remove linear dependent columns
ind <- qr(x)$pivot[seq_len(qr(x)$rank)]
x <- x[, ind]
x.test <- x.test[,ind]
ind1 <- which(is.na(as.vector(lm(y~x)$coef)))-1
x <- x[,-ind1]
x.test <- x.test[,-ind1] 

#################################################
## Fit continuous Bernoulli and plot histogram
#################################################

y.all <- c(y,y.test)
n.all <- length(y.all)
tt <- 1:n.all/n.all

# Find MLE for p

ld<-function(p)
{
  (n.all*p/(2*p-1)-sum(y.all)+n.all/log((1-p)/p))/(p*(p-1))
}
p.hat <- uniroot(ld,c(0.001,0.5))$root

# Plot the histogram

ll <- log((1-p.hat)/p.hat)*p.hat^tt*(1-p.hat)^(1-tt)/(1-2*p.hat)
df <- data.frame(y.all,tt,ll)
colnames(df) <- c("Yield","Probability","Density")

pdf("yield_hist.pdf",width=10,height=10)
ggplot(df, aes(x=Yield*100)) + 
  geom_histogram(aes(y=..density..),color="black",fill="cornflowerblue",breaks=seq(0,100,length=11))+ 
  theme(text = element_text(size = 40)) +labs(title = "Yield density") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_line(aes(x=tt*100,y =ll/100), color = "black",size=2)+ylab("Density")+xlab("Yield")
dev.off()


########################
## Linear regression
########################

fit.lm <- lm(y~x)
pred.lm <- as.matrix(x.test)%*%fit.lm$coef[-1]+fit.lm$coef[1]
rsq.lm <- as.numeric(cor(pred.lm,y.test)^2)
rmse.lm <- sqrt(mean((pred.lm*100-y.test*100)^2))

df.lm=data.frame(y.test,pred.lm)
colnames(df.lm)=c("ObservedYield","PredictedYield")

pdf("linear_model.pdf",width=10,height=10)
ggplot(data = df.lm, aes(x = PredictedYield*100, y = ObservedYield*100)) + 
  theme(text = element_text(size = 40)) +labs(title = "Linear model") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+
  xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-28),y=c(95,89),label=c(paste("R^2:",round(rsq.lm,2)),
                                                  paste("RMSE:",round(rmse.lm,2))), size=10)
dev.off()

#################################################
### Regression with continuous Bernoulli response
#################################################

# Set starting value
beta.start <- lm(y~x)$coef

# Prediction
beta.hat <- fisher.scoring(y,x,beta.start)$beta.hat
pred.cb <- kappa1(as.vector(as.matrix(x.test)%*%beta.hat[-1]+beta.hat[1]))
rsq.cb <- as.numeric(cor(pred.cb,y.test)^2)
rmse.cb <- sqrt(mean((pred.cb*100-y.test*100)^2))

df.cb=data.frame(y.test,pred.cb)

pdf("cb_model.pdf",width=10,height=10)
ggplot(data = df.cb, aes(x = pred.cb*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = "Generalised linear model") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-28),y=c(95,89),label=c(paste("R^2:",round(rsq.cb,2)),
                                                  paste("RMSE:",round(rmse.cb,2))), size=10)
dev.off()



######################################################################
## 2. COMPARISON BETWEEN DUMMY CODING AND CHEMICAL DESCRIPTORS
######################################################################
source("tools.r")
library(ggplot2)

#######################
## Data processing
#######################

# Read data
x <- data.matrix(read.table("Data/XtrainNScaled.txt",header=TRUE))
x.test <- data.frame(read.table("Data/XtestNScaled.txt",header=TRUE))
y <- as.vector(read.table("Data/YTrain.txt",header=FALSE)[,1])
y.test <- as.vector(read.table("Data/YTest.txt",header=FALSE)[,1])

# Add three artificial descriptors for additives
xc <- rbind(x,x.test)
set.seed(2134)
c1 <- factor(xc[,1])
c2 <- factor(xc[,3])
c3 <- factor(xc[,4])
levels(c1) <- runif(22)
levels(c2) <- rnorm(22)
levels(c3) <- runif(22)
xc <- cbind(cbind(as.numeric(c1),as.numeric(c2),as.numeric(c3)),xc)
colnames(xc)[1:3] <- c("add_new1","add_new2","add_new3")

xs <- xc[,c(1,23,50,60)]
colnames(xs) <- c("additive","aryl_halide","base","ligand")

a <- rep(NA,ncol(xs))
for (i in 1:ncol(xs)) a[i] <- length(unique(xs[,i]))
xcf <- matrix(NA,nrow(xs),sum(a))
b <- cumsum(a)
colnames(xcf) <- rep(colnames(xs),times=a)
colnum <- order(unique(xs[,1]))
for (i in 2:ncol(xs)) colnum <- c(colnum,order(unique(xs[,i])))
colnames(xcf) <- paste(colnames(xcf),colnum)

for (i in 1:nrow(xs))
{
  for (j in 1:length(a))
  {
    res <- rep(0, a[j])
    where <- match( xs[i,j], unique(xs[,j]) )
    res[ where ] <- 1 
    xcf[i,(max(b[j-1],0)+1):b[j]] <- res
  }
}

for (i in 1:length(b))
{
  ind <- match(xcf[,b[i]],1)==1
  xcf[ind,(max(b[i-1],0)+1):b[i]] <- -1
}

xcf <- xcf[,-b]

xf <- xcf[1:nrow(x),]
xf.test <- xcf[-c(1:nrow(x)),]

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

xx.train <- xx[1:nrow(x),]
xx.test <- xx[-c(1:nrow(x)),]

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

xxx.train <- xxx[1:nrow(x),]
xxx.test <- xxx[-c(1:nrow(x)),]

# Mixed terms with 4-levels combinations
xxx1 <- xxx[,-c(1:515)]

xxxx <- rep(NA,nrow(xcf))
for (i in 1:21) {
  xxxxp <- xcf[,i]*xxx1[,1597:1680]
  colnames(xxxxp) <- paste(colnames(xcf)[i],colnames(xxx1)[1597:1680],sep=":")
  xxxx <- cbind(xxxx,xxxxp)
}
xxxx <- cbind(xxx,xxxx[,-1])

xxxx.train <- xxxx[1:nrow(x),]
xxxx.test <- xxxx[-c(1:nrow(x)),]

#######################
## Linear models
#######################

# Fit linear models with dummy coding and descriptors matrix
fit1=lm(y.all~as.matrix(scale(xc)))
fit2=lm(y.all~xcf)
plot(fitted(fit1),fitted(fit2),type="l")

# Compare to the model with 19 descriptors
fit3=lm(y.all~scale(as.matrix(rbind(x,x.test))))
plot(fitted(fit1),fitted(fit3),type="l")

#######################
## Random forest
#######################

library(randomForest)

set.seed(1234)
# RF with 120 correlated predictors
fit.rf=randomForest(y~.,data=x)

pred.rf=predict(fit.rf,x.test)
plot(pred.rf,y.test)
abline(0,1,lwd=3,col=2)
cor(pred.rf,y.test)^2
sqrt(mean((pred.rf*100-y.test*100)^2))

#RF with 39 predictors
ind=qr(x)$pivot[seq_len(qr(x)$rank)]
x=x[, ind]
x.test=x.test[,ind]

fit.rf=randomForest(y~.,data=x)

pred.rf=predict(fit.rf,x.test)
plot(pred.rf,y.test)
abline(0,1,lwd=3,col=2)
cor(pred.rf,y.test)^2
sqrt(mean((pred.rf*100-y.test*100)^2))

#RF with dummy coding 
fit.rf=randomForest(y~.,data=data.frame(xf))

pred.rf=predict(fit.rf,data.frame(xf.test))
plot(pred.rf,y.test)
abline(0,1,lwd=3,col=2)
cor(pred.rf,y.test)^2
sqrt(mean((pred.rf*100-y.test*100)^2))

############################
## PLSGLM on 2-level model
############################

nc <- 30
scaling <- FALSE
beta.hat.pls2 <- plsglm.cb.simple(xx.train, y, nc, scaling=scaling)
gc()

eta.hat.pls2 <- as.vector(xx.train%*%beta.hat.pls2[-1]+beta.hat.pls2[1])
fit.pls2 <- kappa1(eta.hat.pls2)
ll.pls2 <- sum(y*eta.hat.pls2-kappa(eta.hat.pls2))

pred.pls2 <- kappa1(as.vector(xx.test%*%beta.hat.pls2[-1]+beta.hat.pls2[1]))
rsq.pls2 <- cor(pred.pls2,y.test)^2
rmse.pls2 <- sqrt(mean((pred.pls2*100-y.test*100)^2))

df.pls2 <- data.frame(y.test,pred.pls2)

filename.pls2 <- "plsglm_model2.pdf"
if (!scaling) filename.pls2 <- "plsglm_model2_NS.pdf"
pdf(filename.pls2,width=10,height=10)
ggplot(data = df.pls2, aes(x = pred.pls2*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = paste("PLSGLM-",nc," with 2 levels",sep="")) +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-31),y=c(95,89),label=c(paste("R^2:",round(rsq.pls2,2)),
                                                  paste("RMSE:",round(rmse.pls2,2))), size=10)
dev.off()

coef.pls2 <- rbind(c("intercept",beta.hat.pls2[1]),cbind(colnames(xx),beta.hat.pls2[-1]))[order(abs(beta.hat.pls2)),]
coef.pls2[496:516,]

############################
## GLM on 2-level model
############################

beta.hat2 <- fisher.scoring(y,xx.train,beta.hat.pls2)$beta.hat #has convergence problems due to many beta=0 weights become close to 0
gc()

eta.hat2 <- as.vector(xx.train%*%beta.hat2[-1]+beta.hat2[1])
fit.hat2 <- kappa1(eta.hat2)
ll.hat2 <- sum(y*eta.hat2-kappa(eta.hat2))

pred.hat2 <- kappa1(as.vector(xx.test%*%beta.hat2[-1]+beta.hat2[1]))
rsq.hat2 <- cor(pred.hat2,y.test)^2
rmse.hat2 <- sqrt(mean((pred.hat2*100-y.test*100)^2))

df.hat2 <- data.frame(y.test,pred.hat2)

pdf("cb_model2.pdf",width=10,height=10)
ggplot(data = df.hat2, aes(x = pred.hat2*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = "GLM with 2 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-28),y=c(95,89),label=c(paste("R^2:",round(rsq.hat2,2)),
                                                  paste("RMSE:",round(rmse.hat2,2))), size=10)
dev.off()

rbind(c("intercept",beta.hat2[1]),cbind(colnames(xx),beta.hat2[-1]))[order(abs(beta.hat2)),]

############################
## PLSGLM on 3-level model
############################

memory.limit(32000)

nc <- 30
scaling <- FALSE
beta.hat.pls3 <- plsglm.cb.simple(xxx.train, y, nc, scaling=scaling)
gc()

eta.hat.pls3 <- as.vector(xxx.train%*%beta.hat.pls3[-1]+beta.hat.pls3[1])
fit.pls3 <- kappa1(eta.hat.pls3)
ll.pls3 <- sum(y*eta.hat.pls3-kappa(eta.hat.pls3))

pred.pls3 <- kappa1(as.vector(xxx.test%*%beta.hat.pls3[-1]+beta.hat.pls3[1]))
rsq.pls3 <- cor(pred.pls3,y.test)^2
rmse.pls3 <- sqrt(mean((pred.pls3*100-y.test*100)^2))

df.pls3 <- data.frame(y.test,pred.pls3)

filename.pls3 <- "plsglm_model3.pdf"
if (!scaling) filename.pls3 <- "plsglm_model3_NS.pdf"
pdf(filename.pls3,width=10,height=10)
ggplot(data = df.pls3, aes(x = pred.pls3*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = paste("PLSGLM-",nc," with 3 levels",sep="")) +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-31),y=c(95,89),label=c(paste("R^2:",round(rsq.pls3,2)),
                                                  paste("RMSE:",round(rmse.pls3,2))), size=10)
dev.off()

coef.pls3 <- rbind(c("intercept",beta.hat.pls3[1]),data.frame(colnames(xxx),beta.hat.pls3[-1]))[order(abs(beta.hat.pls3)),]
coef.pls3[2186:2196,]

# Get all coefficients from the constraints
coef.one <- coef.pls3[2:41,]
for (i in 1:4) {
  ind <- grepl(colnames(xs)[i],coef.one[,1],fixed=TRUE)
  coef.one <- rbind(coef.one,c(paste(colnames(xs[i]),"last"),-sum(as.numeric(coef.one[ind,2]))))
}

coef.two <- coef.pls3[c(42:516),]
for (k in 1:4){
  for (j in 1:3)
    for (i in (max(bb[k-1],0)+1):bb[k])
    { 
      ind <- grepl(paste(colnames(xf)[i],(colnames(xs)[-k])[j],sep=":"),coef.two[,1],fixed=TRUE)
      if(sum(ind)!=0) coef.two <- rbind(coef.two,c(paste(paste(colnames(xf)[i],(colnames(xs)[-k])[j],sep=":"),"last"),-sum(as.numeric(coef.two[ind,2]))))
    }
}

coef.two.rev <- coef.pls3[c(42:516),]
ct.s <- strsplit(coef.two.rev[,1], ":")
for (i in 1:length(ct.s)) {coef.two.rev[i,1]=paste(ct.s[[i]][2],ct.s[[i]][1],sep=":")}

for (k in 1:4){
  for (j in 1:3)
    for (i in (max(bb[k-1],0)+1):bb[k])
    { 
      ind <- grepl(paste(colnames(xf)[i],(colnames(xs)[-k])[j],sep=":"),coef.two.rev[,1],fixed=TRUE)
      if(sum(ind)!=0) coef.two.rev <- rbind(coef.two.rev,c(paste(paste(colnames(xf)[i],(colnames(xs)[-k])[j],sep=":"),"last"),-sum(as.numeric(coef.two.rev[ind,2]))))
    }
}
coef.two <- rbind(coef.two,coef.two.rev[476:502,])

############################
## GLM on 3-level model
############################

beta.hat3 <- fisher.scoring(y,xxx.train,beta.hat.pls3)$beta.hat #same convergence problems
gc()

eta.hat3 <- as.vector(xxx.train%*%beta.hat3[-1]+beta.hat3[1])
fit.hat3 <- kappa1(eta.hat3)
ll.hat3 <- sum(y*eta.hat3-kappa(eta.hat3))

pred.hat3 <- kappa1(as.vector(xxx.test%*%beta.hat3[-1]+beta.hat3[1]))
rsq.hat3 <- cor(pred.hat3,y.test)^2
rmse.hat3 <- sqrt(mean((pred.hat3*100-y.test*100)^2))

df.hat3 <- data.frame(y.test,pred.hat3)

pdf("cb_model3.pdf",width=10,height=10)
ggplot(data = df.hat3, aes(x = pred.hat3*100, y = y.test*100)) + theme(text = element_text(size = 40)) +labs(title = "GLM with 3 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-28),y=c(95,89),label=c(paste("R^2:",round(rsq.hat3,2)),
                                                  paste("RMSE:",round(rmse.hat3,2))), size=10)
dev.off()

rbind(c("intercept",beta.hat3[1]),data.frame(colnames(xxx.train),beta.hat3[-1]))[order(abs(beta.hat3)),]

coef.xxx <- rbind(c("intercept",beta.hat3[1]),coef.one,coef.two,coef.pls3[-c(1:516),])
(coef.xxx[order(abs(as.numeric(coef.xxx[,2]))),])[2300:2320,]

############################
## PLSGLM on 4-level model
############################

nc <- 30
scaling <- FALSE
beta.hat.pls4 <- plsglm.cb.simple(xxxx.train, y, nc, scaling=scaling)
gc()

eta.hat.pls4 <- as.vector(xxxx.train%*%beta.hat.pls4[-1]+beta.hat.pls4[1])
fit.pls4 <- kappa1(eta.hat.pls4)
ll.pls4 <- sum(y*eta.hat.pls4-kappa(eta.hat.pls4))

pred.pls4 <- kappa1(as.vector(xxxx.test%*%beta.hat.pls4[-1]+beta.hat.pls4[1]))
rsq.pls4 <- cor(pred.pls4,y.test)^2
rmse.pls4 <- sqrt(mean((pred.pls4*100-y.test*100)^2))

df.pls4 <- data.frame(y.test,pred.pls4)

filename.pls4 <- "plsglm_model4.pdf"
if (!scaling) filename.pls4 <- "plsglm_model4_NS.pdf"
pdf(filename.pls4,width=10,height=10)
ggplot(data = df.pls4, aes(x = pred.pls4*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = paste("PLSGLM-",nc," with 4 levels",sep="")) +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-29),y=c(95,89),label=c(paste("R^2:",round(rsq.pls4,2)),
                                                  paste("RMSE:",round(rmse.pls4,2))), size=10)
dev.off()

coef.pls4 <- rbind(c("intercept",beta.hat.pls4[1]),cbind(colnames(xxxx),beta.hat.pls4[-1]))[order(abs(beta.hat.pls4)),]
coef.pls4[3940:3960,]



######################################################################
## 3. MODEL SELECTION FOR PLSGLM WITH DIFFERENT LEVELS
######################################################################
source("tools.r")
library(ggplot2)

memory.limit(32000)

#######################
## Data processing
#######################

# Read data
x <- data.matrix(read.table("Data/XtrainNScaled.txt",header=TRUE))
x.test <- data.frame(read.table("Data/XtestNScaled.txt",header=TRUE))
y <- as.vector(read.table("Data/YTrain.txt",header=FALSE)[,1])
y.test <- as.vector(read.table("Data/YTest.txt",header=FALSE)[,1])

# Add three artificial descriptors for additives
xc <- rbind(x,x.test)
set.seed(2134)
c1 <- factor(xc[,1])
c2 <- factor(xc[,3])
c3 <- factor(xc[,4])
levels(c1) <- runif(22)
levels(c2) <- rnorm(22)
levels(c3) <- runif(22)
xc <- cbind(cbind(as.numeric(c1),as.numeric(c2),as.numeric(c3)),xc)
colnames(xc)[1:3] <- c("add_new1","add_new2","add_new3")

xs <- xc[,c(1,23,50,60)]
colnames(xs) <- c("additive","aryl_halide","base","ligand")

a <- rep(NA,ncol(xs))
for (i in 1:ncol(xs)) a[i] <- length(unique(xs[,i]))
xcf <- matrix(NA,nrow(xs),sum(a))
b <- cumsum(a)
colnames(xcf) <- rep(colnames(xs),times=a)
colnum <- order(unique(xs[,1]))
for (i in 2:ncol(xs)) colnum <- c(colnum,order(unique(xs[,i])))
colnames(xcf) <- paste(colnames(xcf),colnum)

for (i in 1:nrow(xs))
{
  for (j in 1:length(a))
  {
    res <- rep(0, a[j])
    where <- match( xs[i,j], unique(xs[,j]) )
    res[ where ] <- 1 
    xcf[i,(max(b[j-1],0)+1):b[j]] <- res
  }
}

for (i in 1:length(b))
{
  ind <- match(xcf[,b[i]],1)==1
  xcf[ind,(max(b[i-1],0)+1):b[i]] <- -1
}

xcf <- xcf[,-b]

xf <- xcf[1:nrow(x),]
xf.test <- xcf[-c(1:nrow(x)),]

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

xx.train <- xx[1:nrow(x),]
xx.test <- xx[-c(1:nrow(x)),]

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

xxx.train <- xxx[1:nrow(x),]
xxx.test <- xxx[-c(1:nrow(x)),]

# Mixed terms with 4-levels combinations
xxx1 <- xxx[,-c(1:515)]

xxxx <- rep(NA,nrow(xcf))
for (i in 1:21) {
  xxxxp <- xcf[,i]*xxx1[,1597:1680]
  colnames(xxxxp) <- paste(colnames(xcf)[i],colnames(xxx1)[1597:1680],sep=":")
  xxxx <- cbind(xxxx,xxxxp)
}
xxxx <- cbind(xxx,xxxx[,-1])

xxxx.train <- xxxx[1:nrow(x),]
xxxx.test <- xxxx[-c(1:nrow(x)),]

#######################
## PLSGLM with CV
#######################

nc <- (4:14)*3
cor2 <- rmse2 <- cor3 <- rmse3 <- cor4 <- rmse4 <- ll2 <- ll3 <- ll4 <- rep(NA,max(nc))
cor2NS <- rmse2NS <- cor3NS <- rmse3NS <- cor4NS <- rmse4NS <- ll2NS <- ll3NS <- ll4NS <- rep(NA,max(nc))

for (i in nc) {
  
  print(paste("Checking",i,"comps..."))
  
  # Scaling
  beta.hat.pls2 <- plsglm.cb.simple(xx.train, y, i, scaling=TRUE)
  gc()
  pred.cb.pls2 <- kappa1(as.vector(as.matrix(xx.test)%*%beta.hat.pls2[-1]+beta.hat.pls2[1]))
  cor2[i] <- cor(pred.cb.pls2,y.test)
  rmse2[i] <- sqrt(mean((pred.cb.pls2*100-y.test*100)^2))
  eta.hat.pls2 <- as.vector(xx.train%*%beta.hat.pls2[-1]+beta.hat.pls2[1])
  ll2[i] <- sum(y*eta.hat.pls2-kappa(eta.hat.pls2))
  
  beta.hat.pls3 <- plsglm.cb.simple(xxx.train, y, i, scaling=TRUE)
  gc()
  pred.cb.pls3 <- kappa1(as.vector(as.matrix(xxx.test)%*%beta.hat.pls3[-1]+beta.hat.pls3[1]))
  cor3[i] <- cor(pred.cb.pls3,y.test)
  rmse3[i] <- sqrt(mean((pred.cb.pls3*100-y.test*100)^2))
  eta.hat.pls3 <- as.vector(xxx.train%*%beta.hat.pls3[-1]+beta.hat.pls3[1])
  ll3[i] <- sum(y*eta.hat.pls3-kappa(eta.hat.pls3))
  
  beta.hat.pls4 <- plsglm.cb.simple(xxxx.train, y, i, scaling=TRUE)
  gc()
  pred.cb.pls4 <- kappa1(as.vector(as.matrix(xxxx.test)%*%beta.hat.pls4[-1]+beta.hat.pls4[1]))
  cor4[i] <- cor(pred.cb.pls4,y.test)
  rmse4[i] <- sqrt(mean((pred.cb.pls4*100-y.test*100)^2))
  eta.hat.pls4 <- as.vector(xxxx.train%*%beta.hat.pls4[-1]+beta.hat.pls4[1])
  ll4[i] <- sum(y*eta.hat.pls4-kappa(eta.hat.pls4))
  
  # Non-scaling
  beta.hat.pls2NS <- plsglm.cb.simple(xx.train, y, i, scaling=FALSE)
  gc()
  pred.cb.pls2NS <- kappa1(as.vector(as.matrix(xx.test)%*%beta.hat.pls2NS[-1]+beta.hat.pls2NS[1]))
  cor2NS[i] <- cor(pred.cb.pls2NS,y.test)
  rmse2NS[i] <- sqrt(mean((pred.cb.pls2NS*100-y.test*100)^2))
  eta.hat.pls2NS <- as.vector(xx.train%*%beta.hat.pls2NS[-1]+beta.hat.pls2NS[1])
  ll2NS[i] <- sum(y*eta.hat.pls2NS-kappa(eta.hat.pls2NS))
  
  beta.hat.pls3NS <- plsglm.cb.simple(xxx.train, y, i, scaling=FALSE)
  gc()
  pred.cb.pls3NS <- kappa1(as.vector(as.matrix(xxx.test)%*%beta.hat.pls3NS[-1]+beta.hat.pls3NS[1]))
  cor3NS[i] <- cor(pred.cb.pls3NS,y.test)
  rmse3NS[i] <- sqrt(mean((pred.cb.pls3NS*100-y.test*100)^2))
  eta.hat.pls3NS <- as.vector(xxx.train%*%beta.hat.pls3NS[-1]+beta.hat.pls3NS[1])
  ll3NS[i] <- sum(y*eta.hat.pls3NS-kappa(eta.hat.pls3NS))
  
  beta.hat.pls4NS <- plsglm.cb.simple(xxxx.train, y, i, scaling=FALSE)
  gc()
  pred.cb.pls4NS <- kappa1(as.vector(as.matrix(xxxx.test)%*%beta.hat.pls4NS[-1]+beta.hat.pls4NS[1]))
  cor4NS[i] <- cor(pred.cb.pls4NS,y.test)
  rmse4NS[i] <- sqrt(mean((pred.cb.pls4NS*100-y.test*100)^2))
  eta.hat.pls4NS <- as.vector(xxxx.train%*%beta.hat.pls4NS[-1]+beta.hat.pls4NS[1])
  ll4NS[i] <- sum(y*eta.hat.pls4NS-kappa(eta.hat.pls4NS))

}
gc()

# Visual results
plot(c(1,max(nc)),c(0,1),col="white",main="R^2 test",xaxt="n")
axis(1, at=nc)
points(nc,cor2[nc]^2,col=2,type="b")
points(nc,cor3[nc]^2,col=3,type="b")
points(nc,cor4[nc]^2,col=4,type="b")
points(nc,cor2NS[nc]^2,col=2,type="b",pch=2)
points(nc,cor3NS[nc]^2,col=3,type="b",pch=2)
points(nc,cor4NS[nc]^2,col=4,type="b",pch=2)
legend("topleft", col=c(2,2,3,3,4,4), 
       lwd=2, cex=1,
       lty=1,
       c("2-levels","2-levelsNS","3-levels","3-levelsNS","4-levels","4-levelsNS"),
       pch=c(1,2,1,2,1,2))

plot(c(1,max(nc)),c(0,50),col="white",main="RSME test",xaxt="n")
axis(1, at=nc)
points(nc,rmse2[nc],col=2,type="b")
points(nc,rmse3[nc],col=3,type="b")
points(nc,rmse4[nc],col=4,type="b")
points(nc,rmse2NS[nc],col=2,type="b",pch=2)
points(nc,rmse3NS[nc],col=3,type="b",pch=2)
points(nc,rmse4NS[nc],col=4,type="b",pch=2)
legend("topleft", col=c(2,2,3,3,4,4), 
       lwd=2, cex=1,
       lty=1,
       c("2-levels","2-levelsNS","3-levels","3-levelsNS","4-levels","4-levelsNS"),
       pch=c(1,2,1,2,1,2))

plot(c(1,max(nc)),c(0,3000),col="white",main="LogL train",xaxt="n")
axis(1, at=nc)
points(nc,ll2[nc],col=2,type="b")
points(nc,ll3[nc],col=3,type="b")
points(nc,ll4[nc],col=4,type="b")
points(nc,ll2NS[nc],col=2,type="b",pch=2)
points(nc,ll3NS[nc],col=3,type="b",pch=2)
points(nc,ll4NS[nc],col=4,type="b",pch=2)
legend("topleft", col=c(2,2,3,3,4,4), 
       lwd=2, cex=1,
       lty=1,
       c("2-levels","2-levelsNS","3-levels","3-levelsNS","4-levels","4-levelsNS"),
       pch=c(1,2,1,2,1,2))

# Save results
write(cor2, file="cv_cor2.txt", sep="\n")
write(cor3, file="cv_cor3.txt", sep="\n")
write(cor4, file="cv_cor4.txt", sep="\n")
write(cor2NS, file="cv_cor2NS.txt", sep="\n")
write(cor3NS, file="cv_cor3NS.txt", sep="\n")
write(cor4NS, file="cv_cor4NS.txt", sep="\n")

write(rmse2, file="cv_rmse2.txt", sep="\n")
write(rmse3, file="cv_rmse3.txt", sep="\n")
write(rmse4, file="cv_rmse4.txt", sep="\n")
write(rmse2NS, file="cv_rmse2NS.txt", sep="\n")
write(rmse3NS, file="cv_rmse3NS.txt", sep="\n")
write(rmse4NS, file="cv_rmse4NS.txt", sep="\n")

write(ll2, file="cv_ll2.txt", sep="\n")
write(ll3, file="cv_ll3.txt", sep="\n")
write(ll4, file="cv_ll4.txt", sep="\n")
write(ll2NS, file="cv_ll2NS.txt", sep="\n")
write(ll3NS, file="cv_ll3NS.txt", sep="\n")
write(ll4NS, file="cv_ll4NS.txt", sep="\n")

# Select best models
nc2 <- which(cor2==max(cor2[nc]))
nc2NS <- which(cor2NS==max(cor2NS[nc]))
nc3 <- which(cor3==max(cor3[nc]))
nc3NS <- which(cor3NS==max(cor3NS[nc]))
nc4 <- which(cor4==max(cor4[nc]))
nc4NS <- which(cor4NS==max(cor4NS[nc]))

plot(c(1,6),c(min(nc),max(nc)),col="white",xaxt="n",yaxt="n",main="Best models R^2")
axis(2,at=nc)
points(c(nc2,nc2NS,nc3,nc3NS,nc4,nc4NS), 
     col=c(2,2,3,3,4,4),
     pch=c(16,17,16,17,16,17),
     cex=2)
text(c(nc2,nc2NS,nc3,nc3NS,nc4,nc4NS)-1,
     labels=c(round(cor2[nc2]^2,4),round(cor2NS[nc2NS]^2,4),
              round(cor3[nc3]^2,4),round(cor3NS[nc3NS]^2,4),
              round(cor4[nc4]^2,4),round(cor4NS[nc4NS]^2,4)),
     col=c(2,2,3,3,4,4))
legend("topright", 
       col=c(2,2,3,3,4,4), 
       lwd=1, cex=1, lty=NA,
       c("2-lev","2-levNS","3-lev","3-levNS","4-lev","4-levNS"),
       pch=c(16,17,16,17,16,17))

# Best 2-levels models
beta.hat.pls2 <- plsglm.cb.simple(xx.train, y, nc2, scaling=TRUE)
gc()

eta.hat.pls2 <- as.vector(xx.train%*%beta.hat.pls2[-1]+beta.hat.pls2[1])
ll.pls2 <- sum(y*eta.hat.pls2-kappa(eta.hat.pls2))
pred.pls2 <- kappa1(as.vector(xx.test%*%beta.hat.pls2[-1]+beta.hat.pls2[1]))
rsq.pls2 <- cor(pred.pls2,y.test)^2
rmse.pls2 <- sqrt(mean((pred.pls2*100-y.test*100)^2))

df.pls2 <- data.frame(y.test,pred.pls2)

pdf("cv_plsglm_model2.pdf",width=10,height=10)
ggplot(data = df.pls2, aes(x = pred.pls2*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = paste("PLSGLM-",nc2," with 2 levels",sep="")) +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-31,-27),y=c(95,89,83),label=c(paste("R^2:",round(rsq.pls2,2)),
                                                  paste("RMSE:",round(rmse.pls2,2)),
                                                  paste("LogL:",round(ll.pls2,2))), size=10)
dev.off()

# Best 3-levels models
beta.hat.pls3 <- plsglm.cb.simple(xxx.train, y, nc3, scaling=TRUE)
gc()

eta.hat.pls3 <- as.vector(xxx.train%*%beta.hat.pls3[-1]+beta.hat.pls3[1])
ll.pls3 <- sum(y*eta.hat.pls3-kappa(eta.hat.pls3))
pred.pls3 <- kappa1(as.vector(xxx.test%*%beta.hat.pls3[-1]+beta.hat.pls3[1]))
rsq.pls3 <- cor(pred.pls3,y.test)^2
rmse.pls3 <- sqrt(mean((pred.pls3*100-y.test*100)^2))

df.pls3 <- data.frame(y.test,pred.pls3)

pdf("cv_plsglm_model3.pdf",width=10,height=10)
ggplot(data = df.pls3, aes(x = pred.pls3*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = paste("PLSGLM-",nc3," with 3 levels",sep="")) +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-29,-27),y=c(95,89,83),label=c(paste("R^2:",round(rsq.pls3,2)),
                                                         paste("RMSE:",round(rmse.pls3,2)),
                                                         paste("LogL:",round(ll.pls3,2))), size=10)
dev.off()

# Best 4-levels models
beta.hat.pls4 <- plsglm.cb.simple(xxxx.train, y, nc4, scaling=TRUE)
gc()

eta.hat.pls4 <- as.vector(xxxx.train%*%beta.hat.pls4[-1]+beta.hat.pls4[1])
ll.pls4 <- sum(y*eta.hat.pls4-kappa(eta.hat.pls4))
pred.pls4 <- kappa1(as.vector(xxxx.test%*%beta.hat.pls4[-1]+beta.hat.pls4[1]))
rsq.pls4 <- cor(pred.pls4,y.test)^2
rmse.pls4 <- sqrt(mean((pred.pls4*100-y.test*100)^2))

df.pls4 <- data.frame(y.test,pred.pls4)

pdf("cv_plsglm_model4.pdf",width=10,height=10)
ggplot(data = df.pls4, aes(x = pred.pls4*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = paste("PLSGLM-",nc4," with 4 levels",sep="")) +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-29,-27),y=c(95,89,83),label=c(paste("R^2:",round(rsq.pls4,2)),
                                                         paste("RMSE:",round(rmse.pls4,2)),
                                                         paste("LogL:",round(ll.pls4,2))), size=10)
dev.off()




