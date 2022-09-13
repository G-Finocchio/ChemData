#######################
## 0. GLOBAL IMPORTS
#######################
source("tools.r")
library(ggplot2)

######################################################################
## 1. COMPARISON OF LINEAR AND CONTINUOUS BERNOULLI MODELS
######################################################################

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

# Mixed terms with 2-level combinations
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


# Mixed terms with 3-level combinations
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

nc <- 28 
beta.hat.pls2 <- plsglm.cb.simple(xx.train, y, nc,scaling=TRUE)
gc()

eta.hat.pls2 <- as.vector(xx.train%*%beta.hat.pls2[-1]+beta.hat.pls2[1])
fit.pls2 <- kappa1(eta.hat.pls2)
ll.pls2 <- sum(y*eta.hat.pls2-kappa(eta.hat.pls2))

pred.pls2 <- kappa1(as.vector(xx.test%*%beta.hat.pls2[-1]+beta.hat.pls2[1]))
rsq.pls2 <- cor(pred.pls2,y.test)^2
rmse.pls2 <- sqrt(mean((pred.pls2*100-y.test*100)^2))

df.pls2 <- data.frame(y.test,pred.pls2)

pdf("plsglm_model2.pdf",width=10,height=10)
ggplot(data = df.pls2, aes(x = pred.pls2*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = "PLSGLM with 2 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-28),y=c(95,89),label=c(paste("R^2:",round(rsq.pls2,2)),
                                                  paste("RMSE:",round(rmse.pls2,2))), size=10)
dev.off()

coef.pls2=rbind(c("intercept",beta.hat.pls2[1]),cbind(colnames(xx),beta.hat.pls2[-1]))[order(abs(beta.hat.pls2)),]
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

sum(y*eta.hat2-kappa(eta.hat2))
rbind(c("intercept",beta.hat2[1]),cbind(colnames(xx),beta.hat2[-1]))[order(abs(beta.hat2)),]

############################
## PLSGLM on 3-level model
############################

memory.limit(32000)

nc=33
beta.hat.pls3 <- plsglm.cb.simple(xxx.train, y, nc,scaling=TRUE)
gc()

eta.hat.pls3 <- as.vector(xxx.train%*%beta.hat.pls3[-1]+beta.hat.pls3[1])
fit.pls3 <- kappa1(eta.hat.pls3)
ll.pls3 <- sum(y*eta.hat.pls3-kappa(eta.hat.pls3))

pred.pls3 <- kappa1(as.vector(xxx.test%*%beta.hat.pls3[-1]+beta.hat.pls3[1]))
rsq.pls3 <- cor(pred.pls3,y.test)^2
rmse.pls3 <- sqrt(mean((pred.pls3*100-y.test*100)^2))

df.pls3 <- data.frame(y.test,pred.pls3)

pdf("plsglm_model3.pdf",width=10,height=10)
ggplot(data = df.pls3, aes(x = pred.pls3*100, y = y.test*100)) + theme(text = element_text(size = 40))+
  labs(title = "PLSGLM with 3 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+
  geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text",x=c(-35,-28),y=c(95,89),label=c(paste("R^2:",round(rsq.pls3,2)),
                                                  paste("RMSE:",round(rmse.pls3,2))), size=10)
dev.off()

coef.pls3=rbind(c("intercept",beta.hat.pls3[1]),data.frame(colnames(xxx),beta.hat.pls3[-1]))[order(abs(beta.hat.pls3)),]
coef.pls3[2186:2196,]

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

#get all coefficients from the constraints
coef.one <- coef.pls3[2:41,]
for (i in 1:4) {ind=grepl(colnames(xs)[i],coef.one[,1],fixed=TRUE); coef.one=rbind(coef.one,c(paste(colnames(xs[i]),"last"),-sum(as.numeric(coef.one[ind,2]))))}

coef.two=coef.pls3[c(42:516),]
for (k in 1:4){
  for (j in 1:3)
for (i in (max(bb[k-1],0)+1):bb[k])
 { 
ind=grepl(paste(colnames(xf)[i],(colnames(xs)[-k])[j],sep=":"),coef.two[,1],fixed=TRUE); if(sum(ind)!=0) coef.two=rbind(coef.two,c(paste(paste(colnames(xf)[i],(colnames(xs)[-k])[j],sep=":"),"last"),-sum(as.numeric(coef.two[ind,2]))))
}
}

coef.two.rev=coef.pls3[c(42:516),]
ct.s=strsplit(coef.two.rev[,1], ":")
for (i in 1:length(ct.s)) {coef.two.rev[i,1]=paste(ct.s[[i]][2],ct.s[[i]][1],sep=":")}

for (k in 1:4){
  for (j in 1:3)
    for (i in (max(bb[k-1],0)+1):bb[k])
    { 
      ind=grepl(paste(colnames(xf)[i],(colnames(xs)[-k])[j],sep=":"),coef.two.rev[,1],fixed=TRUE); if(sum(ind)!=0) coef.two.rev=rbind(coef.two.rev,c(paste(paste(colnames(xf)[i],(colnames(xs)[-k])[j],sep=":"),"last"),-sum(as.numeric(coef.two.rev[ind,2]))))
    }
}
coef.two=rbind(coef.two,coef.two.rev[476:502,])

#need to add three interactions, but they seem to be small as well

coef.xxx=rbind(c("intercept",beta.hat[1]),coef.one,coef.two,coef.pls3[-c(1:516),])
(coef.xxx[order(abs(as.numeric(coef.xxx[,2]))),])[2300:2320,]

####################################
## add four levels combinations
####################################
xxx1=xxx[,-c(1:515)]

xxxx=rep(NA,nrow(xcf))
for (i in 1:21) {xxxxp=xcf[,i]*xxx1[,1597:1680];colnames(xxxxp)=paste(colnames(xcf)[i],colnames(xxx1)[1597:1680],sep=":"); xxxx=cbind(xxxx,xxxxp)}
xxxx=cbind(xxx,xxxx[,-1])

xxxx.train=xxxx[1:nrow(x),]
xxxx.test=xxxx[-c(1:nrow(x)),]

# fit GPLS
nc=29
beta.hat.pls4=plsglm.cb.simple(xxxx.train, y, nc,scaling=TRUE)
gc()
eta.hat.pls4=as.vector(xxxx.train%*%beta.hat.pls4[-1]+beta.hat.pls4[1])
fit.pls4=kappa1(eta.hat.pls4)

plot(y,fit.pls4)
abline(0,1,lwd=3,col=2)

pred.pls4=kappa1(as.vector(xxxx.test%*%beta.hat.pls4[-1]+beta.hat.pls4[1]))
cor(pred.pls4,y.test)^2
sqrt(mean((pred.pls4*100-y.test*100)^2))

df.pls4=data.frame(y.test,pred.pls4)

pdf("plsglm_model4.pdf",width=10,height=10)
ggplot(data = df.pls4, aes(x = pred.pls4*100, y = y.test*100)) + theme(text = element_text(size = 40)) +labs(title = "PLSGLM with 4 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text", x=c(-35,-28), y=c(95,89), label= c(expression(paste(R^2,"=0.15")),"RMSE=28.70"),size=10)
dev.off()

coef.pls4=rbind(c("intercept",beta.hat.pls4[1]),cbind(colnames(xxxx),beta.hat.pls4[-1]))[order(abs(beta.hat.pls4)),]
coef.pls4[3940:3960,]

# calculate the Log-likelihood
LL4=sum(y*eta.hat.pls4-kappa(eta.hat.pls4))


#find best number of components/ CV on test set with and without scaling
xxxx.train=xxxx[1:nrow(x),]
xxxx.test=xxxx[-c(1:nrow(x)),]
xxx.train=xxx[1:nrow(x),]
xxx.test=xxx[-c(1:nrow(x)),]
xx.train=xx[1:nrow(x),]
xx.test=xx[-c(1:nrow(x)),]
nc=(2:12)*3
cor2=rmse2=cor3=rmse3=cor4=rmse4=ll2=ll3=ll4=rep(NA,max(nc))
for (i in nc)
{
  print(paste("Checking",i,"comps..."))
  
  beta.hat.pls2=plsglm.cb.simple(xx.train, y, i,scaling=TRUE)
  gc()
  pred.cb.pls2=kappa1(as.vector(as.matrix(xx.test)%*%beta.hat.pls2[-1]+beta.hat.pls2[1]))
  cor2[i]=cor(pred.cb.pls2,y.test)
  rmse2[i]=sqrt(mean((pred.cb.pls2*100-y.test*100)^2))
  eta.hat.pls2=as.vector(xx.train%*%beta.hat.pls2[-1]+beta.hat.pls2[1])
  ll2[i]=sum(y*eta.hat.pls2-kappa(eta.hat.pls2))
  
  beta.hat.pls3=plsglm.cb.simple(xxx.train, y, i,scaling=TRUE)
  gc()
  pred.cb.pls3=kappa1(as.vector(as.matrix(xxx.test)%*%beta.hat.pls3[-1]+beta.hat.pls3[1]))
  cor3[i]=cor(pred.cb.pls3,y.test)
  rmse3[i]=sqrt(mean((pred.cb.pls3*100-y.test*100)^2))
  eta.hat.pls3=as.vector(xxx.train%*%beta.hat.pls3[-1]+beta.hat.pls3[1])
  ll3[i]=sum(y*eta.hat.pls3-kappa(eta.hat.pls3))
  
  beta.hat.pls4=plsglm.cb.simple(xxxx.train, y, i,scaling=TRUE)
  gc()
  pred.cb.pls4=kappa1(as.vector(as.matrix(xxxx.test)%*%beta.hat.pls4[-1]+beta.hat.pls4[1]))
  cor4[i]=cor(pred.cb.pls4,y.test)
  rmse4[i]=sqrt(mean((pred.cb.pls4*100-y.test*100)^2))
  eta.hat.pls4=as.vector(xxxx.train%*%beta.hat.pls4[-1]+beta.hat.pls4[1])
  ll4[i]=sum(y*eta.hat.pls4-kappa(eta.hat.pls4))

}

plot(c(1,max(nc)),c(0,1),col="white",main="CORR")
points(cor2,col=2)
points(cor3,col=3)
points(cor4,col=4)
legend("topleft", col=c(2,3,4), lwd=2, lty=1, cex=1, c("2-level","3-level","4-level") )

plot(c(1,max(nc)),c(0,50),col="white",main="RSME")
points(rmse2,col=2)
points(rmse3,col=3)
points(rmse4,col=4)
legend("topleft", col=c(2,3,4), lwd=2, lty=1, cex=1, c("2-level","3-level","4-level") )

write(cor2, file="cor2.txt", sep="\n")
write(cor3, file="cor3.txt", sep="\n")
write(cor4, file="cor4.txt", sep="\n")

write(rmse2, file="rmse2.txt", sep="\n")
write(rmse3, file="rmse3.txt", sep="\n")
write(rmse4, file="rmse4.txt", sep="\n")


