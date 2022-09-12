#######################
## Read the data
#######################

x=data.matrix(read.table("Data/Xtrain.txt",header=TRUE))
x.test=data.frame(read.table("Data/Xtest.txt",header=TRUE))
y=as.vector(read.table("Data/YTrain.txt",header=FALSE)[,1])
y.test=as.vector(read.table("Data/YTest.txt",header=FALSE)[,1])

#################################################
## Fit continuous Bernoulli and plot histogram
#################################################
yy=c(y,y.test)
nn=length(yy)
xx=1:nn/nn

#find MLE for p
ld<-function(p)
{
  (nn*p/(2*p-1)-sum(yy)+nn/log((1-p)/p))/(p*(p-1))
}
p.hat=uniroot(ld,c(0.001,0.5))$root

# Plot the historgram
library(ggplot2)
ll=log((1-p.hat)/p.hat)*p.hat^xx*(1-p.hat)^(1-xx)/(1-2*p.hat)
df=data.frame(yy,xx,ll)
colnames(df)=c("Yield","Probability","Density")

### functions for the textbox in ggplot

element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}

pdf("yield_hist.pdf",width=10,height=10)
ggplot(df, aes(x=Yield*100)) + 
  geom_histogram(aes(y=..density..),color="black", fill="cornflowerblue",breaks=seq(0,100,length=11))+ theme(text = element_text(size = 40)) +labs(title = "Yield density") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_line(aes(x=xx*100,y =ll/100), color = "black",size=2)+ylab("Density")+xlab("Yield")
dev.off()


########################
## Run linear regression
########################

#remove linear dependent columns
ind=qr(x)$pivot[seq_len(qr(x)$rank)]
x=x[, ind]
x.test=x.test[,ind]
ind1=which(is.na(as.vector(lm(y~x)$coef)))-1
x=x[,-ind1]
x.test=x.test[,-ind1] 

fit.lm=lm(y~x)
pred.lm=as.matrix(x.test)%*%fit.lm$coef[-1]+fit.lm$coef[1]
cor(pred.lm,y.test)^2
sqrt(mean((pred.lm*100-y.test*100)^2))

df.lm=data.frame(y.test,pred.lm)
colnames(df.lm)=c("ObservedYield","PredictedYield")

pdf("linear_model.pdf",width=10,height=10)
ggplot(data = df.lm, aes(x = PredictedYield*100, y = ObservedYield*100)) + theme(text = element_text(size = 40)) +labs(title = "Linear model") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text", x=c(-35,-28), y=c(95,89), label= c(expression(paste(R^2,"=0.69")),"RMSE=15.1"),size=10)
dev.off()

#################################################
### Regression with continuous Bernoulli response
#################################################

#Fisher scoring
kappa1 <- function(x) {
  k1 <- rep(1/2,length(x))
  tt <- which(x!=0)
  k1[tt] <- 1 - 1/(x[tt]) - 1/(1-exp(x[tt]))
  return(k1)
}

kappa2 <- function(x) {
  k2 <- rep(1/12,length(x))
  tt <- which(x!=0)
  k2[tt] <- 1/(x[tt])^2 + 1/(2-2*cosh(x[tt]))
  return(k2)
}

fisher.scoring<-function(y,x,initial){
  beta0<-initial
  x=cbind(rep(1,nrow(x)),x)
  counter=0
  repeat{
    counter=counter+1
    eta=as.vector(x%*%beta0)
    
    kp=kappa1(eta)
    kpp=kappa2(eta)
  
    kpp[which(kpp<0.01)] <- 0.01
    
    W=kpp
    Z=x%*%beta0+(y-kp)/W
    fit<-lm(Z~x-1,weights=W)
    beta1=fit$coef
    epsilon <- sqrt(sum((beta0-beta1)^2)/sum(beta0^2))
    print(paste("Epsilon: ", epsilon, sep=""))
    if(epsilon<=0.05) break
    if(counter==100) {print("no convergence"); break}
    beta0 <- beta1
  }
  list(beta.hat=as.vector(beta1),zstat=summary(fit)[["coefficients"]][, "t value"])
}

#set starting value
beta.start=lm(y~x)$coef

#estimate beta
beta.hat=fisher.scoring(y,x,beta.start)$beta.hat
eta.hat=as.vector(x%*%beta.hat[-1]+beta.hat[1])
fit=kappa1(eta.hat)

pred.cb=kappa1(as.vector(as.matrix(x.test)%*%beta.hat[-1]+beta.hat[1]))
cor(pred.cb,y.test)^2
sqrt(mean((pred.cb*100-y.test*100)^2))

df.glm=data.frame(y.test,pred.cb)

pdf("cb_model.pdf",width=10,height=10)
ggplot(data = df.glm, aes(x = pred.cb*100, y = y.test*100)) + theme(text = element_text(size = 40)) +labs(title = "Generalised linear model") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text", x=c(-35,-28), y=c(95,89), label= c(expression(paste(R^2,"=0.78")),"RMSE=12.6"),size=10)
dev.off()

##############################################################
## Compare linear models with dummy coding and chemical descriptors
#############################################################

x=data.matrix(read.table("Data/XtrainNScaled.txt",header=TRUE))
x.test=data.frame(read.table("Data/XtestNScaled.txt",header=TRUE))

#add three artificial descriptors for additives
xc=rbind(x,x.test)
set.seed(2134)
c1=factor(xc[,1])
c2=factor(xc[,3])
c3=factor(xc[,4])
levels(c1)=runif(22)
levels(c2)=rnorm(22)
levels(c3)=runif(22)
xc=cbind(cbind(as.numeric(c1),as.numeric(c2),as.numeric(c3)),xc)
colnames(xc)[1:3]=c("add_new1","add_new2","add_new3")

xs=xc[,c(1,23,50,60)]
colnames(xs)=c("additive","aryl_halide","base","ligand")

a=rep(NA,ncol(xs))
for (i in 1:ncol(xs)) a[i]=length(unique(xs[,i]))
xcf=matrix(NA,nrow(xs),sum(a))
b=cumsum(a)
colnames(xcf)=rep(colnames(xs),times=a)
colnum=order(unique(xs[,1]))
for (i in 2:ncol(xs)) colnum=c(colnum,order(unique(xs[,i])))
colnames(xcf)=paste(colnames(xcf),colnum)

for (i in 1:nrow(xs))
{
  for (j in 1:length(a))
  {
    res <- rep(0, a[j])
    where <- match( xs[i,j], unique(xs[,j]) )
    res[ where ] <- 1 
    xcf[i,(max(b[j-1],0)+1):b[j]]=res
  }
}

for (i in 1:length(b))
{
  ind=match(xcf[,b[i]],1)==1
  xcf[ind,(max(b[i-1],0)+1):b[i]]=-1
}

xcf=xcf[,-b]

xf=xcf[1:nrow(x),]
xf.test=xcf[-c(1:nrow(x)),]


#fit linear models with dummy coding and descriptors matrix

fit1=lm(yy~as.matrix(scale(xc)))
fit2=lm(yy~xcf)
plot(fitted(fit1),fitted(fit2),type="l")

#compare to the model with 19 descriptors

fit3=lm(yy~scale(as.matrix(rbind(x,x.test))))
plot(fitted(fit1),fitted(fit3),type="l")

################################################
### Compare to random forest
################################################

library(randomForest)

set.seed(1234)
#RF with 120 correlated predictors
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

################################################
### mixed terms with continuous Bernoulli
################################################

########################################
## generalized PLS function
source("tools.r")
plsglm.cb.simple <- function(x, y, m.plsglm, 
                            beta0=NULL,
                            centering=TRUE, scaling=TRUE, intercept=TRUE,
                            maxit=20, tol=0.05,
                            verbose=FALSE){
  
  return(plsglm.cb(X=x, Y=y, ncomp=m.plsglm, 
                       beta0=beta0,
                       centering=centering, scaling=scaling, intercept=intercept,
                       maxit=maxit, tol=tol,
                       verbose=verbose)$BETA[,,m.plsglm])
}
#######################################

###############################################
#build mixed terms with two level combinations
################################################

xx=rep(1,nrow(xcf))
bb=cumsum(a-1)
for (j in 1:3)
  for (i in (max(bb[j-1],0)+1):bb[j]) {xxp=xcf[,i]*xcf[,-c(1:bb[j])]; colnames(xxp)=paste(colnames(xcf)[i],colnames(xcf[,-c(1:bb[j])]),sep=":");  xx=cbind(xx,xxp)}
xx=cbind(xcf,xx[,-1])


#fit GPLS
nc=28 
beta.hat.pls2=plsglm.cb.simple(xx, yy, nc,scaling=TRUE)
gc()
eta.hat.pls2=as.vector(xx%*%beta.hat.pls2[-1]+beta.hat.pls2[1])
fit.pls2=kappa1(eta.hat.pls2)

plot(yy,fit.pls2)
abline(0,1,lwd=3,col=2)

pred.pls2 <- fit.pls2[-(1:length(y))]
cor(pred.pls2,y.test)^2
sqrt(mean((pred.pls2*100-y.test*100)^2))

df.pls2=data.frame(y.test,pred.pls2)

pdf("plsglm_model2.pdf",width=10,height=10)
ggplot(data = df.pls2, aes(x = pred.pls2*100, y = y.test*100)) + theme(text = element_text(size = 40)) +labs(title = "PLSGLM with 2 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text", x=c(-35,-28), y=c(95,89), label= c(expression(paste(R^2,"=0.92")),"RMSE=7.55"),size=10)
dev.off()

coef.pls2=rbind(c("intercept",beta.hat.pls2[1]),cbind(colnames(xx),beta.hat.pls2[-1]))[order(abs(beta.hat.pls2)),]
coef.pls2[496:516,]

#fit GLM
beta.hat=fisher.scoring(yy,xx,beta.hat.pls2)$beta.hat #has convergence problems due to many beta=0 weights become close to 0
eta.hat=as.vector(xx%*%beta.hat[-1]+beta.hat[1])
fit=kappa1(eta.hat)

plot(yy,fit)
abline(0,1,lwd=3,col=2)

pred.hat2 <- fit[-(1:length(y))]
cor(pred.hat2,y.test)^2
sqrt(mean((pred.hat2*100-y.test*100)^2))

df.hat2=data.frame(y.test,pred.hat2)

pdf("cb_model2.pdf",width=10,height=10)
ggplot(data = df.hat2, aes(x = pred.hat2*100, y = y.test*100)) + theme(text = element_text(size = 40)) +labs(title = "GLM with 2 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text", x=c(-35,-28), y=c(95,89), label= c(expression(paste(R^2,"=0.92")),"RMSE=7.54"),size=10)
dev.off()

sum(yy*eta.hat-kappa(eta.hat))
rbind(c("intercept",beta.hat[1]),cbind(colnames(xx),beta.hat[-1]))[order(abs(beta.hat)),]

# calculate the Log-likelihood
kappa<-function(x) log((exp(x)-1)/x)
LL2=sum(yy*eta.hat.pls2-kappa(eta.hat.pls2))

#############################
#add three level combinations
###############################

xcf1=xcf
colnames(xcf1)=c(rep("additive",21),rep("aryl_halide",14),rep("base",2),rep("ligand",3))
xx1=xx[,-c(1:40)]

xxx=rep(1,nrow(xcf))
ind=rep(TRUE,ncol(xx1))
for (j in 1:2) 
{
  ind=as.logical((!grepl(colnames(xcf1)[bb[j]],colnames(xx1)))*(ind));
  for (i in (max(bb[j-1],0)+1):bb[j]) {xxxp=xcf[,i]*xx1[,ind]; colnames(xxxp)=paste(colnames(xcf)[i],colnames(xx1)[ind],sep=":"); xxx=cbind(xxx,xxxp)
  }
}
xxx=cbind(xx,xxx[,-1])

# fit GPLS
memory.limit(32000)
nc=33
beta.hat.pls3=plsglm.cb.simple(xxx, yy, nc,scaling=TRUE)
gc()
eta.hat.pls3=as.vector(xxx%*%beta.hat.pls3[-1]+beta.hat.pls3[1])
fit.pls3=kappa1(eta.hat.pls3)

plot(yy,fit.pls3)
abline(0,1,lwd=3,col=2)

pred.pls3 <- fit.pls3[-(1:length(y))]
cor(pred.pls3,y.test)^2
sqrt(mean((pred.pls3*100-y.test*100)^2))

df.pls3=data.frame(y.test,pred.pls3)

pdf("plsglm_model3.pdf",width=10,height=10)
ggplot(data = df.pls3, aes(x = pred.pls3*100, y = y.test*100)) + theme(text = element_text(size = 40)) +labs(title = "PLSGLM with 3 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text", x=c(-35,-28), y=c(95,89), label= c(expression(paste(R^2,"=0.98")),"RMSE=3.36"),size=10)
dev.off()

coef.pls3=rbind(c("intercept",beta.hat.pls3[1]),data.frame(colnames(xxx),beta.hat.pls3[-1]))[order(abs(beta.hat.pls3)),]
coef.pls3[2186:2196,]

# calculate the Log-likelihood
LL3=sum(yy*eta.hat.pls3-kappa(eta.hat.pls3))

#fit GLM
beta.hat=fisher.scoring(yy,xxx,beta.hat.pls3)$beta.hat #same convergence problems
gc()
eta.hat=as.vector(xxx%*%beta.hat[-1]+beta.hat[1])
fit=kappa1(eta.hat)
rbind(c("intercept",beta.hat[1]),data.frame(colnames(xxx),beta.hat[-1]))[order(abs(beta.hat)),]

plot(yy,fit)
abline(0,1,lwd=3,col=2)

pred.hat3 <- fit[-(1:length(y))]
cor(pred.hat3,y.test)^2
sqrt(mean((pred.hat3*100-y.test*100)^2))

df.hat3=data.frame(y.test,pred.hat3)

pdf("cb_model3.pdf",width=10,height=10)
ggplot(data = df.hat3, aes(x = pred.hat3*100, y = y.test*100)) + theme(text = element_text(size = 40)) +labs(title = "GLM with 3 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text", x=c(-35,-28), y=c(95,89), label= c(expression(paste(R^2,"=0.99")),"RMSE=2.99"),size=10)
dev.off()

#get all coefficients from the constraints
coef.one=coef.pls3[2:41,]
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

# fit GPLS
nc=29
beta.hat.pls4=plsglm.cb.simple(xxxx, yy, nc,scaling=TRUE)
gc()
eta.hat.pls4=as.vector(xxxx%*%beta.hat.pls4[-1]+beta.hat.pls4[1])
fit.pls4=kappa1(eta.hat.pls4)

plot(yy,fit.pls4)
abline(0,1,lwd=3,col=2)

pred.pls4 <- fit.pls4[-(1:length(y))]
cor(pred.pls4,y.test)^2
sqrt(mean((pred.pls4*100-y.test*100)^2))

df.pls4=data.frame(y.test,pred.pls4)

pdf("plsglm_model4.pdf",width=10,height=10)
ggplot(data = df.pls4, aes(x = pred.pls4*100, y = y.test*100)) + theme(text = element_text(size = 40)) +labs(title = "PLSGLM with 4 levels") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text", x=c(-35,-28), y=c(95,89), label= c(expression(paste(R^2,"=0.98")),"RMSE=3.36"),size=10)
dev.off()

coef.pls4=rbind(c("intercept",beta.hat.pls4[1]),cbind(colnames(xxxx),beta.hat.pls4[-1]))[order(abs(beta.hat.pls4)),]
coef.pls4[3940:3960,]

# calculate the Log-likelihood
LL4=sum(yy*eta.hat.pls4-kappa(eta.hat.pls4))


#find best number of components/ CV on test set with and without scaling
xxxx.train=xxxx[1:nrow(x),]
xxxx.test=xxxx[-c(1:nrow(x)),]
xxx.train=xxx[1:nrow(x),]
xxx.test=xxx[-c(1:nrow(x)),]
xx.train=xx[1:nrow(x),]
xx.test=xx[-c(1:nrow(x)),]
nc=42
cor2=rmse2=cor3=rmse3=cor4=rmse4=ll2=ll3=ll4=cor2S=rmse2S=cor3S=rmse3S=cor4S=rmse4S=ll2S=ll3S=ll4S=rep(NA,nc)
for (i in 5:nc)
{
  beta.hat.pls2=plsglm.cb.simple(xx.train, y, i,scaling=FALSE)
  pred.cb.pls2=kappa1(as.vector(as.matrix(xx.test)%*%beta.hat.pls2[-1]+beta.hat.pls2[1]))
  cor2[i]=cor(pred.cb.pls2,y.test)
  rmse2[i]=sqrt(mean((pred.cb.pls2*100-y.test*100)^2))
  eta.hat.pls2=as.vector(xx.train%*%beta.hat.pls2[-1]+beta.hat.pls2[1])
  ll2[i]=sum(y*eta.hat.pls2-kappa(eta.hat.pls2))
  
  beta.hat.pls3=plsglm.cb.simple(xxx.train, y, i,scaling=FALSE)
  pred.cb.pls3=kappa1(as.vector(as.matrix(xxx.test)%*%beta.hat.pls3[-1]+beta.hat.pls3[1]))
  cor3[i]=cor(pred.cb.pls3,y.test)
  rmse3[i]=sqrt(mean((pred.cb.pls3*100-y.test*100)^2))
  eta.hat.pls3=as.vector(xxx.train%*%beta.hat.pls3[-1]+beta.hat.pls3[1])
  ll3[i]=sum(y*eta.hat.pls3-kappa(eta.hat.pls3))
  
  beta.hat.pls4=plsglm.cb.simple(xxxx.train, y, i,scaling=FALSE)
  pred.cb.pls4=kappa1(as.vector(as.matrix(xxxx.test)%*%beta.hat.pls4[-1]+beta.hat.pls4[1]))
  cor4[i]=cor(pred.cb.pls4,y.test)
  rmse4[i]=sqrt(mean((pred.cb.pls4*100-y.test*100)^2))
  eta.hat.pls4=as.vector(xxxx.train%*%beta.hat.pls4[-1]+beta.hat.pls4[1])
  ll4[i]=sum(y*eta.hat.pls4-kappa(eta.hat.pls4))
  
  
  beta.hat.pls2=plsglm.cb.simple(xx.train, y, i)
  pred.cb.pls2=kappa1(as.vector(as.matrix(xx.test)%*%beta.hat.pls2[-1]+beta.hat.pls2[1]))
  cor2S[i]=cor(pred.cb.pls2,y.test)
  rmse2S[i]=sqrt(mean((pred.cb.pls2*100-y.test*100)^2))
  eta.hat.pls2=as.vector(xx.train%*%beta.hat.pls2[-1]+beta.hat.pls2[1])
  ll2S[i]=sum(y*eta.hat.pls2-kappa(eta.hat.pls2))
  
  beta.hat.pls3=plsglm.cb.simple(xxx.train, y, i)
  pred.cb.pls3=kappa1(as.vector(as.matrix(xxx.test)%*%beta.hat.pls3[-1]+beta.hat.pls3[1]))
  cor3S[i]=cor(pred.cb.pls3,y.test)
  rmse3S[i]=sqrt(mean((pred.cb.pls3*100-y.test*100)^2))
  eta.hat.pls3=as.vector(xxx.train%*%beta.hat.pls3[-1]+beta.hat.pls3[1])
  ll3S[i]=sum(y*eta.hat.pls3-kappa(eta.hat.pls3))
  
  beta.hat.pls4=plsglm.cb.simple(xxxx.train, y, i)
  pred.cb.pls4=kappa1(as.vector(as.matrix(xxxx.test)%*%beta.hat.pls4[-1]+beta.hat.pls4[1]))
  cor4S[i]=cor(pred.cb.pls4,y.test)
  rmse4S[i]=sqrt(mean((pred.cb.pls4*100-y.test*100)^2))
  eta.hat.pls4=as.vector(xxxx.train%*%beta.hat.pls4[-1]+beta.hat.pls4[1])
  ll4S[i]=sum(y*eta.hat.pls4-kappa(eta.hat.pls4))
}

