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
ggplot(df, aes(x=Yield)) + 
  geom_histogram(aes(y=..density..),color="black", fill="cornflowerblue",breaks=seq(0,1,length=11))+ theme(text = element_text(size = 40)) +labs(title = "Yield density") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_line(aes(x=xx,y =ll), color = "black",size=2)+ylab("Density")+xlab("Yield")
dev.off()

############################
## Plot the link function g
############################

# Generate data
xx <- seq(0.065, 0.935, length.out = 1000)
g.link <- function(x) {
  tt <- apply(as.matrix(x), 1, FUN=function(v) min(max(v,0.001),0.999))
  g <- 3.5*tan(pi*(2*tt-1)/2) 
  return(g)
}
df <- data.frame(x = xx, y = g.link(xx))

# Create the plot
pdf("link_function_g.pdf", width = 10, height = 10)
ggplot(df, aes(x = x , y = y)) +
  geom_line(color = "cornflowerblue", size = 2) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),  # Force tick marks
    limits = c(0, 1)                    # Force visual limits
  ) +scale_y_continuous(
    breaks = c(-10,0,10),  # Force tick marks
    limits = c(-17,17)                    # Force visual limits
  ) +
  theme(text = element_text(size = 40)) +
  labs(title = "Link function") +
  theme(
    plot.title = element_textbox(
      hjust = 0.5,
      margin = margin(t = 5, b = 5)
    )
  ) +
  ylab("Systematic Component") +
  xlab("Expected Yield")
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
kappa1<-function(x) 1/(1-exp(-x))-1/x
kappa2<-function(x) 1/x^2-exp(x)/(1-exp(x))^2

fisher.scoring<-function(y,x,initial){
  beta0<-initial
  x=cbind(rep(1,nrow(x)),x)
  counter=0
  repeat{
    counter=counter+1
    eta=as.vector(x%*%beta0)
    W=kappa2(eta)
    W1=diag(1/kappa2(eta))
    Z=x%*%beta0+W1%*%(y-kappa1(eta))
    fit<-lm(Z~x-1,weights=W)
    beta1=fit$coef
    epsilon <- sqrt(sum((beta0-beta1)^2)/sum(beta0^2))
    if(epsilon<=1e-6) break
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
ggplot(data = df.lm, aes(x = pred.cb*100, y = y.test*100)) + theme(text = element_text(size = 40)) +labs(title = "Generalised linear model") +
  theme(plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5)))+
  geom_point(color="cornflowerblue",size=3)+xlab("Predicted Yield")+ylab("Observed Yield")+geom_abline(intercept = 0, slope = 1, size = 1.5,linetype = "dashed")+xlim(-50,100)+
  annotate("text", x=c(-35,-28), y=c(95,89), label= c(expression(paste(R^2,"=0.78")),"RMSE=12.6"),size=10)
dev.off()

##############################################################
## Compare linear models with dummy coding and chemical descriptors
#############################################################

x=data.matrix(read.table("Data/Xtrain.txt",header=TRUE))
x.test=data.frame(read.table("Data/Xtest.txt",header=TRUE))

#add three artificial descriptors for additives
xc=rbind(x,x.test)
set.seed(2134)
c1=factor(xc[,1])
c2=factor(xc[,3])
c3=factor(xc[,4])
levels(c1)=runif(22)
levels(c2)=rnorm(22)
levels(c3)=runif(22)
xc=cbind(scale(cbind(as.numeric(c1),as.numeric(c2),as.numeric(c3))),xc)
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
yy=c(y,y.test)
fit1=lm(yy~as.matrix(xc))
fit2=lm(yy~xcf)
plot(fitted(fit1),fitted(fit2),type="l")

#compare to the model with 19 descriptors

fit3=lm(yy~as.matrix(rbind(x,x.test)))
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

