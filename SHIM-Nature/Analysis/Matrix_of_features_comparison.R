#######################
## Read the data
#######################

x=data.matrix(read.table("Data/Xtrain.txt",header=TRUE))
x.test=data.frame(read.table("Data/Xtest.txt",header=TRUE))
y=as.vector(read.table("Data/YTrain.txt",header=FALSE)[,1])
y.test=as.vector(read.table("Data/YTest.txt",header=FALSE)[,1])

##############################################################
## Compare linear models with dummy coding and chemical descriptors
#############################################################

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


