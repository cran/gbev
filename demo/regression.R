

cat(" Linear regression example \n")
n<-1000
varX<-1
varME<-.5
x1<-rnorm(n,sd=sqrt(varX))
x2<-rnorm(n,sd=sqrt(varX))
w1<-x1+rnorm(n,sd=sqrt(varME))
w2<-x2+rnorm(n,sd=sqrt(varME))
y<-x1+2*x2+rnorm(n,sd=.3)

dat<-data.frame(y=y,w1=w1,w2=w2,z1=rnorm(n),z2=rnorm(n),z3=rnorm(n))

cat(" Specify measurement error model \n")
pME<-2                         ## 2 covariates measured with error, the rest error free

numComp<-1                     ## number of components in gaussian mixture for X-distribution
SigmaME<-diag(varME,pME)
SigmaX<-array(dim=c(numComp,pME,pME))         
mu<-array(dim=c(numComp,pME))
pComp<-array(1/numComp,dim=c(numComp,1))
for(i in 1:numComp)
{
SigmaX[i,,]<-diag(varX,pME)
mu[i,]<-rep(0,pME)
}
### list required by "gbev" for measurement error model
meModel<-list(SigmaX=SigmaX,mu=mu,SigmaME=SigmaME,pComp=pComp,numComp=numComp)

fit<-gbev(y~w1+w2+z1+z2+z3,data=dat,
          measErrorModel=meModel,     
          method="L2",              ## Squared error loss
          nboost=2000,              ## 1000 boosting iterations
          lambda=5,                 ## regularization of regression tree
          maxDepth=2,               ## maximum tree depth, 2 corresponds stumps
          mc=2,                     ## number of monte-carlo samples per tree build 
          minSplit=3,               ## minimum number of obs in node to split
          minBucket=0,              ## minimum number of obs in nodes
          sPoints=1,                ## number of sampled candidate split points
          intermPred=5)             ## increments of iterations to store predictions 


### 5-fold cross-validation
hcv<-cvLoss(object=fit,k=5,random=FALSE,loss="L2")
plot(hcv$iters,hcv$cvLoss,type="l")

##### Plot estimated marginal effect of covariates on response ####
split.screen(c(1,2))
screen(1)
hp1<-part.dep(fit,varIndx=1,ngrid=200,plt=FALSE,firstTree=1,lastTree=hcv$estIter)
x<-hp1$x[hp1$x<2&hp1$x>-2]
f1<-hp1$pred[hp1$x<2&hp1$x>-2]
plot(x,f1,type="l")
abline(0,1)                           ## true regression line between "y" and "x1"
abline(0,1*(varX/(varX+varME)),lty=5) ## true regression line between "y" and "w1"

screen(2)
hp2<-part.dep(fit,varIndx=2,ngrid=200,plt=FALSE,firstTree=1,lastTree=hcv$estIter)
x<-hp2$x[hp2$x<2&hp2$x>-2]
f2<-hp2$pred[hp2$x<2&hp2$x>-2]
plot(x,f2,type="l")
abline(0,2)                            ## true regression line between "y" and "x2"
abline(0,2*(varX/(varX+varME)),lty=5)  ## true regression line between "y" and "w2"

close.screen(all=T)

## Plot variable importance
var.imp(fit)


