# source("C:/fun/gbev/R/gbev.R")



.First.lib <- function(lib, pkg)
{
     library.dynam("gbev", pkg, lib)
     require(mvtnorm)

}



gbev<-function(formula = formula(data),
                     data = list(),
                     weights=NULL,
                     measErrorModel=NULL,
                     method="L2", 
                     indepFitUpdate=1,
                     nboost=100,
                     lambda=100,
                     maxDepth=2,
                     minSplit=10, 
                     minBucket=0,
                     sPoints=10,
                     mc=2,
                     intermPred=10,
                     maxSplitAttempts=10)

{
  
if(sum(is.na(data))>0)
stop(" 'gbev' does not handle NA's \n");



nBoosts<-nboost


   call <- match.call()
   m <- match.call(expand = FALSE)
 
   m$measErrorModel<- m$method <-m$indepFitUpdate<- m$nboost<- m$lambda <- NULL
   m$maxDepth <- NULL
   m$minSplit <- m$minBucket <-m$sPoints<- m$mc <-m$intermPred <-m$maxSplitAttempts <- NULL
  
   m[[1]] <- as.name("model.frame")
   m$na.action <- na.pass
   m.keep <- m

   m <- eval(m, parent.frame())
  
   Terms <- attr(m, "terms")
   a <- attributes(Terms)

   y <- model.extract(m, response)

   x <- model.frame(delete.response(Terms),
                    data,
                    na.action=na.pass)

   names.covar<-names(x)

   name.response <- dimnames(attr(terms(formula),"factors"))[[1]][1]

   var.type <- rep(0,ncol(x))
   var.levels <- vector("list",ncol(x))

   for(i in 1:ncol(x))
   {
      if(is.ordered(x[,i]))
      {
         var.levels[[i]] <- levels(x[,i])
         x[,i] <- as.numeric(x[,i])-1
         var.type[i] <- 0
      }
      else if(is.factor(x[,i]))
      {
         var.levels[[i]] <- levels(x[,i])
         x[,i] <- as.numeric(x[,i])-1
         var.type[i] <- max(x[,i],na.rm=TRUE)+1
      }
      else if(is.numeric(x[,i]))
      {
         var.levels[[i]] <- quantile(x[,i],prob=(0:10)/10,na.rm=TRUE)
      }
      else
      {
         stop("variable ",i,": ",names.covar[i]," is not of type numeric, ordered, or factor.")
      }
      
      # check for some variation in each variable
      if(length(unique(var.levels[[i]])) == 1)
      {
         warning("variable ",i,": ",var.names[i]," has no variation.")
      }
   }

## "ms" can be used to sub-sample number of covariates used to split a node. 
##      here it is set to all covariates 
ms<-length(x[1,]) 


#print(names.covar)
#print(name.response)
x<-as.matrix(x)

## Check to see if x is numeric
## 
for(i in 1:ncol(x)){
if(!is.numeric(x[,i]))
stop(" 'gbev' only handles numeric covariates. \n");
}


hfit<-gbev.fit(w=x,y=y,weights=weights,measErrorModel=measErrorModel,method=method, 
                     indepFitUpdate=indepFitUpdate,nboost=nboost,lambda=lambda,
                     maxDepth=maxDepth,
                     m=ms, minSplit=minSplit,minBucket=minBucket,sPoints=sPoints,mc=mc,
                     intermPred=intermPred,maxSplitAttempts=maxSplitAttempts)

tn<-hfit$totNumNodes

numNodeParam<-4  # number of node-parameters (i.e doubles)
numNodePosit<-5  # number of node-positions  (i.e integers)


h<-list(nTrees=nBoosts,totNumNodes=hfit$totNumNodes,oobError=hfit$oobError,
    treeStart=hfit$treeStart,treeEnd=hfit$treeEnd,nodeParameters=hfit$nodeParameters[1:(numNodeParam*tn)],
    nodePositions=hfit$nodePositions[1:(numNodePosit*tn)],iPred=hfit$iPred,recordPred=hfit$recordPred,
    namesCovar=names.covar,nameResponse=name.response,Terms=Terms)



gbev.obj<-h
gbev.obj$nboost<-nboost
gbev.obj$call<-call
gbev.obj$measErrorModel<-measErrorModel
gbev.obj$varLevels<-var.levels
gbev.obj$varType<-var.type
gbev.obj$ncovar<-length(names.covar)
gbev.obj$method<-method

gbev.obj$indepFitUpdate<-indepFitUpdate
gbev.obj$lambda<-lambda

gbev.obj$maxDepth<-maxDepth
gbev.obj$ms<-ms
gbev.obj$minSplit<-minSplit
gbev.obj$minBucket<-minBucket
gbev.obj$sPoints<-sPoints
gbev.obj$mc<-mc
gbev.obj$y<-y
gbev.obj$w<-x

gbev.obj$weights<-weights
gbev.obj$intermPred<-intermPred
gbev.obj$n<-length(y)
gbev.obj$ncovar<-length(names.covar)
gbev.obj$maxSplitAttempts<-maxSplitAttempts


class(gbev.obj)<-"gbev"

 


return(gbev.obj)

}





gbev.fit<-function(w,y,
                     weights=NULL,
                     measErrorModel=NULL,
                     method="L2", 
                     indepFitUpdate=1,
                     nboost=1000,
                     lambda=100,
                     maxDepth=2,
                     m=1,
                     minSplit=10, 
                     minBucket=0,
                     sPoints=10,
                     mc=2,
                     intermPred=10,
                     maxSplitAttempts=10)
{

nBoosts<-nboost


if(is.null(measErrorModel))
{
   stop("Need to specify measurement error model.")
}


checkMeasErrorModel(measErrorModel=measErrorModel);



n<-length(y)
p<-ncol(w)
pME<-length(measErrorModel$mu[1,])

if(pME>p)
{
   stop("Dimension of measurement error model larger than number of 'numeric' covariates.")
}


if(!((indepFitUpdate==0)|(indepFitUpdate==1)))
{
   stop("indepFitUpdate must be 0 or 1")
} 

if(nrow(w)!=n)
{
   stop("Number of rows in w doesn't match number of y's")
}

if(as.integer(mc)<1)
{
   stop("'mc' must be 1 or larger ")
}

if(is.null(weights))
{
weights<-rep(1,n)
}

   if(maxDepth < 2)
   {
      stop("maxDepth must be at least 2, corresponds a single split.")
   }

   if(minBucket > minSplit/2)
   {
      stop("minSplit must be atleast 2 time minBucket.")
   }

   supported.methods <-
      c("L2","logLike")  
   if(!is.element(method,supported.methods))
   {
      stop("Method ",method," is not supported")
   }
    
   if(method=="L2")
   {
    splitFunction<-1    
   }
   if(method=="logLike")
   {
    splitFunction<-2  
    }
   
  



               
 
wVec<-matrix(w,nrow=n*p,ncol=1)
maxNodes<-(2^(maxDepth)-1)*nBoosts

#
#########################################
### COMPUTE CONDITIONAL DENSITIES #######
#########################################
 
if(ncol(w)>1)
{
h<-mixModelMatrices(w=as.matrix(w[,1:pME]),numComp=measErrorModel$numComp,
                    mu=measErrorModel$mu,SigmaME=measErrorModel$SigmaME,
                    SigmaJ=measErrorModel$SigmaX,pComp=measErrorModel$pComp)
}
if(ncol(w)==1)
{
h<-mixModelMatrices(w=w,numComp=measErrorModel$numComp,
                    mu=measErrorModel$mu,SigmaME=measErrorModel$SigmaME,
                    SigmaJ=measErrorModel$SigmaX,pComp=measErrorModel$pComp)
}

 

if(is.null(intermPred))
{
recordBoosts<-nBoosts-1 ## only final fitted values are stored
}
if(!is.null(intermPred))
{
recordBoosts<-seq(1,nBoosts,intermPred) ## will record all predictions at iterations given by "recordBoosts"
}
nr<-length(recordBoosts)

ranSeed<-9583939 ## this has no effect, since now use R's RNG-functions.
#print(Sys.time());

hfit<-.C("gbevLWC",as.double(y),as.double(wVec[,1]),as.double(weights),
       as.integer(n),as.integer(p),as.integer(pME),as.integer(nboost),
       as.double(lambda*n),as.integer(maxDepth),as.integer(m),as.integer(minSplit),
       as.integer(minBucket),as.integer(sPoints),as.integer(splitFunction),
       as.integer(indepFitUpdate),as.integer(mc),
       as.integer(h$numComp),as.double(h$compProb),as.double(h$condMean),
       as.double(h$cholCov),numNodes=integer(1),nodePositions=integer(maxNodes*6),
       nodeParameters=double(maxNodes*6),
       treeStart=integer(nBoosts),treeEnd=integer(nBoosts),totNumNodes=integer(1),
       oobError=double(nBoosts),as.integer(recordBoosts),iPred=double(n*nr),
       as.integer(maxSplitAttempts),
       as.integer(ranSeed),PACKAGE = "gbev")
     
#print(Sys.time());
  
list(totNumNodes=hfit$totNumNodes,oobError=hfit$oobError,treeStart=hfit$treeStart,treeEnd=hfit$treeEnd,
nodeParameters=hfit$nodeParameters,nodePositions=hfit$nodePositions,iPred=hfit$iPred,recordPred=recordBoosts)

}




treeMatrix<-function(object,deci=3)
{
numNodeParam<-4  # number of node-parameters (i.e doubles)
numNodePosit<-5  # number of node-positions  (i.e integers)



if(is.null(object))
{
}

tM<-array(dim=c(object$totNumNodes,(numNodeParam+numNodePosit)))

tM[,1:numNodePosit]<-matrix(object$nodePositions,ncol=object$totNumNodes,nrow=numNodePosit,byrow=F)
tM[,(numNodePosit+1):(numNodeParam+numNodePosit)]<-matrix(round(object$nodeParameters,deci),ncol=object$totNumNodes,nrow=numNodeParam,byrow=F)


for(i in 2:numNodePosit)
tM[tM[,i]!=-99,i]<-tM[tM[,i]!=-99,i]+1

#print(tMList)


colnames(tM)<-c("Node","ParentNode","LeftChild","RightChild","SplitVar","SplitPoint","theta","dRSS","Eobs")



tM[tM==-99]<-NA

tM[1,2]<-0
tM[,2]<-tM[,2]-1
tM[tM[,2]==0,2]<--1
tM[is.na(tM[,3]),3]<--1
tM[is.na(tM[,4]),4]<--1
tM[is.na(tM[,5]),5]<--1
tM[is.na(tM[,8]),8]<-0

tM
}



checkMeasErrorModel<-function(measErrorModel)
{

if(is.null(measErrorModel$numComp))
stop(" Need specify number of mixture components (in 'numComp').")

if(is.null(measErrorModel$mu))
stop(" Need specify means of mixture components (in 'mu').")

if(is.null(measErrorModel$SigmaME))
stop(" Need specify covariance matrix of measurement error (in 'SigmaME').")


if(is.null(measErrorModel$SigmaX))
stop(" Need specify covariance matrix of mixtures (in 'SigmaX').")

if(is.null(measErrorModel$pComp))
stop(" Need specify mixture probabilities (in 'pComp').")


if(!is.null(measErrorModel$beta))
stop(" Dependence of distribution of X on other covariates, not yet implemented.")


}







var.imp<-function(object,nboosts=NULL,std=FALSE,plt=TRUE)
{
if(class(object)!="gbev")
stop("'object' not of class 'gbev' ");

if(is.null(nboosts)){
nboosts<-object$nTrees
}

treeEnd<-object$treeEnd+1
treeEnd<-treeEnd[nboosts]
splitVar<-object$nodePositions[(1+(4*object$totNumNodes)):(5*object$totNumNodes)]+1
splitVar<-splitVar[1:treeEnd]
dRSS<-object$nodeParameters[(1+(2*object$totNumNodes)):(3*object$totNumNodes)]
dRSS<-dRSS[1:treeEnd]

varImp<-array(dim=c(object$ncovar,1))

for(i in 1:object$ncovar){
varImp[i]<-sum(dRSS[splitVar==i])/object$n

}

if(std){

varImp<-(varImp/var(object$y)) 
}

rownames(varImp)<-object$namesCovar

colnames(varImp)<-"Importance"
print(varImp)
#if(std)
#cat("Importance is reducution in MSE due to covariates relative to variance of response.\n");

if(plt)
r <- barplot(t(varImp), col='gray',beside=T,ylab="Importance")

}








predict_gbev<-function(object,newdata,latent=TRUE,firstTree=NULL,lastTree=NULL)
{

###################################################################################
## Can compute multiple predictions for each individual if "length(lastTree)>1"
###################################################################################



if(is.null(lastTree))
{
lastTree<-object$nTrees
}
num_EndTrees<-length(lastTree)



if(is.data.frame(newdata))
{
  
      x <- model.frame(delete.response(object$Terms),
                       newdata,
                       na.action=na.pass)
   }
   else
   {  

     if(ncol(newdata)!=length(attr(object$Terms,"term.labels")))
       stop(" number of columns in 'newdata' doesn't match predictors in model.");      

      x <- newdata
   }

 
if(is.null(firstTree)||is.null(lastTree))
{
firstTree<-1
lastTree<-object$nTrees
}
if(firstTree<1)
{
firstTree<-1
}
for(i in 1:num_EndTrees)
{
if(lastTree[i]>object$nTrees)
{
lastTree[i]<-object$nTrees
warning("'lastTree' exceeded number fit, set to 'object$nTrees'");
}
}
########################################
# This is version where "x" is matrix. 
# Need extend to when x is data frame
#########################################

nPred<-length(x[,1])


p<-length(x[1,])

if(p!=length(object$namesCovar))
stop(" Incorrect dimension for x ");
 

wh<-array(dim=c(nPred*p,1))

for(i in 1:nPred)
wh[(1+(i-1)*p):((i)*p)]<-x[i,];


#h1<-.C("predictEnsembleR",as.double(wh),as.integer(nPred),as.integer(p),as.integer(object$nTrees),as.integer(object$nodePositions),as.double(object$nodeParameters),as.integer(object$treeStart),as.integer(object$treeEnd),as.integer(object$totNumNodes),res=double(nPred))
if(latent)
{


h1<-.C("predictEnsembleR2_mult",as.double(wh),as.integer(nPred),as.integer(p),as.integer(firstTree),
      as.integer(lastTree),as.integer(num_EndTrees),as.integer(object$nodePositions),as.double(object$nodeParameters),
      as.integer(object$treeStart),as.integer(object$treeEnd),as.integer(object$totNumNodes),
      res=double(nPred*num_EndTrees),PACKAGE = "gbev")


res<-matrix(h1$res,ncol=num_EndTrees,nrow=nPred,byrow=FALSE)
} # matches if(latent) 
if(!latent)
{

cat(" Observed data model prediction is not yet implemented.\n");


}  # matches if(!latent) 



list(pred=res)
}



































#############################################################################################################
# FUNCTION:   "mixModelMatrices" 
# 
# ARGUMENTS:  w           (dim=c(n,p)) contains error measured covariates
#             numComp     number of mixture compontents
#             mu          (dim(p,numComp)), mu[,j] mean vector of j-th mixture 
#             SigmaME     Var-Covar of measurement error matrix
#             SigmaJ      SigmaJ[j,,] is variance-covariance matrix of j-th mixture component
#             pComp       pComp[j] is mixture probability of j-th mixture
#
# RETURNS:  list:  "condMean" a vector containing conditional means for all w[i,]'s from each mixture component 
#                  "compProb" compProb[i,j] is "probability of x[i,] coming from j-th component"
#                  "cholCov"  cholesky factorization of the conditional var-covar of j-th component density 
#
# PURPOSE:  This function produces vectors that are passed on to the C-function "treeToR_MC"
#           and are used for simulating MC-data for estimating (using crude MC) node belonging probabilities
#############################################################################################################
mixModelMatrices<-function(w,numComp,mu,SigmaME,SigmaJ,pComp)
{

p<-length(w[1,])
n<-length(w[,1])




### VECTORS TO STOR RESULTS ######
condMean<-array(dim=c(n*numComp*p,1))
compProb<-array(dim=c(n*numComp,1))
cholCov<-array(dim=c(numComp*p*p,1))



Ainv<-array(dim=c(numComp,p,p))
SigmaJinv<-array(dim=c(numComp,p,p))
alphaJ<-array(dim=c(numComp,p))
SigmaMEinv<-solve(SigmaME)


for(i in 1:numComp)
{
SigmaJinv[i,,]<-solve(SigmaJ[i,,])
Ainv[i,,]<-solve(SigmaMEinv+SigmaJinv[i,,])
alphaJ[i,]<-SigmaJinv[i,,]%*%mu[i,]
}


#### WRITE "condMean" ####
for(i in 1:n)
{
 
for(j in 1:numComp)
{
startIndx<-(i-1)*numComp*p+(j-1)*p+1
endIndx<-(i-1)*numComp*p+j*p


condMean[startIndx:endIndx]<-Ainv[j,,]%*%(SigmaMEinv%*%w[i,]+alphaJ[j,])
}
}

#### WRITE "compProb" ####

for(i in 1:n)
{
for(j in 1:numComp)
{
compProb[(i-1)*numComp+j]<-dmvnorm(w[i,], mean=mu[j,], sigma=(SigmaME+SigmaJ[j,,]), log=FALSE)
}
startIndx<-(i-1)*numComp+1
endIndx<-i*numComp
compProb[startIndx:endIndx]<-compProb[startIndx:endIndx]*pComp
hs<-sum(compProb[startIndx:endIndx])
compProb[startIndx:endIndx]<-compProb[startIndx:endIndx]/hs
}

#### WRITE "cholCov" ####
for(j in 1:numComp)
{
hc<-chol(Ainv[j,,])
hc<-t(hc)
for(i in 1:p)
{
startIndx<-(j-1)*p*p+(i-1)*p+1
endIndx<-(j-1)*p*p+i*p
cholCov[startIndx:endIndx]<-hc[i,]
}
}

list(condMean=condMean,cholCov=cholCov,compProb=compProb,numComp=numComp)

}



plotLoss<-function(object,loss="L2",startIter=1,plt=TRUE)
{
if(length(object$recordPred)<=1)
{
stop("No intermediate predictions were saved, nothing to plot.");
}




if(loss=="L2")
{
np<-length(object$recordPred)
ninc<-object$recordPred[2]-object$recordPred[1]
res<-array(dim=c(np,1))

for(i in 1:np)
{
res[i,1]<-mean((object$y-object$iPred[(1+(i-1)*object$n):(i*object$n)])^2)
}

}

if(loss=="logLike")
{

np<-length(object$recordPred)
ninc<-object$recordPred[2]-object$recordPred[1]
res<-array(dim=c(np,1))

for(i in 1:np)
{
#res[i,1]<-mean((object$y-object$iPred[(1+(i-1)*object$n):(i*object$n)])^2)
res[i,1]<-mean(object$y*log(object$iPred[(1+(i-1)*object$n):(i*object$n)]))
res[i,1]<-res[i,1]+mean((1-object$y)*log(1-object$iPred[(1+(i-1)*object$n):(i*object$n)]))
}
res<--res
}

iters<-c(1,(c(2:np)*ninc))
if(plt)
{
xaxis<-iters
hindx<-xaxis>=startIter
plot(xaxis[hindx],res[hindx],type="l",xlab="number of boosts",ylab="loss")
points(xaxis[hindx],res[hindx],type="p")
}
list(iters=iters,loss=res)
}





















# part.dep(object=fit,varIndx=1,firstTree=NULL,lastTree=NULL)
part.dep<-function(object,varIndx,ngrid=50,firstTree=NULL,lastTree=NULL,plt=TRUE)
{





 
if(is.null(firstTree)||is.null(lastTree))
{
firstTree<-1
lastTree<-object$nTrees
}
if(firstTree<1)
{
firstTree<-1
}
if(lastTree>object$nTrees)
{
lastTree<-object$nTrees
warning("'lastTree' exceeded number fit, set to 'object$nTrees'");
}


if(length(varIndx)==1)
{
########################################
# This is version where "x" is matrix. 
# Need extend to when x is data frame
#########################################

min.x<-as.numeric(object$varLevels[[varIndx]])[1]
max.x<-as.numeric(object$varLevels[[varIndx]])[11]

#x<-seq(max(min.x,xlim[1]),min(max.x,xlim[2]),by=((max.x-min.x)/ngrid))
 

x<-seq(min.x,max.x,by=((max.x-min.x)/ngrid))
 

nx<-length(x)
x<-matrix(x,ncol=1,nrow=nx)
nPred<-nx


p<-length(x[1,])



wh<-array(dim=c(nPred*p,1))

for(i in 1:nPred)
wh[(1+(i-1)*p):((i)*p)]<-x[i,];

h1<-.C("partialDependenceR",as.double(wh),as.integer(nPred),as.integer(p),as.integer(object$ncovar),
      as.integer(varIndx-1),as.integer(firstTree),
      as.integer(lastTree),as.integer(object$nodePositions),as.double(object$nodeParameters),
      as.integer(object$treeStart),as.integer(object$treeEnd),as.integer(object$totNumNodes),
      res=double(nPred),PACKAGE = "gbev")


###partialDependenceR(double *x,int *nPred,int *pd,int *p,int *partDepVar,int *startEnsemble,int *endEnsemble,int *nodePositions,double *nodeParameters,int *treeStart,int *treeEnd,int *totNumNodes,double *res);

if(plt)
plot(x,h1$res,ylab="partial dependence",xlab=fit$namesCovar[varIndx],type="l")

dat.bv<-0
}

if(length(varIndx)==2)
{


min.x1<-as.numeric(object$varLevels[[varIndx[1]]])[1]
max.x1<-as.numeric(object$varLevels[[varIndx[1]]])[11]

min.x2<-as.numeric(object$varLevels[[varIndx[2]]])[1]
max.x2<-as.numeric(object$varLevels[[varIndx[2]]])[11]


x1<-seq(min.x1,max.x1,by=((max.x1-min.x1)/ngrid))
x2<-seq(min.x2,max.x2,by=((max.x2-min.x2)/ngrid))



nx1<-length(x1)
nx2<-length(x2)


eval.grid<-array(dim=c(nx1*nx2,2))
for(i in 1:nx1)
{for(j in 1:nx2){eval.grid[(i-1)*nx2+j,1]<-x1[i];eval.grid[(i-1)*nx2+j,2]<-x2[j];}}

x<-eval.grid
nPred<-length(x[,1])

p<-length(x[1,])

vI<-varIndx-1

wh<-array(dim=c(nPred*p,1))

for(i in 1:nPred)
wh[(1+(i-1)*p):((i)*p)]<-x[i,];


h1<-.C("partialDependenceR",as.double(wh),as.integer(nPred),as.integer(p),as.integer(object$ncovar),
      as.integer(vI),as.integer(firstTree),
      as.integer(lastTree),as.integer(object$nodePositions),as.double(object$nodeParameters),
      as.integer(object$treeStart),as.integer(object$treeEnd),as.integer(object$totNumNodes),
      res=double(nPred),PACKAGE = "gbev")


dat.bv<-data.frame(x=eval.grid[,1],y=eval.grid[,2],z=h1$res)
#wireframe(z~x*y,dat.bv,aspect=c(1,0.5),drape=T,screen=list(z=-150,x=-60),colorkey=list(space="right",height=0.6))
 contour(x=x1,y=x2,
        z=matrix(dat.bv$z,ncol=nx2,nrow=nx1,byrow=TRUE),
        xlab=object$namesCovar[varIndx[1]],ylab=object$namesCovar[varIndx[2]])

}

if(length(varIndx)>2)
{
stop(" partial dependence plots are implemented only for 1 and 2 variables, not more. ");
}


list(pred=h1$res,x=x,dat=dat.bv)


}

cvLoss<-function(object,k,random=F,loss="logLike")
{


supported.methods <-c("logLike","L2")  
   if(!is.element(loss,supported.methods))
   {
      stop("loss ",loss," is not supported. Choose logLike or L2.")
   }

nP<-length(object$recordPred)

res<-array(dim=c(nP,k))
nk<-round(object$n/k)

cvWeightsMatrix<-array(1,dim=c(object$n,k))

if(random==FALSE)
{
for(i in 1:k)
{
if(i<k)
cvWeightsMatrix[(1+(i-1)*nk):(nk*i),i]<-0

if(i==k)
cvWeightsMatrix[(1+(i-1)*nk):(object$n),i]<-0
}
}

if(random==TRUE)
{

for(i in 1:k)
{
indx<-c(1:object$n)


if(i==1){

setToZero<-sample(indx,size=nk,replace=F)
takeOut<-setToZero
cvWeightsMatrix[setToZero,i]<-0

}

if((i>1)&(i<k)){
indx<-indx[-takeOut]
setToZero<-sample(indx,size=nk,replace=F)
takeOut<-c(takeOut,setToZero)
cvWeightsMatrix[setToZero,i]<-0
}



if(i==k){
indx<-indx[-takeOut]
cvWeightsMatrix[indx,i]<-0
}

}
}



for(i in 1:k)
{
#cat(" Cross-validation iteration: ",i," of ",k," \n");
cvWeights<-cvWeightsMatrix[,i]



#print(cvWeights)

pfit<-gbev.fit(w=object$w,y=object$y,
               weights=cvWeights,
               measErrorModel=object$measErrorModel,
               method=object$method,
               indepFitUpdate=object$indepFitUpdate,
               nboost=object$nboost,
              lambda=(object$lambda*sum(cvWeights)/length(object$y)),
               maxDepth=object$maxDepth,
               m=object$ms, 
               minSplit=object$minSplit,
               minBucket=object$minBucket,
               sPoints=object$sPoints,
               mc=object$mc,
               intermPred=object$intermPred,
               maxSplitAttempts=object$maxSplitAttempts)





for(j in 1:nP)
{

hind<-(1+(j-1)*object$n):(j*object$n)

if(loss=="logLike")
{
hhh1<-sum(object$y*log(pfit$iPred[hind])*as.numeric(cvWeights==0))
hhh1<-hhh1/sum(cvWeights==0)
res[j,i]<-hhh1
hhh1<-sum((1-object$y)*log(1-pfit$iPred[hind])*as.numeric(cvWeights==0))
hhh1<-hhh1/sum(cvWeights==0)
res[j,i]<-res[j,i]+hhh1
}
if(loss=="L2")
{
hhh2<-(object$y-pfit$iPred[hind])*as.numeric(cvWeights==0)
hhh1<-sum(hhh2^2)
hhh1<-hhh1/sum(cvWeights==0)
res[j,i]<-hhh1
}


}



}
resFinal<-apply(res,1,mean)
#plot(c(1:nP)*object$intermPred,resFinal,type="l")

iters<-c(1:nP)*object$intermPred
if(loss=="logLike")
{
mL<-max(resFinal)
estIter<-iters[resFinal==mL][1]
resFinal<--resFinal
}
if(loss=="L2")
{
mL<-min(resFinal)
estIter<-iters[resFinal==mL][1]
}
list(iters=iters,cvLoss=resFinal,estIter=estIter)

}




















