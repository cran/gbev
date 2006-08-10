\name{gbev}
\alias{gbev}
\title{Boosted regression trees with errors-in-variables}
\description{Fits boosted regression trees with errors-in-variables.
}
\usage{

gbev(formula = formula(data),
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
}
\arguments{
\item{formula}{A symbolic description of the model
to be fit. NOTE: the error measured covariates must be placed first in the \code{formula}.}
\item{data}{A data frame containing variables in the model.}
\item{weights}{Weights applied to observations. Defaults to 1.}
\item{measErrorModel}{This is a list specifying the distribution of the latent 
covariates and the measurement error. Here it is assumed that the latent covariates 
are a mixture of normals (possibly multivariate), and that the measurement error is 
normally distributed. See examples below for details. Note that at least one 
covariate must be specified as measured with error.}
\item{method}{ Can be \code{L2} for squared error loss, \code{logLike} for binary regression with negative log-likelihood loss.}
\item{indepFitUpdate}{ If \code{indepFitUpdate=1} then model fit is updated using independent MC-sample,
                          else the fit is updated using the same MC-sample as used in tree fitting.}

\item{nboost}{ Number of boosting iterations to perform. }
\item{lambda}{ Regularization parameter. }
\item{maxDepth}{Determines maximum interaction depth of trees fit. \code{maxDepth=2} fits stumps.}
\item{minSplit}{Minimum expected number of observations in a node for it to be split. }
\item{minBucket}{Minimum expected number of observations in a node. }
\item{sPoints}{Number of points sampled at random from which to choose split. }
\item{mc}{Number of Monte-Carlo samples drawn to compute node probabilities. }
\item{intermPred}{Increments of iterations at which intermediate predictions are saved, 
required for \code{cvLoss} function.}
\item{maxSplitAttempts}{Maximum number of attempts at finding a valid split point. When 
splitting a node, \code{sPoints} candidate splits are supposed 
to be found for each covariate, however, each randomly sampled split point 
does not necessarily give 
a valid split point (i.e. a point satisfying \code{minBucket} and \code{minSplit}), and 
\code{maxSplitAttempts} is the maximum number of attempts at finding such a point.   
.}

}
\details{This function performs non-parametric regression when some or all 
covariates are measured 
with error. It is assumed that the latent (error measured) covariate vector, \eqn{X}, 
is observed via a random variable \eqn{W}, the relation between the two being
\deqn{W=X+U}
where \eqn{U} is the measurement error. This function assumes that the density of \eqn{X} 
is known and representable as a finite mixture of multivariate normal densities, while the measurement 
error \eqn{U} is assumed to be multivariate normal. These two densities are specified 
in the \code{measErrorModel} argument of the \code{gbev} function, see examples 
below for further details on  \code{measErrorModel}. In practice the densities of \eqn{X} 
and \eqn{U} will not be known, but must be estimated in some manner. This is typically 
done using replicate measurements or validation data, see Carroll et.al. (1995) for further 
details. 

It is assumed that in the \code{formula} argument the error measured covariates occur 
first. That is, for example, if one has covariates \eqn{W1}, \eqn{W2}, \eqn{Z1} and 
\eqn{Z2}, where the first two are error contaminated and the last two error-free then 
\code{formula} must be \code{y~w1+w2+z1+z2}.

The function \code{gbev} estimates \eqn{E(Y|X)} where \eqn{Y} is the response, 
using boosted regression trees. That is, a model \eqn{F(x)} for \eqn{E(Y|x)} is estimated 
such that  \eqn{F(x)=\sum_{k=1}^{K}T_k(x)} where  each \eqn{T_k(x)} is a regression tree, 
and each regression tree is estimated by fitting it to the gradient of the previous model 
estimate using least squares, as described in Friedman (2001). 

The model \eqn{F(x)} is built up iteratively, such that at each boosting iteration 
a tree model is added to \eqn{F(x)} to minimize the squared error between the 
observed data model, given by \eqn{E(F(X)|w)}, and the current model residuals. More 
specifically, on the \eqn{k}-th iteration, the observed data residuals are 
\eqn{\tilde{y_i}=(y_i-E(F_k(X)|w))} for \eqn{i=1,..,n}, and a tree 
\eqn{T_{k}(x)=\sum_{j=1}^{J}c_{kj}I(x\in R_{kj})} is fit by fitting \eqn{E(T_{k}(X)|w)}, 
using least squares and recursive partitioning, to the current residuals \eqn{\tilde{y_i}}. 
The model estimate is then 
updated by setting \eqn{F_{k+1}(x)=F_{k}(x)+T_{k}(x)}, and the procedure iterated.

 
The tree fitting requires evaluating 
probabilities of the form \eqn{Pr(X\in R|w)=E(I(X\in R)|w)}, i.e. the probability that the covariate 
\eqn{X} is in a rectangular region \eqn{R}, given that \eqn{W=w}. These probabilities 
are here estimated using Monte Carlo sampling, and the argument \code{mc} in \code{gbev} 
regulates how many samples are drawn, at each iteration, for each observation, 
the default being \code{mc=2}. 
Boosting is known to perform best when the amount of fitting at each iteration is small, 
here this is regulated by the size of the regression trees, the \code{maxDepth} argument, 
and the regularization parameter \code{lambda}. The larger values of \code{lambda} the 
less fitting is done per boosting iteration, and the more iterations are typically 
required to achieve adequate fitting. The \code{lambda} parameter is also used to control 
the Monte Carlo error in the function estimate, due to the sampling in the tree fitting. 
It turns out that the larger value of \code{lambda} the smaller Monte Carlo error is in 
the regression function estimate, for the same value of \code{mc}. Some experimentation 
is required to find an appropriate \code{lambda} value such that the Monte Carlo 
error is acceptable.    

Two loss functions are available for estimation, specified in the \code{method} argument. 
Setting \code{method="L2"} squared error loss is used, while \code{method="logLike"} 
negative log-likelihood loss is used for performing binary, \eqn{Y} is \eqn{0} or \eqn{1}, 
regression. 

The arguments \code{maxDepth}, \code{minSplit}, \code{minBucket}, and 
\code{sPoints} control the tree fitting. \code{maxDepth} controls the depth of the 
regression tree to be fit, with \code{maxDepth=2} fitting a tree containing a single 
split (two terminal nodes), \code{maxDepth=3} fits a tree with 4 terminal nodes, obtained 
by splitting the nodes of a single split tree. Interestingly, \code{maxDepth=2} fits 
an additive model, \code{maxDepth=3} a model with all second order interactions, 
\code{maxDepth=4} a model with all third order interactions, and so on. The argument 
\code{sPoints} is the number of candidate split points sampled for each covariate when 
splitting a node. In the absense of measurement error node splitting is typically done 
by examining all values of the observed covariates falling in the node being split, however 
with measurement error the realized values of \eqn{X} are not observed, and sampling 
is here used to generate candidate split points, with \code{sPoints} governing the number 
of such split points sampled for each covariate. Interestingly, smaller values of 
\code{sPoints} often work the best. The argument \code{minSplit} determines the 
smallest (expected) number of observations in a node for a split to be attempted, 
while the argument 
\code{minBucket} the smallest (expected) number of observations in a terminal node. Note 
that using non-zero values of \code{lambda}, \code{minBucket=0} will run without error. 
  


}
\value{
\code{gbev} returns
\code{gbev.object}.
}
\references{

J.H. Friedman (2001). "Greedy Function Approximation: A Gradient Boosting
Machine," Annals of Statistics 29(5):1189-1232.

T. Hastie, R. Tibshirani and J.H. Friedman (2001). "The Elements of Statistical Learning" 
Springer. 

R. Carroll, D. Ruppert and L. Stenfanski (1995). "Measurement Error in Nonlinear Models," 
Chapman and Hall.


}

\author{Joe Sexton \email{j.a.sexton@medisin.uio.no}}
\examples{

### Univariate regression example
n<-500
varX<-1
varME<-0.25
varNoise<-0.3^2

### Data 
x<-rnorm(n,sd=sqrt(varX))                              ### Error free covariate
w<-x+rnorm(n,sd=sqrt(varME))                           ### Error contaminated version
fx<-sin(pi*x/2)/(1+2*(x^2)*((2*as.numeric(x>=0)-1)+1)) ### True regression function  
y<-fx+rnorm(n,sd=sqrt(varNoise))                       ### Response                           
dat<-data.frame(y=y,w=w)

### Measurement error model ####
###  
### The measurement error model is a list of the following components:
###
### SigmaX:    the covariance matrices of the mixture model for the error free covariates 
###            SigmaX[i,,] is the covariance matrix of the i-th mixture density
### mu:        the means of the mixture model for the error free covariates 
###            mu[i,] is the mean-vector of the i-th mixture density
### SigmaME:   the covariance matrix of the measurment error
### pComp:     the weights of the mixture distribution, pComp[i] is the weight of the 
###            i-th mixture density
### numComp:   the number of components in the mixture 
##
p<-1
pME<-1

numComp<-3                                    ## number of components in gaussian mixture for X-distribution
SigmaME<-diag(varME,pME)
SigmaJ<-array(dim=c(numComp,pME,pME))         
mu<-array(dim=c(numComp,pME))
pComp<-array(1/numComp,dim=c(numComp,1))
for(i in 1:numComp)
{
SigmaJ[i,,]<-diag(varX,pME)
mu[i,]<-rep(0,pME)
}
### list required by "gbev" for measurement error model
meModel<-list(SigmaX=SigmaJ,mu=mu,SigmaME=SigmaME,pComp=pComp,numComp=numComp)


fit<-gbev(y~w,data=dat,
          measErrorModel=meModel,     
          method="L2",              ## Squared error loss
          nboost=1000,              ## 1000 boosting iterations
          lambda=5,                 ## regularization of regression tree
          maxDepth=2,               ## maximum tree depth, 2 corresponds stumps
          mc=2,                     ## number of monte-carlo samples per tree build 
          minSplit=3,               ## minimum number of obs in node to split
          minBucket=0,              ## minimum number of obs in nodes
          sPoints=10,               ## number of sampled candidate split points
          intermPred=5)             ## increments of iterations to store predictions 

### 5-fold cross-validation
hcv<-cvLoss(object=fit,k=5,random=FALSE,loss="L2")
plot(hcv$iters,hcv$cvLoss,type="l")

hp<-part.dep(object=fit,varIndx=1,firstTree=1,lastTree=hcv$estIter)

x<-seq(-2,2,by=.02)
fx<-sin(pi*x/2)/(1+2*(x^2)*((2*as.numeric(x>=0)-1)+1)) 
points(x,fx,type="l",lty=5)




## Simulated binary regression example, 
## with: Y=I( X1*X2+X2*X3+X1*X3>0), with measurement error on X's
n<-1000
p<-3
varX<-1     ## 
varME<-0.5  ## measurement error variance

x<-rnorm(p*n)
x<-matrix(x,ncol=p,nrow=n)
## add measurement error
w<-x+matrix(rnorm(p*n,sd=sqrt(varME)),ncol=p,nrow=n)   

x<-x[,c(1:p)]*x[,c(2:p,1)]
x<-apply(x,1,sum)
threshold<-0
y<-as.numeric(x>threshold)
dat<-data.frame(y=y,w1=w[,1],w2=w[,2],w3=w[,3])  ##  must be modified if(p!=3)


#### Measurement error model ######
numComp<-1                               ##  Number of components in mixture 
SigmaME<-diag(varME,p)                   ##  Covariance matrix of measurement error
SigmaJ<-array(dim=c(numComp,p,p))        ##  Covariance matices for mixture
mu<-array(dim=c(numComp,p))              ##  Mean vectors for mixture components
pComp<-array(1/numComp,dim=c(numComp,1)) ##  Mixture probabilities
for(i in 1:numComp)
{                                        ## filling in mixture model for X-distribution
SigmaJ[i,,]<-diag(varX,p)
mu[i,]<-rep(0,p)
}
## The list for measurement error model 
meModel<-list(SigmaX=SigmaJ,mu=mu,SigmaME=SigmaME,pComp=pComp,numComp=numComp)

fit<-gbev(y~w1+w2+w3,data=dat,
         measErrorModel=meModel,   
         method="logLike",       ## loss function
         nboost=1000,            ## number of boosting iterations
         lambda=40,              ## regularization parameter used in regression tree
         maxDepth=3,             ## maximum depth of regression tree 
         minSplit=10,             ## minimum number of observations in  node to  split
         minBucket=0,            ## minimum number in split node to allow split
         sPoints=2,             ## number of sampled canditate split points 
         mc=2,                   ## monte-carlo sample size used in each regression tree
         intermPred=10)          ## Increments of iterations to store loss function
        

## plot loss function as function of iterations
hp<-plotLoss(fit,loss="logLike",startIter=10)

## bivariate partial dependence plot
hdp<-part.dep(object=fit,varIndx=c(1,2),firstTree=1,
lastTree=1000,ngrid=50)
dpp<-data.frame(x1=hdp$dat$x,x2=hdp$dat$y,prob=hdp$dat$z)
library(lattice)
wireframe(prob~x1*x2,dpp,aspect=c(1,0.5),drape=TRUE,screen=list(z=50,x=-60),
scales=list(arrows=FALSE),xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))

}
\keyword{nonparametric}
\keyword{tree}
