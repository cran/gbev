\name{var.imp}
\alias{var.imp}
\title{Variable importance.}
\description{
Computes and plots the variable importance of variables used in model. 
}
\usage{
var.imp(object,nboosts=NULL,std=FALSE,plt=TRUE)
}
\arguments{
\item{object}{A fitted object of type \code{gbev.object}.}
\item{nboosts}{Number of boosting iterations, if \code{NULL} set to maximum number used.}
\item{std}{If FALSE then the reduction in mean squared error due to covariates is given, 
if TRUE then this reduction is standardized with respect to variance of response.}
\item{plt}{If TRUE the importances are plotted.}
}
\details{   
Variable importance is defined as the reduction in mean squared error due to the different 
covariates. For a 
given covariate, this is computed by idenifying all nodes split using that covariate 
and adding together the reductions in mean squared error due to splitting these nodes. 
}
\keyword{ hplot }
