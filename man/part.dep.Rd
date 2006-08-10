\name{part.dep}
\alias{part.dep}
\title{ Partial dependence}
\description{ Partial dependence function for boosted tree model.
}
\usage{
part.dep(object,varIndx,ngrid=50,
firstTree=NULL,lastTree=NULL,plt=TRUE)
}
\arguments{
\item{object}{A fitted object of type \code{gbev.object}.}
\item{varIndx}{Index of variables to compute partial dependence, 
               can be of length 1 or 2, corresponding univariate and bivariate 
               partial dependence functions.}
\item{ngrid}{Number of grid points along each axis at which to compute partial dependence function.}
\item{firstTree}{Index of tree in boosted tree sequence to start computations, defaults to 1.}
\item{lastTree}{Index of last tree in boosted tree sequence to use for computations, defaults to number of boosting iterations.}
\item{plt}{If \code{TRUE} then the partial dependence is plotted.}
}
\details{   
Computes partial dependence functions used to summarized the marginal effect of 
covariates on the response. These were introduced in Friedman (2001) and are also 
described in Hastie et.al. (2001) in chapter 10.13.2. The \code{part.dep} function, however, 
computes the partial dependence of the latent covariates, the \eqn{X}'s, on the response, 
and not the error-contaminated \eqn{W}'s. This requires a slight modification of the 
procedure described in the references. Specifically, Friedman (2001) describes a 
procedure where to compute the partial dependence using boosted tree models one 
needs to know the proportion of the observations falling in the various terminal nodes 
of the trees. With measurement error, however, this proportion is not observable and 
must be estimated, which here is done using the Monte Carlo samples used in the tree 
fitting of each boosting iteration. 

}
\value{ A list with elements \code{pred}, \code{x} and \code{dat}.  
        If univariate partial dependence, \code{pred} is the partial dependence at the 
        points \code{x}. If bivariate partial dependence,  
        then \code{dat} is a data-frame with the
        partial dependence given in variable \code{z} evaluated at the variables 
        \code{x} and \code{y}.}
\references{ 
J.H. Friedman (2001). "Greedy Function Approximation: A Gradient Boosting
Machine," Annals of Statistics 29(5):1189-1232.

T. Hastie, R. Tibshirani and J.H. Friedman (2001). "The Elements of Statistical Learning" 
Springer. 

}
\keyword{ hplot }

