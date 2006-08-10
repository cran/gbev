\name{gbev.fit}
\alias{gbev.fit}
\title{Helper function for gbev-function.}
\description{Help function for fitting boosted regression trees with errors-in-variables.}
\usage{

gbev.fit(w,y,
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
}
\arguments{
\item{w}{Matrix of covariates.}
\item{y}{A vector of responses.}
\item{weights}{Weights applied to observations. Defaults to 1.}
\item{measErrorModel}{This is a list specifying the distribution of the latent 
covariates and the measurement error. Here it is assumed that the latent covariates 
are a mixture of normals (possibly multivariate), and that the measurement error is 
normally distributed. See examples below for details.}
\item{method}{ Can be \code{L2} for squared error loss, \code{logLike} for binary regression with negative log-likelihood loss.}
\item{indepFitUpdate}{ If \code{indepFitUpdate=1} then model fit is updated using independent MC-sample,
                          else the fit is updated using the same MC-sample as used in tree fitting.}

\item{nboost}{ Number of boosting iterations to perform. }
\item{lambda}{ Regularization parameter. }
\item{maxDepth}{Determines maximum interaction depth of trees fit. \code{maxDepth=2} fits stumps.}
\item{m}{Number of randomly sampled covariates used to split node, usually set to 
number of covariates.}
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
\details{This function is called by \code{gbev} but is not intended for 
end-user use.
}
\keyword{nonparametric}
\keyword{tree}
