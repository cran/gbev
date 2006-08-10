\name{cvLoss}
\alias{cvLoss}
\title{Cross-validation of boosting iteration.}
\description{Performs k-fold cross-validation of number of boosting iterations.}
\usage{
cvLoss(object,k,random=F,loss="logLike")
}
\arguments{
\item{object}{A fitted model of type \code{gbev.object}.}
\item{k}{Is the \code{k} in \code{k}-fold cross-validation.}
\item{random}{If \code{TRUE} then cross-validation is done by randomly sampling (without replacement) 
 the validation group, else the validation groups are determined from order 
appearing in data.}
\item{loss}{Can be \code{logLike} for binary regression (negative log-likelihood loss), or \code{L2} for squared error loss.}
}
\value{
Returns a list containing, 
\code{iters} which are iterations the cross-validation loss is evaluated,  
\code{cvLoss} which is the cross-validated loss-function, 
\code{estIter} which is the iteration minimizing the cross-validated loss.
}
\keyword{hplot}

