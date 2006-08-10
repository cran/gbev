\name{plotLoss}
\alias{plotLoss}
\title{Plot loss function versus iteration.}
\description{
The loss function at various boosting iterations is plotted. 
}
\usage{
plotLoss(object,loss="L2",startIter=1,plt=TRUE)
}
\arguments{
\item{object}{A fitted object of type \code{gbev.object}.}
\item{loss}{Can be \code{"L2"} or \code{"logLike"} for binary regression.}
\item{startIter}{Iteration at which to start plotting.}
\item{plt}{If \code{TRUE} then the loss is plotted, otherwise it is only returned.}

}
\details{   
The loss function at increments of iterations equal to \code{object$intermPred} is plotted 
and returned.

}
\value{
The function returns 
\code{iters} which are the iterations at which loss is evaluated, and \code{loss} the 
loss function value.  
}

\keyword{ hplot }
