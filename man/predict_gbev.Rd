\name{predict_gbev}
\alias{predict_gbev}
\title{Predictions from boosted tree model.}
\description{
Returns value of fitted function at given points.
}
\usage{
predict_gbev(object,newdata,latent=TRUE,firstTree=NULL,lastTree=NULL)
}
\arguments{
\item{object}{A fitted object of type \code{gbev.object}.}
\item{newdata}{A data frame from which predictions are made.}
\item{latent}{If \code{TRUE} the value of the fitted function between the error-free
(but unobserved) covariates and the response is returned, if \code{FALSE} 
nothing happens.}
\item{firstTree}{Index of first tree from which predictions are computed. Defaults to 1.}
\item{lastTree}{Index of last tree from which predictions are computed. Defaults to number 
of trees in fitting.}

}
\value{
\code{pred} the predictions/function evaluations.
}
\keyword{ hplot }

