\name{mixModelMatrices}
\alias{mixModelMatrices}
\title{Helper function for gbev-function.}
\description{
A help function for creating measurement error model, not intended for end-user use.
}
\usage{
mixModelMatrices(w,numComp,mu,SigmaME,SigmaJ,pComp)
}
\arguments{
\item{w}{A matrix of covariates.}
\item{numComp}{Number of components in mixture of normals model.}
\item{mu}{Mean vectors of mixture of normal model.}
\item{SigmaME}{Covariance matrix of measurement error.}
\item{SigmaJ}{Covariance matrices of mixture of normal model.}
\item{pComp}{Vector of weights of mixture.}

}
\keyword{hplot}

