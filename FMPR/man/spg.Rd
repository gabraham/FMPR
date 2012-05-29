\name{spg}
\alias{spg}
\title{Fit an L1-penalised linear model with an L1 fusion penalty (graph fused
lasso), using the smoothed proximal gradient method}
\description{
   Fit an L1 linear model.
  }
\usage{
spg(X, Y, C=NULL, lambda=0, gamma=0, C=NULL, tol=1e-6,
    mu=1e-04, maxiter=10000, simplify=FALSE, verbose=FALSE) 
}

\arguments{
  \item{X}{input matrix, of dimension nobs x nvars.}
  \item{Y}{response variable.}
  \item{lambda}{vector of non-negative penalties}
  \item{gamma}{vector of non-negative fusion penalties}
  \item{C}{K by K matrix}
  \item{tol}{tolerance for SPG method}
  \item{mu}{solution accuracy for SPG method}
  \item{maxiter}{maximum number of iterations}
  \item{simplify}{logical, if the returned list contains one matrix, return it instead
  of a list}
  \item{verbose}{logical, print verbose information}
}
\details{
}
\value{
If simplify is FALSE, or any of the number of penalties is more than one,
a nested list, one level for each of lambda, and gamma, with the final
level being the estimated weights B. Otherwise, a matrix.
}
\references{
Adapted from MATLAB code by Xi Chen,
\url{http://www.cs.cmu.edu/~xichen/Code/SPG_Multi_Graph.zip}
}
\author{Gad Abraham\cr
Maintainer: Gad Abraham \email{gad.abraham@unimelb.edu.au}}
\seealso{
\code{optim.spg}, \code{crossval.spg},
\code{gennetwork}
}
\examples{
N <- 100 # samples
p <- 50  # inputs (variables, predictors)
K <- 10  # tasks
X <- scale(matrix(rnorm(N * p), N, p))
B <- matrix(rnorm(p * K), p, K) # dense weights
Y <- scale(X \%*\% B + rnorm(N * K))
C <- gennetwork(Y)
l <- max(maxlambda1(X, Y)) * 2^seq(-3, 0, length=10)
g <- 2^seq(-3, 3, length=10)
s <- spg(X=X, Y=Y, C=C, lambda=l, gamma=g)
B <- s[[3]][[1]] # weights for lambda[3], gamma[1]
}
\keyword{models}
\keyword{regression}

 
