\name{fmpr}
\alias{fmpr}
\title{Fit an L1-penalised linear model with an L2 fusion penalty}
\description{
   Fit an L1 linear model.
  }
\usage{
fmpr(X, Y, lambda1=0, lambda2=0, lambda3=0, G=NULL,
      maxiter=1e5, eps=1e-6, verbose=FALSE, simplify=FALSE,
      sparse=TRUE, nzmax=nrow(X))
}

\arguments{
  \item{X}{input matrix, of dimension nobs x nvars.}
  \item{Y}{response variable.}
  \item{lambda1}{vector of non-negative penalties}
  \item{lambda2}{vector of non-negative penalties}
  \item{lambda3}{vector of non-negative penalties}
  \item{G}{K by K matrix}
  \item{maxiter}{maximum number of coordinate descent iterations}
  \item{eps}{non-negative tolerance for convergence}
  \item{verbose}{logical, print verbose information}
  \item{simplify}{logical, if the returned list contains one matrix, return it instead
  of a list}
  \item{sparse}{logical, return a sparseMatrix object or not}
  \item{nzmax}{maximum number of variables to be allowed in the model,
  overrides the L1 penalty.}
}
\details{
}
\value{
If simplify is FALSE, or any of the number of penalties is more than one,
a nested list, one level for each of lambda1, lambda2, lambda3, with the final
level being the estimated weights B. Otherwise, a matrix.
}
\references{
}
\author{Gad Abraham\cr
Maintainer: Gad Abraham \email{gad.abraham@unimelb.edu.au}}
\seealso{
\code{optim.fmpr}, \code{crossval.fmpr},
\code{graph.sqr}, \code{graph.abs}
}
\examples{
N <- 100 # samples
p <- 50  # inputs (variables, predictors)
K <- 10  # tasks
X <- scale(matrix(rnorm(N * p), N, p))
B <- matrix(rnorm(p * K), p, K) # dense weights
Y <- scale(X \%*\% B + rnorm(N * K))
G <- graph.sqr(cor(Y))
l1 <- max(maxlambda1(X, Y)) * 2^seq(-3, 0, length=10)
l3 <- 2^seq(-3, 3, length=10)
f <- fmpr(X=X, Y=Y, G=G, lambda1=l1, lambda3=l3)
B <- f[[2]][[1]][[3]] # weights for lambda1[2], lambda2[1], lambda3[3]
}
\keyword{models}
\keyword{regression}

 
