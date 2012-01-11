library(glmnet)
library(gdata)

dyn.load("~/Code/L1L2/groupridge.so")

# lambda1: scalar or K-vector
# lambda2: scalar or K-vector
# lambda3: scalar 
groupridge3 <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, grp=NULL,
      maxiter=1e5, eps=1e-6, verbose=FALSE)
{
   p <- ncol(X)
   Y <- cbind(Y)
   K <- ncol(Y)

   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   if(length(lambda2) == 1)
      lambda2 <- rep(lambda2, K)

   if(is.null(grp))
      grp <- 1:K

   r <- .C("groupridge3", as.numeric(X), as.numeric(Y), 
      numeric(p * K), nrow(X), ncol(X), K,
      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
      as.integer(grp), as.integer(maxiter),
      as.double(eps), as.integer(verbose), integer(1))

   status <- r[[14]]
   if(!status)
      warning("groupridge failed to converge")
   matrix(r[[3]], p, K)
}

# lambda1: scalar or K-vector
# lambda2: scalar or K-vector
# lambda3: scalar 
groupridge4 <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, G=NULL,
      maxiter=1e5, eps=1e-6, verbose=FALSE)
{
   p <- ncol(X)
   Y <- cbind(Y)
   K <- ncol(Y)

   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   if(length(lambda2) == 1)
      lambda2 <- rep(lambda2, K)

   if(is.null(G))
      G <- matrix(0, K, K)

   r <- .C("groupridge4", as.numeric(X), as.numeric(Y), 
      numeric(p * K), nrow(X), ncol(X), K,
      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
      as.integer(G), as.integer(maxiter),
      as.double(eps), as.integer(verbose), integer(1))
   status <- r[[14]]
   if(!status)
      warning("groupridge failed to converge")
   matrix(r[[3]], p, K)
}

lasso3 <- function(X, y, lambda1=0,
      maxiter=1e5, eps=1e-6, verbose=FALSE)
{
   p <- ncol(X)

   r <- .C("lasso3", as.numeric(X), as.numeric(y), numeric(p), 
      nrow(X), ncol(X), as.numeric(lambda1),
      as.integer(maxiter), as.double(eps),
      as.integer(verbose)
   )
   matrix(r[[3]], p)
}

ridge <- function(X, Y, lambda2=0)
{
   XX <- crossprod(X)
   XY <- crossprod(X, Y)
   p <- ncol(X)
   qr.solve(XX + diag(p) * lambda2, XY)
}

# zero mean, unit norm (not unit variance)
standardise <- function(x)
{
   x1 <- sweep(x, 2, colMeans(x))
   v <- apply(x1, 2, function(z) sqrt(sum(z^2)))
   sweep(x1, 2, v, FUN="/")
}

#source("lasso.R")

maxlambda1 <- function(X, Y)
{
   Y <- cbind(Y)
   r <- .C("maxlambda1", as.numeric(X), as.numeric(Y),
      numeric(ncol(Y)), nrow(X), ncol(X), ncol(Y))
   r[[3]]
}

blockX <- function(X, p, C)
{
   N <- nrow(X)
   
   Xblock <- matrix(0, N * C, p * C) 
   s1 <- seq(1, N * C, N)
   s2 <- seq(1, p * C, p)

   for(i in 1:C)
   {
      j <- s1[i]:(s1[i] + N - 1)
      k <- s2[i]:(s2[i] + p - 1)
      Xblock[j, k] <- X
   }

   Xblock
}

## Threshold the correlation matrix and convert into groups, finding the graph
## that spans the connected vertices. Each connected subgraph is a group.
#cor2grp <- function(Y, thresh=0.7)
#{
#   R <- abs(cor(Y))
#   diag(R) <- 0
#   w <- which(R >= thresh, arr.ind=TRUE)
#   K <- numeric(length(unique(w)))
#}

