# Test groupridge

library(MASS)

dyn.load("groupridge.so")

set.seed(39827431)

N <- 10
p <- 10
K <- 5
grp <- 1:K

groupridge <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g,
      maxiter=1e6, eps=1e-6)
{
   p <- ncol(X)
   K <- ncol(Y)

   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   if(length(lambda2) == 1)
      lambda2 <- rep(lambda2, K)

   if(length(lambda3) == 1)
      lambda3 <- rep(lambda3, K)

   r <- .C("groupridge", as.numeric(X), as.numeric(Y), 
      numeric(p * K), nrow(X), ncol(X), ncol(Y),
      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
      as.integer(g), as.integer(maxiter), as.double(eps))
   matrix(r[[3]], p, K)
}

X <- matrix(rnorm(N * p), N, p)
Y <- matrix(rnorm(N * K), N, K)

g <- groupridge(X, Y, lambda1=0, lambda2=0, lambda3=0, g=grp, eps=1e-10)
(err <- mean((X %*% g - Y)^2))

b <- ginv(X) %*% Y
(err <- mean((X %*% b - Y)^2))

diag(cor(b, g))

