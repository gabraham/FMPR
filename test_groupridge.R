# Test groupridge

library(MASS)
library(glmnet)

dyn.load("groupridge.so")

options(error=dump.frames)

#set.seed(3827431)


groupridge <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g,
      maxiter=1e5, eps=1e-6, verbose=FALSE)
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
      as.integer(g), as.integer(maxiter), as.double(eps),
      as.integer(verbose)
   )
   matrix(r[[3]], p, K)
}

groupridge_simple <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g,
      maxiter=1e3, eps=1e-6, verbose=FALSE)
{
   p <- ncol(X)
   K <- ncol(Y)

   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   if(length(lambda2) == 1)
      lambda2 <- rep(lambda2, K)

   if(length(lambda3) == 1)
      lambda3 <- rep(lambda3, K)

   r <- .C("groupridge_simple", as.numeric(X), as.numeric(Y), 
      numeric(p * K), nrow(X), ncol(X), ncol(Y),
      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
      as.integer(g), as.integer(maxiter), as.double(eps),
      as.integer(verbose)
   )
   matrix(r[[3]], p, K)
}


groupridge2 <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g,
      maxiter=1e5, eps=1e-6, verbose=FALSE)
{
   p <- ncol(X)
   K <- ncol(Y)

   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   if(length(lambda2) == 1)
      lambda2 <- rep(lambda2, K)

   if(length(lambda3) == 1)
      lambda3 <- rep(lambda3, K)

   r <- .C("groupridge2", as.numeric(X), as.numeric(Y), 
      numeric(p * K), nrow(X), ncol(X), ncol(Y),
      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
      as.integer(g), as.integer(maxiter), as.double(eps), as.integer(verbose))
   matrix(r[[3]], p, K)
}

# a very simple version of lasso, for sanity checking
simpler <- function(x, y, lambda, maxiter=200)
{
   N <- length(y)
   p <- ncol(x)
   beta <- numeric(p)
   v <- diag(crossprod(x))
   for(iter in 1:maxiter)
   {
      for(j in 1:p)
      {
	 d1 <- crossprod(x[,j], x %*% beta - y)
	 d2 <- v[j]
	 d <- beta[j] - d1 / d2 - lambda * sign(beta[j])

	 # The solution is zero when:
	 # 1. beta[j] starts as zero, and Newton step is 
	 #    less than lambda
	 # or
	 # 2. beta[j] doesn't start as zero, but changes sign
	 if(beta[j] == 0 && abs(d1 / d2) <= lambda) {
	    beta[j] <- 0
	 } else if(d * beta[j] < 0) {
	    beta[j] <- 0
	 } else {
	    beta[j] <- d
	 }
      }
   }
   beta
}


source("lasso.R")

N <- 100
p <- 1500
K <- 1
grp <- 1

# scale to zero mean and unit norm (not unit variance)
X0 <- matrix(rnorm(N * p), N, p)
X1 <- sweep(X0, 2, colMeans(X0))
X <- sweep(X1, 2, sqrt(diag(crossprod(X1))), FUN="/")
B <- matrix(rnorm(p * K), p, K)
B[sample(p, p - 10),] <- 0
Y <- X %*% B

system.time({
   g1 <- groupridge(X, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
      eps=1e-4, verbose=TRUE)
})
(err1 <- mean((X %*% g1 - Y)^2))
system.time({
   g2 <- groupridge(X, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
      eps=1e-10, verbose=TRUE)
})
(err2 <- mean((X %*% g2 - Y)^2))
system.time({
   g3 <- groupridge2(X, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
      eps=1e-10, verbose=TRUE)
})
(err3 <- mean((X %*% g3 - Y)^2))
g4 <- lasso(X, drop(Y), lambda1=1e-2, maxepoch=100)
(err4 <- mean((X %*% g4 - Y)^2))
g5 <- lasso.activeset(X, drop(Y), lambda1=1e-2)
(err5 <- mean((X %*% g5 - Y)^2))
g <- glmnet(X, drop(Y), lambda=1e-2)
b <- as.matrix(coef(g))[-1, ]
g6 <- groupridge_simple(X, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
      eps=1e-16,  verbose=TRUE)
(err6 <- mean((X %*% g6 - Y)^2))
g7 <- simpler(X, drop(Y), lambda=1e-2)
(err7 <- mean((X %*% g7 - Y)^2))

(r <- cbind(True=drop(B), glmnet=b, gr1=drop(g1), gr2=drop(g2), gr3=drop(g3),
      l1=g4, l2=g5, simple=drop(g6), simpler=drop(g7)))
cor(r)

