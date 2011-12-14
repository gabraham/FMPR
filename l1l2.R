  
options(error=dump.frames)

library(MASS)

#set.seed(393218)

N <- 100
p <- 20
K <- 1
X <- scale(matrix(rnorm(N * p), N, p))
b <- rnorm(p) #* sample(0:1, p, replace=TRUE)
# same B for all groups
B <- matrix(rep(b, K), p, K)
Y <- X %*% B + rnorm(N * K)

soft <- function(a, b)
{
   sign(a) * pmax(0, abs(a) - b) 
}

# ridge regression using coordinate descent, proportional shrinkage
# in each variable
ridge <- function(X, Y, lambda2=0, maxiter=1000)
{
   n <- nrow(X)
   p <- ncol(X)
   B <- numeric(p)
   L <- numeric(N)
   E <- mean((L - Y)^2)
   Eold <- Inf

   for(iter in 1:maxiter)
   {
      # over all variables
      for(j in 1:p)
      {
	 s <- drop((crossprod(X[,j], L - Y) + lambda2 * B[j]) /
	       (crossprod(X[,j]) + lambda2))
	 Bj <- B[j] - s
	 D <- Bj - B[j]
	 B[j] <- Bj
	 L <- L + D * X[,j]
      }
      E <- mean((L - Y)^2)
      if(abs(E - Eold) <= 1e-6)
	 break
      Eold <- E
   }
   cat("ridge done in", iter, "iterations\n")
   B
}

l1l2 <- function(X, Y, lambda1=0, lambda2=0, maxiter=100)
{
   n <- nrow(X)
   p <- ncol(X)
   K <- ncol(Y)
   Bhat <- matrix(0, p, K)
   L <- matrix(0, N, K)

   for(iter in 1:maxiter)
   {
      # over all variables
      for(j in 1:p)
      {
	 #nrmj <- sum(B[j,] * B[j,])# / sqrt(length())

	 Bhat2 <- Bhat
	 L2 <- L

	 # first we apply the l1 penalty to get B sparse in each column k=1:K,
	 # i.e., K univariable steps
	 for(k in 1:K)
	 {
	    s <- sum(X[,j] * (L2[,k] - Y[,k])) / sum(X[,j]^2)
	    Bold <- Bhat2[j, k]
	    Bhat2[j, k] <- soft(Bhat2[j, k] - s, lambda1)
	    browser()
	    D <- Bhat2[j, k] - Bold
	    L2[,k] <- L2[,k] + D * X[,j]
	 }

	 L <- L + X[, j] %*% S

      }
   }
   Bhat
}

# grp: a
groupridge <- function(X, Y, grp, lambda=0)
{
   N <- nrow(X)
   p <- ncol(X)
   K <- ncol(Y)
   B <- matrix(0, p, K)
   L <- matrix(0, N, K)

   for(k in 1:K)
   {
      for(j in 1:p)
      {
	 
      }
      L <- L + X[,]
   }
   B
}

# When lambda1=0, should be identical to least squares
lambda2 <- 10
#r.svd <- drop(ginv(X) %*% Y)
r.ridge1 <- ridge(X, Y, lambda2=lambda2)
r.ridge2 <- drop(qr.solve(crossprod(X) + lambda2 * diag(p), crossprod(X, Y)))
diag(cor(cbind(r.ridge1, r.ridge2)))
#mean((r.svd - r.ridge1)^2)
#mean((r.svd - r.ridge2)^2)
mean((r.ridge1 - r.ridge2)^2)
#r1 <- l1l2(X, Y)
#mean((r1 - r.ridge)^2)
#apply(r1, 1, crossprod)
#
## Just shrink across outputs, no sparsity
#r2 <- l1l2(X, Y, lambda2=1000)
#mean((r2 - r.ridge)^2)
#apply(r2, 1, crossprod)
#
