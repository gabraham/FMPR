
library(MASS)


seed <- sample(1e9, 1)
set.seed(seed)

N <- 100
p <- 10
K <- 5
X <- matrix(rnorm(N * p), N, p)
B <- matrix(rnorm(p * K), p, K)
Y <- scale(X %*% B + rnorm(N * K), scale=FALSE)

g <- rep(1, K)
lambda <- 1e3
B0 <- ginv(X) %*% Y
#B1 <- groupridge(X, Y, lambda=1e1, g=g)
#
#diag(cor(B, B0))
#diag(cor(B, B1))
#diag(cor(B0, B1))

# g: a vector of length K partitioned into M values
#groupridge <- function(X, Y, lambda, g)
#{
   #N <- nrow(X)
   #p <- ncol(X)
   #K <- ncol(Y)
   M <- max(g)
   B1 <- matrix(0, p, K)
   S <- table(g)

   for(iter in 1:50)
   {
      for(k in 1:K)
      {
         for(j in 1:p)
         {
            d1 <- crossprod(X[, j], X %*% B1[, k] - Y[, k])
            d2 <- crossprod(X[, j]) 

	    # Penalise by the differences from all other variables
	    # in the same partition
	    w <- which(g == g[k])
	    # sum((B[j, k] - B[j, w])^2) / 2
	    d1l <- d1 + lambda * sum(B1[j, k] - B1[j, w])
	    d2l <- d2 + lambda * N

            B1[j, k] <- B1[j, k] - d1l / d2l
         }
      }
   }
   
#   B
#}

diag(cor(B0, B))
diag(cor(B1, B))
diag(cor(B0, B1))
summary(diag(tcrossprod(B0)))
summary(diag(tcrossprod(B1)))

