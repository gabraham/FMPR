#

R2 <- function(pr, y) 1 - sum((pr - y)^2) / sum((y - mean(y))^2)

crossval.groupridge <- function(X, Y, nfolds=5, G,
      lambda1, lambda2, lambda3, maxiter)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- groupridge4(X[folds != fold, ], Y[folds != fold, ], G=G,
	 maxiter=maxiter, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3)
      p <- X[folds == fold, ] %*% g
      R2(as.numeric(p), as.numeric(Y[folds == fold, ]))
   })
   mean(s)
}

crossval.lasso <- function(X, Y, nfolds=5, lambda2, maxiter)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- ridge(X[folds != fold, ], Y[folds != fold, ], lambda2=lambda2)
      p <- X[folds == fold, ] %*% g
      R2(as.numeric(p), as.numeric(Y[folds == fold, ]))
   })
   mean(s)
}

crossval.lasso <- function(X, Y, nfolds=5, lambda1, maxiter)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- lasso3(X[folds != fold, ], Y[folds != fold, ], maxiter=maxiter,
	    lambda1=lambda1)
      p <- X[folds == fold, ] %*% g
      R2(as.numeric(p), as.numeric(Y[folds == fold, ]))
   })
   mean(s)
}

crossval.ridge <- function(X, Y, nfolds=5, lambda2)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- ridge(X[folds != fold, ], Y[folds != fold, ], lambda2=lambda2)
      p <- X[folds == fold, ] %*% g
      R2(as.numeric(p), as.numeric(Y[folds == fold, ]))
   })
   mean(s)
}

optim.groupridge <- function(X, Y, nfolds, G, grid=20,
   L1=seq(0, 10, length=grid),
   L2=seq(0, 10, length=grid),
   L3=seq(0, 10, length=grid), maxiter)
{
   r <- array(0, dim=c(length(L1), length(L2), length(L3)))

   for(i in seq(along=L1)) {
      for(j in seq(along=L2)) {
	 for(k in seq(along=L3)) {
	    cat(i, j, k, "\n")
	    r[i, j, k] <- crossval.groupridge(X=X, Y=Y, nfolds=nfolds,
		  lambda1=L1[i], lambda2=L2[j], lambda3=L3[k],
		  G, maxiter=maxiter)
	 }
      }
   }

   w <- which(r == max(r, na.rm=TRUE), arr.ind=TRUE)
   list(
      R2=r,
      opt=c(L1=L1[w[1]], L2=L2[w[2]], L3=L3[w[3]])
   )
}

optim.lasso <- function(X, Y, nfolds, grid=20,
   L1=seq(0, 10, length=grid), maxiter)
{
   r <- numeric(length(L1))

   for(i in seq(along=L1)) {
      r[i] <- crossval.lasso(X=X, Y=Y, nfolds=nfolds,
	    lambda1=L1[i], maxiter=maxiter)
   }

   w <- which(r == max(r, na.rm=TRUE))
   list(
      R2=r,
      opt=c(L1=L1[w[1]])
   )
}

optim.ridge <- function(X, Y, nfolds, grid=20, L2=seq(0, 10, length=grid))
{
   r <- numeric(length(L2))

   for(i in seq(along=L2)) {
      r[i] <- crossval.ridge(X=X, Y=Y, nfolds=nfolds, lambda2=L2[i])
   }

   w <- which(r == max(r, na.rm=TRUE))
   list(
      R2=r,
      opt=c(L2=L2[w[1]])
   )
}

