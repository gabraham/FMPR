#

R2 <- function(pr, y) 1 - sum((pr - y)^2) / sum((y - mean(y))^2)
mse <- function(pr, y) mean((pr - y)^2)
center <- function(x) scale(x, center=TRUE, scale=FALSE)

crossval.groupridge <- function(X, Y, nfolds=5, G=NULL,
      lambda1, lambda2, lambda3, maxiter=1e5)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   res <- array(0,
      c(nfolds, length(lambda1), length(lambda2), length(lambda3))
   )

   g <- lapply(1:nfolds, function(fold) {
      cat("fold", fold, "\n")
      groupridge4(scale(X[folds != fold, ]),
	    center(Y[folds != fold, ]), G=G, maxiter=maxiter,
	    lambda1=lambda1, lambda2=lambda2, lambda3=lambda3)
   })

   for(fold in 1:nfolds)
   {
      for(i in 1:length(lambda1))
      {
	 for(j in 1:length(lambda2))
	 {
	    for(k in 1:length(lambda3))
	    {
	       p <- scale(X[folds == fold, ]) %*% g[[fold]][[i]][[j]][[k]]
	       r <- R2(as.numeric(p), as.numeric(center(Y[folds == fold, ])))
	       res[fold, i, j, k] <- r
	    }
	 }
      }
   }
   apply(res, c(2, 3, 4), mean, na.rm=TRUE)
}

crossval.lasso <- function(X, Y, nfolds=5, fun=R2, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- lasso3(scale(X[folds != fold, ]), center(Y[folds != fold, ]), ...)
      p <- scale(X[folds == fold, ]) %*% g
      apply(p, 2, fun, y=center(Y[folds == fold, ]))
   })
   rowMeans(s)
}

crossval.glmnet <- function(X, Y, nfolds=5, fun=R2, ...)
{
   N <- nrow(X)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- glmnet(X[folds != fold, ], y[folds != fold], ...)
      p <- cbind(1, X[folds == fold, ]) %*% coef(g)
      apply(p, 2, fun, y=Y[folds == fold, ])
   })
   rowMeans(s)
}

crossval.ridge <- function(X, Y, nfolds=5, fun=R2, lambda2)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- try(ridge(scale(X[folds != fold, ]), center(Y[folds != fold, ]),
	    lambda2=lambda2))
      if(is(g, "try-error")) {
	 -Inf
      } else {
	 p <- scale(X[folds == fold, ]) %*% g
	 apply(p, 2, fun, y=center(Y[folds == fold, ]))
      }
   })
   rowMeans(s)
}

optim.groupridge <- function(X, Y, nfolds, G=NULL, grid=20,
   L1=seq(0, 10, length=grid),
   L2=seq(0, 10, length=grid),
   L3=seq(0, 10, length=grid), maxiter=1e5, verbose=FALSE)
{
   if(is.null(G) || nrow(G) == 1)
      L3 <- 0

   r <- crossval.groupridge(X=X, Y=Y, nfolds=nfolds,
      lambda1=L1, lambda2=L2, lambda3=L3,
      G=G, maxiter=maxiter)

   w <- which(r == max(r, na.rm=TRUE), arr.ind=TRUE)
   list(
      R2=r,
      opt=c(L1=L1[w[1]], L2=L2[w[2]], L3=L3[w[3]])
   )
}

optim.lasso <- function(X, Y, nfolds, grid=20,
   L1=seq(0, 10, length=grid), maxiter=1e5, verbose=FALSE)
{
   r <- crossval.lasso(X=X, Y=Y, nfolds=nfolds,
	    lambda1=L1, maxiter=maxiter, verbose=verbose)

   w <- which(r == max(r, na.rm=TRUE))
   list(
      R2=r,
      opt=c(L1=L1[w[1]])
   )
}

optim.ridge <- function(X, Y, nfolds, grid=20, L2=seq(0, 10, length=grid))
{
   r <- crossval.ridge(X=X, Y=Y, nfolds=nfolds, lambda2=L2)

   w <- which(r == max(r, na.rm=TRUE))
   list(
      R2=r,
      opt=c(L2=L2[w[1]])
   )
}

