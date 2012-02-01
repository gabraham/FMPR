#

R2 <- function(pr, y) 1 - sum((pr - y)^2) / sum((y - mean(y))^2)
mse <- function(pr, y) mean((pr - y)^2)
center <- function(x) scale(x, center=TRUE, scale=FALSE)

crossval.lasso <- function(X, Y, nfolds=5, fun=R2, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- lasso(scale(X[folds != fold, ]), center(Y[folds != fold, ]), ...)
      p <- scale(X[folds == fold, ]) %*% g

      # one R^2 for each penalty
      apply(p, 2, fun, y=center(Y[folds == fold, ]))
   })
   rowMeans(s)
}

crossval.ridge <- function(X, Y, nfolds=5, fun=R2, lambda2)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- ridge(scale(X[folds != fold, ]), center(Y[folds != fold, ]),
	    lambda2=lambda2)

      sapply(g, function(m) {
	 p <- scale(X[folds == fold, ]) %*% m
	 R2(p, as.numeric(center(Y[folds == fold, ])))
      })
   })
   rowMeans(s)
}

crossval.fmpr <- function(X, Y, nfolds=5, G=NULL,
      lambda1, lambda2, lambda3, type, maxiter=1e5, verbose=FALSE)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   res <- array(0,
      c(nfolds, length(lambda1), length(lambda2), length(lambda3))
   )

   g <- lapply(1:nfolds, function(fold) {
      cat("fold", fold, " ")
      fmpr(scale(X[folds != fold, ]),
	    center(Y[folds != fold, ]), G=G, maxiter=maxiter,
	    lambda1=lambda1, lambda2=lambda2, lambda3=lambda3,
	    type=type, verbose=verbose)
   })
   cat("\n")

   for(fold in 1:nfolds)
   {
      Xtest <- scale(X[folds == fold, ])
      Ytest <- as.numeric(center(Y[folds == fold, ]))

      for(i in 1:length(lambda1))
      {
	 for(j in 1:length(lambda2))
	 {
	    for(k in 1:length(lambda3))
	    {
	       p <- Xtest %*% g[[fold]][[i]][[j]][[k]]
	       r <- R2(as.numeric(p), Ytest)
	       res[fold, i, j, k] <- r
	    }
	 }
      }
   }
   apply(res, c(2, 3, 4), mean, na.rm=TRUE)
}

crossval.elnet.glmnet <- function(X, Y, nfolds=5, fun=R2, lambda1, alpha, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)

   res <- array(NA, c(nfolds, length(lambda1), length(alpha)))

   for(fold in 1:nfolds)
   {
      for(i in 1:length(alpha))
      {
         g <- glmnet(X[folds != fold, ], Y[folds != fold], lambda=lambda1,
	       alpha=alpha[i], ...)
	 b <- as.matrix(coef(g))
         p <- cbind(1, X[folds == fold, ]) %*% b

	 # one R^2 for each penalty
	 res[fold, 1:ncol(b), i] <- apply(p, 2, fun, y=Y[folds == fold, ])
      }
   }
   apply(res, c(2, 3), mean, na.rm=TRUE)
}

optim.fmpr <- function(X, Y, nfolds, G=NULL, grid=20,
   L1=seq(0, 10, length=grid),
   L2=seq(0, 10, length=grid),
   L3=seq(0, 10, length=grid),
   type, maxiter=1e5, verbose=FALSE)
{
   # The L2 and L3 penalties must be in increasing order.
   # L1 should be in decreasing order.
   # This allows for certain optimisations in fmpr.
   L1 <- sort(L1, decreasing=TRUE)
   L2 <- sort(L2)
   L3 <- sort(L3)

   if(is.null(G) || nrow(G) == 1)
      L3 <- 0

   r <- crossval.fmpr(X=X, Y=Y, nfolds=nfolds,
      lambda1=L1, lambda2=L2, lambda3=L3,
      G=G, maxiter=maxiter, type=type, verbose=verbose)

   w <- which(r == max(r, na.rm=TRUE), arr.ind=TRUE)
   list(
      R2=r,
      opt=c(L1=L1[w[1]], L2=L2[w[2]], L3=L3[w[3]])
   )
}

optim.fmpr.2 <- function(X, Y, nfolds, G=NULL, grid=20,
   L1, L2, L3,
   type, maxiter=1e5, verbose=FALSE)
{
   cv <- function()
   {
      r <- crossval.fmpr(X=X, Y=Y, nfolds=nfolds,
	 lambda1=L1, lambda2=L2, lambda3=L3,
	 G=G, maxiter=maxiter, type=type, verbose=verbose)
      max(r, na.rm=TRU)
   }
   opt <- optim(c(L1[1], L2[1], L3[1]), fr=cv, method="L-BFGS-B",
      lower=c(L1[1], L2[1], L3[1]), upper=c(L1[2], L2[2], L3[2]))
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

optim.elnet.glmnet <- function(X, Y, nfolds, L1, alpha)
{  
   r <- crossval.elnet.glmnet(X=X, Y=Y, nfolds=nfolds,
	 lambda1=L1, alpha=alpha)
   w <- which(r == max(r, na.rm=TRUE), arr.ind=TRUE)

   list(
      R2=r,
      opt=c(L1=L1[w[1]], alpha=alpha[w[2]])
   )
}
