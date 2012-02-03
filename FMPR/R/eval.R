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

   # fmpr() is capable of handling lambda1,lambda2,lambda3 of varying lengths,
   # but we don't use that feature here because we want to run different
   # penalties on the cross-validation folds rather than different folds every
   # time
   for(fold in 1:nfolds)
   {
      cat("inner fold", fold, "\n")
      Xtest <- scale(X[folds == fold, , drop=FALSE])
      Ytest <- as.numeric(center(Y[folds == fold, , drop=FALSE]))

      f <- fmpr(scale(X[folds != fold, , drop=FALSE]),
	    center(Y[folds != fold, , drop=FALSE]),
	    G=G, maxiter=maxiter,
	    lambda1=lambda1, lambda2=lambda2, lambda3=lambda3,
	    type=type, verbose=verbose)

      for(i in 1:length(lambda1))
      {
	 for(j in 1:length(lambda2))
	 {
	    for(k in 1:length(lambda3))
	    {
	       p <- Xtest %*% f[[i]][[j]][[k]]
	       res[fold, i, j, k] <- cbind(R2(as.numeric(p), Ytest))
	    }
	 } 
      }
   }

   apply(res, c(2, 3, 4), mean, na.rm=TRUE)
}

crossval.spg <- function(X, Y, nfolds=5, C=NULL,
      lambda, gamma, maxiter=1e5, verbose)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   res <- array(0,
      c(nfolds, length(lambda), length(gamma))
   )

   for(fold in 1:nfolds)
   {
      cat("inner fold", fold, "\n")
      Xtest <- scale(X[folds == fold, , drop=FALSE])
      Ytest <- as.numeric(center(Y[folds == fold, , drop=FALSE]))

      f <- spg(scale(X[folds != fold, , drop=FALSE]),
	    center(Y[folds != fold, , drop=FALSE]),
	    C=C, maxiter=maxiter, verbose=verbose,
	    lambda=lambda, gamma=gamma)

      for(i in 1:length(lambda1))
      {
	 for(j in 1:length(lambda2))
	 {
	    p <- Xtest %*% f[[i]][[j]]
	    res[fold, i, j] <- cbind(R2(as.numeric(p), Ytest))
	 } 
      }
   }

   apply(res, c(2, 3), mean, na.rm=TRUE)
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

optim.fmpr <- function(..., method="grid")
{
   if(method == "search") {
      optim.fmpr.search(...)
   } else if(method == "grid") {
      optim.fmpr.grid(...)
   } else {
      stop("unknown optim method: ", method)
   }
}

optim.fmpr.grid <- function(X, Y, nfolds=5, G=NULL, grid=20,
   lambda1=seq(0, 10, length=grid),
   lambda2=seq(0, 10, length=grid),
   lambda3=seq(0, 10, length=grid),
   type, maxiter=1e5, verbose=FALSE)
{
   # The lambda2 and lambda3 penalties must be in increasing order.
   # lambda1 should be in decreasing order.
   # This allows for certain optimisations in fmpr.
   lambda1 <- sort(lambda1, decreasing=TRUE)
   lambda2 <- sort(lambda2)
   lambda3 <- sort(lambda3)

   if(is.null(G) || nrow(G) == 1)
      lambda3 <- 0

   r <- crossval.fmpr(X=X, Y=Y, nfolds=nfolds,
      lambda1=lambda1, lambda2=lambda2, lambda3=lambda3,
      G=G, maxiter=maxiter, type=type, verbose=verbose)

   w <- which(r == max(r, na.rm=TRUE), arr.ind=TRUE)
   list(
      R2=r,
      opt=c(lambda1=lambda1[w[1]], lambda2=lambda2[w[2]], lambda3=lambda3[w[3]])
   )
}

optim.fmpr.search <- function(X, Y, nfolds=5, G=NULL,
   lambda1, lambda2, lambda3,
   type, maxiter=1e5, verbose=FALSE)
{
  
   cv <- function(par)
   {
      r <- crossval.fmpr(X=X, Y=Y, nfolds=nfolds,
         lambda1=par[1], lambda2=par[2], lambda3=par[3],
         G=G, maxiter=maxiter, type=type, verbose=verbose)
      -max(r, na.rm=TRUE)
   }

   opt <- optim(par=c(min(lambda1), min(lambda2), min(lambda3)),
      fn=cv, method="L-BFGS-B",
      lower=c(min(lambda1), min(lambda2), min(lambda3)),
      upper=c(max(lambda1), max(lambda2), max(lambda3)),
      control=list(factr=1e4))

   list(
      R2=opt$value * -1,
      opt=opt$par
   )
}

optim.spg <- function(X, Y, nfolds=5, C=NULL, grid=20,
   lambda=seq(0, 10, length=grid),
   gamma=seq(0, 10, length=grid),
   maxiter=1e5, verbose=FALSE)
{
   lambda <- sort(lambda)
   gamma <- sort(gamma)

   if(is.null(C) || nrow(C) == 1)
      gamma <- 0

   r <- crossval.spg(X=X, Y=Y, nfolds=nfolds,
      lambda=lambda, gamma=gamma, C=C,
      maxiter=maxiter, verbose=verbose)

   w <- which(r == max(r, na.rm=TRUE), arr.ind=TRUE)
   list(
      R2=r,
      opt=c(lambda=lambda[w[1]], gamma=gamma[w[2]])
   )
}

optim.lasso <- function(X, Y, nfolds=5, grid=20,
   lambda1=seq(0, 10, length=grid), maxiter=1e5, verbose=FALSE)
{
   r <- crossval.lasso(X=X, Y=Y, nfolds=nfolds,
	    lambda1=lambda1, maxiter=maxiter, verbose=verbose)

   w <- which(r == max(r, na.rm=TRUE))
   list(
      R2=r,
      opt=c(lambda1=lambda1[w[1]])
   )
}

optim.ridge <- function(X, Y, nfolds=5, grid=20, lambda2=seq(0, 10, length=grid))
{
   r <- crossval.ridge(X=X, Y=Y, nfolds=nfolds, lambda2=lambda2)

   w <- which(r == max(r, na.rm=TRUE))
   list(
      R2=r,
      opt=c(lambda2=lambda2[w[1]])
   )
}

optim.elnet.glmnet <- function(X, Y, nfolds=5, lambda1, alpha)
{  
   r <- crossval.elnet.glmnet(X=X, Y=Y, nfolds=nfolds,
	 lambda1=lambda1, alpha=alpha)
   w <- which(r == max(r, na.rm=TRUE), arr.ind=TRUE)

   list(
      R2=r,
      opt=c(lambda1=lambda1[w[1]], alpha=alpha[w[2]])
   )
}

