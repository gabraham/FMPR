#

R2 <- function(pr, y) 
{
   pr <- cbind(pr)
   y <- cbind(y)
   if(ncol(y) != ncol(pr))
      stop("ncol(y) doesn't match ncol(pr)")

   s <- sapply(1:ncol(y), function(k) {
      1 - sum((pr[,k] - y[,k])^2) / sum((y[,k] - mean(y[,k]))^2)
   })
   mean(pmax(0, s))
}

mse <- function(pr, y) mean((pr - y)^2)
center <- function(x) scale(x, center=TRUE, scale=FALSE)

crossval.lasso <- function(X, Y, nfolds=5, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- lasso(scale(X[folds != fold, ]), scale(Y[folds != fold, ]), ...)
      p <- scale(X[folds == fold, ]) %*% g

      # one R^2 for each penalty
      apply(p, 2, R2, y=scale(Y[folds == fold, ]))
   })
   
   list(
      R2=rowMeans(s),
      lambda1=list(...)$lambda1
   )
}

crossval.ridge <- function(X, Y, nfolds=5, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      cat("inner fold", fold, "\n")
      g <- ridge(scale(X[folds != fold, ]), scale(Y[folds != fold, ]), ...)

      sapply(g, function(m) {
	 p <- scale(X[folds == fold, ]) %*% m
	 R2(p, scale(Y[folds == fold, ]))
      })
   })
   list(
      R2=rowMeans(s),
      lambda2=list(...)$lambda2
   )
}

crossval.fmpr <- function(X, Y, nfolds=5,
      graph.fun=graph.sqr, graph.thresh=0, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)

   l1 <- max(length(list(...)$lambda1), 1)
   l2 <- max(length(list(...)$lambda2), 1)
   l3 <- max(length(list(...)$lambda3), 1)
   verbose <- list(...)$verbose
   res <- array(0, c(nfolds, l1, l2, l3))

   # fmpr() is capable of handling lambda1,lambda2,lambda3 of varying lengths,
   # but we don't use that feature here because we want to run different
   # penalties on the cross-validation folds rather than different folds every
   # time
   for(fold in 1:nfolds)
   {
      cat("inner fold", fold, "\n")
      Xtest <- scale(X[folds == fold, , drop=FALSE])
      Ytest <- scale(Y[folds == fold, , drop=FALSE])
      Xtrain <- scale(X[folds != fold, , drop=FALSE])
      Ytrain <- scale(Y[folds != fold, , drop=FALSE])
      G <- 0
      if(ncol(Y) > 1)
      {
	 R <- cor(Ytrain)
	 G <- graph.fun(R, graph.thresh)
      }

      f <- fmpr(X=scale(X[folds != fold, , drop=FALSE]),
	    Y=scale(Y[folds != fold, , drop=FALSE]), G=G, ...)

      for(i in 1:l1)
      {
	 for(j in 1:l2)
	 {
	    for(k in 1:l3)
	    {
	       p <- as.matrix(Xtest %*% f[[i]][[j]][[k]])
	       #res[fold, i, j, k] <- cbind(R2(as.numeric(p), Ytest))
	       res[fold, i, j, k] <- cbind(R2(p, Ytest))

	    }
	 } 
      }
      #if(verbose)
	 cat("best R2:", m <- max(res[fold, , ,], na.rm=TRUE), "\n")
   }

   list(
      R2=apply(res, c(2, 3, 4), mean, na.rm=TRUE),
      lambda1=list(...)$lambda1,
      lambda2=list(...)$lambda2,
      lambda3=list(...)$lambda3
   )
}

crossval.spg <- function(X, Y, nfolds=5, cortype=2, corthresh=0, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)

   ll <- max(length(list(...)$lambda), 1)
   lg <- max(length(list(...)$gamma), 1)
   res <- array(0, c(nfolds, ll, lg))

   for(fold in 1:nfolds)
   {
      cat("inner fold", fold, "\n")

      Xtrain <- scale(X[folds != fold, , drop=FALSE])
      Ytrain <- scale(Y[folds != fold, , drop=FALSE])
      Xtest <- scale(X[folds == fold, , drop=FALSE])
      Ytest <- scale(Y[folds == fold, , drop=FALSE])
      C <- gennetwork(Ytrain, corthresh, cortype)
      f <- spg(Xtrain, Ytrain, C=C, ...)

      for(i in 1:ll)
      {
	 for(j in 1:lg)
	 {
	    p <- Xtest %*% f[[i]][[j]]
	    res[fold, i, j] <- cbind(R2(p, Ytest))
	 } 
      }
   }

   list(
      R2=apply(res, c(2, 3), mean, na.rm=TRUE),
      lambda=list(...)$lambda,
      gamma=list(...)$gamma
   )
}

# glmnet produces the solution path for multiple lambdas, for a given alpha
crossval.elnet.glmnet <- function(X, Y, nfolds=5, lambda1=NULL, alpha=1, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)

   l1 <- length(lambda1)
   la <- length(alpha)
   res <- array(NA, c(nfolds, l1, la))

   g <- vector("list", nfolds)
   lambda1 <- sort(lambda1, decreasing=TRUE)

   for(fold in 1:nfolds)
   {
      g[[fold]] <- foreach(i=1:la) %dopar% {
         r <- glmnet(X[folds != fold, ], Y[folds != fold],
	       lambda=lambda1, alpha=alpha[i])
	 as.matrix(coef(r))
      }
   }

   for(fold in 1:nfolds)
   {
      for(i in 1:la)
      {
	 b <- g[[fold]][[i]]
         p <- cbind(1, X[folds == fold, ]) %*% b

	 # one R^2 for each penalty
	 res[fold, 1:ncol(b), i] <- apply(p, 2, R2, y=Y[folds == fold, ])
      }
   }
   
   r2 <- apply(res, c(2, 3), mean, na.rm=TRUE)
  
   list(
      R2=r2,
      lambda1=lambda1,
      alpha=alpha
   )
}

optim.fmpr <- function(method="grid", ...)
{
   if(method == "search") {
      optim.fmpr.search(...)
   } else if(method == "grid") {
      optim.fmpr.grid(...)
   } else {
      stop("unknown optim method: ", method)
   }
}

optim.fmpr.grid <- function(...)
{
   r <- crossval.fmpr(...)

   w <- rbind(which(r$R2 == max(r$R2, na.rm=TRUE), arr.ind=TRUE))[1,]
   list(
      R2=r$R2,
      opt=c(
	 lambda1=r$lambda1[w[1]],
	 lambda2=r$lambda2[w[2]],
	 lambda3=r$lambda3[w[3]])
   )
}

optim.fmpr.search <- function(X, Y, nfolds=5,
   graph.thresh=0.5, graph.fun=sqr,
   lambda1=c(10, 0), lambda2=c(10, 0), lambda3=c(10, 0),
   maxiter=1e5, verbose=FALSE)
{
  
   cv <- function(par)
   {
      r <- crossval.fmpr(X=X, Y=Y, nfolds=nfolds,
         lambda1=par[1], lambda2=par[2], lambda3=par[3],
	 graph.thresh=graph.thresh, graph.fun=graph.fun,
         maxiter=maxiter, verbose=verbose)
      -max(r$R2, na.rm=TRUE)
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

optim.spg <- function(...)
{
   r <- crossval.spg(...)

   w <- rbind(which(r$R2 == max(r$R2, na.rm=TRUE), arr.ind=TRUE))[1,]
   list(
      R2=r$R2,
      opt=c(lambda=r$lambda[w[1]], gamma=r$gamma[w[2]])
   )
}

optim.lasso <- function(...)
{
   r <- crossval.lasso(...)

   w <- rbind(which(r$R2 == max(r$R2, na.rm=TRUE)))[1,]
   list(
      R2=r$R2,
      opt=c(lambda1=r$lambda1[w[1]])
   )
}

optim.ridge <- function(...)
{
   r <- crossval.ridge(...)

   w <- rbind(which(r$R2 == max(r$R2, na.rm=TRUE)))[1,]
   list(
      R2=r$R2,
      opt=c(lambda2=r$lambda2[w[1]])
   )
}

optim.elnet.glmnet <- function(...)
{  
   r <- crossval.elnet.glmnet(...)

   w <- rbind(which(r$R2 == max(r$R2, na.rm=TRUE), arr.ind=TRUE))[1,]

   list(
      R2=r$R2,
      opt=c(lambda1=r$lambda1[w[1]], alpha=r$alpha[w[2]])
   )
}

