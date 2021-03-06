
which.max2 <- function(x, dup=FALSE, dup.idx=1)
{
   w <- which(max(x) == x, arr.ind=TRUE)
   if(!dup) {
      w[dup.idx, ]
   } else {
      w
   }
}

R2 <- function(pr, y, average=TRUE) 
{
   pr <- cbind(pr)
   y <- cbind(y)
   if(ncol(y) != ncol(pr))
      stop("ncol(y) doesn't match ncol(pr)")

   s <- sapply(1:ncol(y), function(k) {
      1 - sum((pr[,k] - y[,k])^2) / sum((y[,k] - mean(y[,k]))^2)
   })
   s[is.nan(s)] <- 0
   if(average) {
      mean(s)
   } else {
      s
   }
}

mse <- function(pr, y) mean((pr - y)^2)
center <- function(x) scale(x, center=TRUE, scale=FALSE)

# scale and set NAs caused by zero variance to zero
scalefix <- function(X)
{
   s <- scale(X)
   s[is.na(s)] <- 0
   s
}

# FMPR can do ridge regression but it's faster to do it using specialised code
crossval.ridge <- function(X, Y, nfolds=5, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   
   l <- foreach(fold=1:nfolds) %dopar% {
      cat("inner fold", fold, "\n")
      g <- ridge(scalefix(X[folds != fold, ]),
	 scalefix(Y[folds != fold, , drop=FALSE]), ...)
      sapply(g, function(m) {
	 if(is.null(m)) {
	    0
	 } else {
	    p <- scalefix(X[folds == fold, , drop=FALSE]) %*% m
	    R2(p, scalefix(Y[folds == fold, , drop=FALSE]))
	 }
      })
   }

   s <- do.call(cbind, l)

   list(
      R2=rowMeans(s),
      lambda=list(...)$lambda
   )
}

crossval.fmpr <- function(X, Y, nfolds=5, cortype=2,
   cormethod="pearson", ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   K <- ncol(Y)
   folds <- sample(1:nfolds, N, TRUE)

   l <- max(length(list(...)$lambda), 1)
   l2 <- max(length(list(...)$lambda2), 1)
   g <- max(length(list(...)$gamma), 1)
   verbose <- list(...)$verbose
   if(is.null(verbose)) {
      verbose <- FALSE
   }
   res <- array(0, c(nfolds, l, l2, g, K))
   nz <- array(0, c(nfolds, l, l2, g, K))

   # fmpr() is capable of handling lambda,gamma of varying lengths,
   # but we don't use that feature here because we want to run different
   # penalties on the same cross-validation folds rather than different
   # folds every time
   resl <- foreach(fold=1:nfolds) %dopar% {
      r <- array(0, c(l, l2, g, K))
      nzl <- array(0, c(l, l2, g, K))

      if(verbose) {
	 cat("inner fold", fold, "\n")
      }
      Xtest <- scalefix(X[folds == fold, , drop=FALSE])
      Ytest <- scalefix(Y[folds == fold, , drop=FALSE])
      Xtrain <- scalefix(X[folds != fold, , drop=FALSE])
      Ytrain <- scalefix(Y[folds != fold, , drop=FALSE])
      C <- NULL
      if(ncol(Y) > 1) {
	 C <- gennetwork(Ytrain, cortype=cortype,
	    cormethod=cormethod)
      }

      f <- fmpr(X=Xtrain, Y=Ytrain, C=C, ...)
      for(i in 1:l)
      {
	 for(m in 1:l2)
	 {
	    for(j in 1:g)
	    {
	       B <- f[[i]][[m]][[j]]
	       if(is.null(B)) {
		  r[i, m, j, ] <- 0
		  nzl[i, m, j, ] <- 0
	       } else {
		  p <- as.matrix(Xtest %*% B)
		  r[i, m, j, ] <- cbind(R2(p, Ytest, average=FALSE))
		  nzl[i, m, j, ] <- colSums(B != 0)
	       }
	    }
	 }
      }
      list(r=r, nzl=nzl)
   }
 
   for(fold in 1:nfolds) {
      res[fold, , , , ] <- resl[[fold]]$r
      nz[fold, , , , ] <- resl[[fold]]$nzl
   }

   list(
      R2=apply(res, c(2, 3, 4), mean, na.rm=TRUE),
      res=res,
      nz=nz,
      lambda=list(...)$lambda,
      lambda2=list(...)$lambda2,
      gamma=list(...)$gamma
   )
}

crossval.spg <- function(X, Y, nfolds=5, cortype=2, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)

   ll <- max(length(list(...)$lambda), 1)
   lg <- max(length(list(...)$gamma), 1)
   res <- array(0, c(nfolds, ll, lg))

   resl <- foreach(fold=1:nfolds) %dopar% {
      r <- matrix(0, ll, lg)
      cat("inner fold", fold, "\n")

      Xtrain <- scalefix(X[folds != fold, , drop=FALSE])
      Ytrain <- scalefix(Y[folds != fold, , drop=FALSE])
      Xtest <- scalefix(X[folds == fold, , drop=FALSE])
      Ytest <- scalefix(Y[folds == fold, , drop=FALSE])
      C <- gennetwork(Ytrain, cortype)
      f <- spg(Xtrain, Ytrain, C=C, ...)

      for(i in 1:ll)
      {
	 for(j in 1:lg)
	 {
	    p <- Xtest %*% f[[i]][[j]]
	    r[i, j] <- cbind(R2(p, Ytest))
	 } 
      }
      r
   }

   for(fold in 1:nfolds) {
      res[fold, , ] <- resl[[fold]]
   }

   list(
      R2=apply(res, c(2, 3), mean, na.rm=TRUE),
      lambda=list(...)$lambda,
      gamma=list(...)$gamma
   )
}

optim.fmpr <- function(method="grid", ...)
{
   if(method == "search") {
      optim.fmpr.search(...)
   } else if(method == "grid") {
      optim.fmpr.grid(...)
   } else if(method == "random") {
      optim.fmpr.random(...)
   } else {
      stop("unknown optim method: ", method)
   }
}

optim.fmpr.random2 <- function(...)
{
   r <- crossval.fmpr(...)

   w <- rbind(which(r$R2 == max(r$R2, na.rm=TRUE), arr.ind=TRUE))[1,]
   list(
      R2=r$R2,
      opt=c(
	 lambda=r$lambda[w[1]],
	 lambda2=r$lambda2[w[2]],
	 gamma=r$gamma[w[3]]
      )
   )
}

optim.fmpr.grid <- function(...)
{
   r <- crossval.fmpr(...)

   K <- ncol(list(...)$Y)

   res <- apply(r$res, 2:5, mean)
   nz <- apply(r$nz, 2:5, mean)
   idxmax <- apply(res, 4, which.max2)
   colnames(idxmax) <- paste("Task", 1:K, sep="")
   rownames(idxmax) <- c("lambda", "lambda2", "gamma")
   nzmax <- nz[t(rbind(idxmax, 1:K))]
   resmax <- res[t(rbind(idxmax, 1:K))]

   w <- rbind(which(r$R2 == max(r$R2, na.rm=TRUE), arr.ind=TRUE))[1,]
   list(
      lambda=r$lambda,
      lambda2=r$lambda2,
      gamma=r$gamma,
      R2=r$R2,
      res=r$res,
      nz=r$nz,
      opt=c(
	 lambda=r$lambda[w[1]],
	 lambda2=r$lambda2[w[2]],
	 gamma=r$gamma[w[3]]
      ),
      opt2=list(
	 res=res,
	 nz=nz,
	 idxmax=idxmax,
	 nzmax=nzmax,
	 resmax=resmax
      )
   )
}

optim.fmpr.search <- function(X, Y, nfolds=5,
   cortype=1, lambda=c(10, 0), gamma=c(10, 0), lambda2=c(10, 0),
   maxiter=1e5, verbose=FALSE)
{
  
   cv <- function(par)
   {
      r <- crossval.fmpr(X=X, Y=Y, nfolds=nfolds,
         lambda=par[1], gamma=par[2], lambda2=par[3],
	 cortype=cortype, maxiter=maxiter, verbose=verbose)
      -max(r$R2, na.rm=TRUE)
   }

   opt <- optim(par=c(min(lambda), min(gamma)),
      fn=cv, method="L-BFGS-B",
      lower=c(min(lambda), min(gamma), min(lambda2)),
      upper=c(max(lambda), max(gamma), max(lambda2)),
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

optim.ridge <- function(...)
{
   r <- crossval.ridge(...)

   w <- rbind(which(r$R2 == max(r$R2, na.rm=TRUE)))[1,]
   list(
      R2=r$R2,
      opt=c(lambda=r$lambda[w[1]])
   )
}

