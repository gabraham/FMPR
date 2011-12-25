
library(MASS)
library(glmnet)

options(error=dump.frames)

#seed <- sample(1e9, 1)
#set.seed(seed)
#set.seed(4398790)

dyn.load("groupridge.so")

groupridge <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g,
      maxiter=1e6, eps=1e-4)
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

maxlambda1 <- function(X, Y)
{
   r <- .C("maxlambda1", as.numeric(X), as.numeric(Y),
      numeric(ncol(Y)), nrow(X), ncol(X), ncol(Y))
   r[[3]]
}

blockX <- function(X, p, C)
{
   N <- nrow(X)
   
   Xblock <- matrix(0, N * C, p * C) 
   s1 <- seq(1, N * C, N)
   s2 <- seq(1, p * C, p)

   for(i in 1:C)
   {
      j <- s1[i]:(s1[i] + N - 1)
      k <- s2[i]:(s2[i] + p - 1)
      Xblock[j, k] <- X
   }

   Xblock
}

crossval <- function(X, Y, grp, nfolds=10)
{
   N <- nrow(X)
   K <- ncol(Y)
   p <- ncol(X)

   folds <- sample(nfolds, N, TRUE)
   lapply(1:nfolds, function(fold) {
      Xtrain <- X[folds != fold, ]
      Ytrain <- Y[folds != fold, ]
      Xtest <- X[folds == fold, ]
      Ytest <- Y[folds == fold, ]

      maxL <- maxlambda1(Xtrain, Ytrain)
      #B1 <- groupridge(Xtrain, Ytrain, lambda1=0, lambda2=0, lambda3=0, g=g)
      #mean((B0 - B1)^2)
      L2 <- 10^seq(-6, 6, length=20)
      # highest L1 penalty that makes all coefs zero, can also select
      # per task
      L1 <- seq(1, 0.1, length=20) * max(maxL)

      R <- lapply(L1, function(l1) {
         lapply(L2, function(l2) {
            groupridge(Xtrain, Ytrain,
		  lambda1=l1, lambda2=0, lambda3=l2, g=grp)
         })
      })
      
      r <- lapply(R, function(r) {
         sapply(r, function(x) rowMeans(x[1:p, , drop=FALSE]))
      })
      
      P <- lapply(R, function(r) {
         lapply(r, function(m) {
            Xtest %*% m
         })
      })
      
      R2 <- lapply(P, function(p) {
         sapply(p, function(r) {
            mean(sapply(1:K, function(k) {
               1 - mean((r[,k] - Ytest[,k])^2) / var(Ytest[, k])
            }))
         })
      })

      #m <- do.call(cbind, R2)
      #w <- which(m == max(m), arr.ind=TRUE)
      #c(L2=L2[w[1]], L1=L1[w[2]], R2=m[w[1], w[2]] )
      list(L1=L1, L2=L2, R2=R2)
   })
}

run <- function(N=100, p=200, K=5, M=sample(K, 1))
{
   Xtrain <- scale(matrix(sample(0:2, N * p, replace=TRUE), N, p))
   Xtest <- scale(matrix(sample(0:2, N * p, replace=TRUE), N, p))
   grp <- sample(M, K, TRUE)
   B <- matrix(sample(0:1, p * K, TRUE, prob=c(0.9, 0.1)), p, K)
   Ytrain <- scale(Xtrain %*% B + rnorm(N, 0, 1), scale=FALSE)
   Ytest <- scale(Xtest %*% B + rnorm(N, 0, 1), scale=FALSE)
   Ynet <- outer(grp, grp, function(x, y) as.integer(x == y))
   
   write.table(Xtrain, file="X.txt",
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(Ytrain, file="Y.txt",
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(Ynet, file="Ynetwork.txt",
         col.names=FALSE, row.names=FALSE, quote=FALSE)

   XtrainB <- blockX(Xtrain, p, K)
   XtestB <- blockX(Xtest, p, K)
   ytrain <- as.numeric(Ytrain)
   ytest <- as.numeric(Ytest)
   g <- cv.glmnet(XtrainB, ytrain, nfolds=10)
   B.lasso <- as.matrix(coef(g))[-1, ]
   P.lasso <- XtestB %*% B.lasso
   R2.lasso <- 1 - mean((P.lasso - ytest)^2) / var(ytest)
   
   res <- crossval(Xtrain, Ytrain, grp, nfolds=10)
   
   b <- lapply(res, function(r) {
      m <- do.call(cbind, r$R2)
      w <- which(m == max(m), arr.ind=TRUE)
      c(L2=r$L2[w[1]], L1=r$L1[w[2]], R2=m[w[1], w[2]])
   })
   
   b2 <- colMeans(do.call(rbind, b))
   
   Btrain <- groupridge(Xtrain, Ytrain, lambda1=b2[2],
         lambda2=0, lambda3=b2[1], g=grp)
   Ptest <- Xtest %*% Btrain
   R2test <- mean(sapply(1:K, function(k) {
      1 - mean((Ptest[,k] - Ytest[,k])^2) / var(Ytest[, k])
   }))
   
   system("./gflasso.sh")
   B.gfl <- matrix(scan("gflasso_betas.txt", what=numeric()), byrow=TRUE,
      nrow=p, ncol=K)
   P.gfl <- Xtest %*% B.gfl
   
   R2.gfl <- mean(sapply(1:K, function(k) {
      1 - mean((P.gfl[,k] - Ytest[,k])^2) / var(Ytest[, k])
   }))
   
   c(lasso=R2.lasso, l1l2=R2test, gflasso=R2.gfl)
}

res <- replicate(10, run())

