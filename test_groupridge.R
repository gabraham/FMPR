# Test groupridge

#library(MASS)
library(glmnet)
library(gdata)
library(ROCR)

dyn.load("groupridge.so")

options(error=dump.frames)

s <- sample(1e6L, 1)
set.seed(s)

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

groupridge_simple <- function(X, Y, lambda1=0, lambda2=0, lambda3=0,
      g=NULL, maxiter=1e3, eps=1e-6, verbose=FALSE)
{
   p <- ncol(X)
   K <- ncol(Y)

   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   if(length(lambda2) == 1)
      lambda2 <- rep(lambda2, K)

   if(length(lambda3) == 1)
      lambda3 <- rep(lambda3, K)

   if(is.null(g))
      g <- 1:K

   r <- .C("groupridge_simple", as.numeric(X), as.numeric(Y), 
      numeric(p * K), nrow(X), ncol(X), ncol(Y),
      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
      as.integer(g), as.integer(maxiter), as.double(eps),
      as.integer(verbose)
   )
   matrix(r[[3]], p, K)
}


groupridge2 <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g=NULL,
      maxiter=1e5, eps=1e-6, verbose=TRUE)
{
   p <- ncol(X)
   K <- ncol(Y)

   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   if(length(lambda2) == 1)
      lambda2 <- rep(lambda2, K)

   if(length(lambda3) == 1)
      lambda3 <- rep(lambda3, K)
      
   if(is.null(g))
      g <- 1:K

   r <- .C("groupridge2", as.numeric(X), as.numeric(Y), 
      numeric(p * K), nrow(X), ncol(X), ncol(Y),
      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
      as.integer(g), as.integer(maxiter),
      as.double(eps), as.integer(verbose))
   matrix(r[[3]], p, K)
}

# lambda1: scalar or K-vector
# lambda2: scalar or K-vector
# lambda3: scalar 
groupridge3 <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, grp=NULL,
      maxiter=1e5, eps=1e-6, verbose=FALSE)
{
   p <- ncol(X)
   Y <- cbind(Y)
   K <- ncol(Y)

   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   if(length(lambda2) == 1)
      lambda2 <- rep(lambda2, K)

   if(is.null(grp))
      grp <- 1:K

   r <- .C("groupridge3", as.numeric(X), as.numeric(Y), 
      numeric(p * K), nrow(X), ncol(X), K,
      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
      as.integer(grp), as.integer(maxiter),
      as.double(eps), as.integer(verbose), integer(1))
   matrix(r[[3]], p, K)
}

lasso2 <- function(X, y, lambda1=0,
      maxiter=1e5, eps=1e-6, verbose=TRUE)
{
   p <- ncol(X)

   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   r <- .C("lasso2", as.numeric(X), as.numeric(y), numeric(p), 
      nrow(X), ncol(X), as.numeric(lambda1),
      as.integer(maxiter), as.double(eps),
      as.integer(verbose)
   )
   matrix(r[[3]], p)
}

lasso3 <- function(X, y, lambda1=0,
      maxiter=1e5, eps=1e-6, verbose=FALSE)
{
   p <- ncol(X)

   r <- .C("lasso3", as.numeric(X), as.numeric(y), numeric(p), 
      nrow(X), ncol(X), as.numeric(lambda1),
      as.integer(maxiter), as.double(eps),
      as.integer(verbose)
   )
   matrix(r[[3]], p)
}

ridge <- function(X, Y, lambda2=0)
{
   XX <- crossprod(X)
   XY <- crossprod(X, Y)
   p <- ncol(X)
   qr.solve(XX + diag(p) * lambda2, XY)
}

# a very simple version of lasso, for sanity checking
simpler <- function(x, y, lambda, maxiter=2000)
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

# zero mean, unit norm (not unit variance)
standardise <- function(x)
{
   x1 <- sweep(x, 2, colMeans(x))
   v <- apply(x1, 2, function(z) sqrt(sum(z^2)))
   sweep(x1, 2, v, FUN="/")
}

#source("lasso.R")

maxlambda1 <- function(X, Y)
{
   Y <- cbind(Y)
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

R2 <- function(pr, y) 1 - sum((pr - y)^2) / sum((y - mean(y))^2)

#################################################################################
## Test groupridge in standard lasso setting, single task
#
#N <- 100
#p <- 500
#K <- 1
#grp <- 1
#
## scale to zero mean and unit norm (not unit variance)
#X0 <- matrix(rnorm(N * p), N, p)
#X <- standardise(X0)
#B <- matrix(rnorm(p * K), p, K)
#B[sample(p, p - 10),] <- 0
#Y <- X %*% B
#
#maxL1 <- maxlambda1(X, Y)
#
#g <- glmnet(X, Y)
#gl <- as.matrix(coef(g))[-1,]
#
#gr <- sapply(g$lambda * sqrt(N), function(l1) {
#   groupridge3(X, Y, lambda1=l1)
#})
#
#summary(diag(cor(cbind(gl, gr))))
#
#stop()

crossval <- function(X, Y, nfolds=5, fun, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- fun(X[folds != fold, ], Y[folds != fold, ], ...)
      p <- X[folds == fold, ] %*% g
      R2(as.numeric(p), as.numeric(Y[folds == fold, ]))
   })
   mean(s)
}

makedata <- function(rep, dir=".", N=100, p=50, K=5, nfolds=10)
{
   cat("rep", rep, "\n")

   Xtrain <- standardise(matrix(sample(0:2, N * p, replace=TRUE), N, p))
   Xtest <- standardise(matrix(sample(0:2, N * p, replace=TRUE), N, p))
   grp <- rep(1, K)
   XtrainB <- blockX(Xtrain, p, K)
   XtestB <- blockX(Xtest, p, K)

   XtrainBs <- standardise(XtrainB)
   XtestBs <- standardise(XtestB)
   
   b <- 1 * sample(0:1, p, TRUE, prob=c(0.8, 0.2))
   B <- sapply(1:K, function(k) b)
   write.table(B, file=sprintf("%s/B_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   
   Ytrain <- scale(Xtrain %*% B + rnorm(N * K, 0, 0.1), scale=FALSE)
   Ytest <- scale(Xtest %*% B + rnorm(N * K, 0, 0.1), scale=FALSE)
   ytrain <- as.numeric(Ytrain)
   ytest <- as.numeric(Ytest)
   Ynet <- outer(grp, grp, function(x, y) as.integer(x == y))
   lowerTriangle(Ynet) <- 0  
   diag(Ynet) <- 0
   
   write.table(Xtrain, file=sprintf("%s/Xtrain_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(Ytrain, file=sprintf("%s/Ytrain_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(Ynet, file=sprintf("%s/Ynetwork_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(XtrainBs, file=sprintf("%s/XtrainB_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(cbind(grp), file=sprintf("%s/grp_%s.txt", dir, rep),
	 col.names=FALSE, row.names=FALSE, quote=FALSE)

   write.table(Xtest, file=sprintf("%s/Xtest_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(Ytest, file=sprintf("%s/Ytest_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(XtestBs, file=sprintf("%s/XtestB_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
}

run.spg <- function(rep, dir=".", nfolds=10, r=25)
{
   oldwd <- getwd()
   setwd(dir)

   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep), header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep), header=FALSE))

   system(sprintf("%s/spg.sh %d %d %d", oldwd, rep, nfolds, r))
   B.spg <- as.matrix(read.table(sprintf("spg_beta_%s.txt", rep), header=FALSE))
   P.spg <- Xtest %*% B.spg
   R2.spg <- R2(as.numeric(P.spg), as.numeric(Ytest))
   cat("R2.spg:", R2.spg, "\n")

   B.spg.l <- as.matrix(read.table(sprintf("spg_beta_lasso_%s.txt", rep),
	 header=FALSE))
   P.spg.l <- Xtest %*% B.spg.l
   R2.spg.l <- R2(as.numeric(P.spg.l), as.numeric(Ytest))
   cat("R2.spg.l:", R2.spg.l, "\n")

   setwd(oldwd)

   list(
      R2=c(spg=R2.spg, lasso=R2.spg.l),
      beta=list(spg=B.spg, lasso=B.spg.l)
   )
}

run.lasso <- function(rep, dir=".", nfolds=10, r=25)
{
   oldwd <- getwd()
   setwd(dir)

   Xtrain <- as.matrix(read.table(sprintf("Xtrain_%s.txt", rep),
	    header=FALSE))
   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep), header=FALSE))
   Ytrain <- as.matrix(read.table(sprintf("Ytrain_%s.txt", rep),
	    header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep), header=FALSE))
   XtrainB <- as.matrix(read.table(sprintf("XtrainB_%s.txt", rep),
	    header=FALSE))
   XtestB <- as.matrix(read.table(sprintf("XtestB_%s.txt", rep),
	    header=FALSE))

   ytrain <- as.numeric(Ytrain)
   ytest <- as.numeric(Ytest)

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   L1mult <- seq(1, 1e-9, length=r) # a multiplier on maxL1, not penalty itself
   maxL1B <- maxlambda1(XtrainB, ytrain)
   r.gl <- sapply(L1mult, function(m) {
      crossval(XtrainB, ytrain, fun=lasso3, lambda1=m * maxL1B,
	    nfolds=nfolds)
   })
   m.gl <- which(r.gl == max(r.gl, na.rm=TRUE))
   gl <- lasso3(XtrainB, ytrain, lambda1=L1mult[m.gl] * maxL1B)
   P.gl <- XtestB %*% gl
   R2.gl <- R2(as.numeric(P.gl), ytest)
   cat("R2.gl:", R2.gl, "\n")

   setwd(oldwd)

   list(
      R2=R2.gl,
      beta=matrix(gl, p, K)
   )
}

run.ridge <- function(rep, dir=".", nfolds=10, r=25)
{
   oldwd <- getwd()
   setwd(dir)

   Xtrain <- as.matrix(read.table(sprintf("Xtrain_%s.txt", rep),
	    header=FALSE))
   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep), header=FALSE))
   Ytrain <- as.matrix(read.table(sprintf("Ytrain_%s.txt", rep),
	    header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep), header=FALSE))
   XtrainB <- as.matrix(read.table(sprintf("XtrainB_%s.txt", rep),
	    header=FALSE))
   XtestB <- as.matrix(read.table(sprintf("XtestB_%s.txt", rep),
	    header=FALSE))

   ytrain <- as.numeric(Ytrain)
   ytest <- as.numeric(Ytest)

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)
   
   L2 <- seq(0, 3, length=r)
   r <- sapply(L2, function(l2) {
      crossval(XtrainB, ytrain, fun=ridge, lambda2=l2, nfolds=nfolds)
   })
   m <- which(r == max(r, na.rm=TRUE))
   l <- ridge(XtrainB, ytrain, lambda2=L2[m])
   P <- XtestB %*% l
   R2.r <- R2(as.numeric(P), as.numeric(ytest))
   cat("R2.r:", R2.r, "\n")

   setwd(oldwd)

   list(
      R2=R2.r,
      beta=matrix(l, p, K)
   )
}

run.groupridge <- function(rep, dir=".", nfolds=10, r=25)
{
   oldwd <- getwd()
   setwd(dir)

   Xtrain <- as.matrix(read.table(sprintf("Xtrain_%s.txt", rep),
	    header=FALSE))
   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep), header=FALSE))
   Ytrain <- as.matrix(read.table(sprintf("Ytrain_%s.txt", rep), header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep), header=FALSE))
   grp <- scan(sprintf("grp_%s.txt", rep))

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   # different L1 on rows, different L3 on columns
   L1mult <- seq(1, 1e-6, length=r) # a multiplier on maxL1, not penalty itself
   maxL1 <- maxlambda1(Xtrain, Ytrain)
   L3 <- seq(0, 3, length=r)
   r <- sapply(L1mult, function(m) {
      cat(m, "; ")
      s <- sapply(L3, function(l3) {
	 cat(l3, " ")
	 crossval(Xtrain, Ytrain, fun=groupridge3,
	       lambda1=maxL1 * m, lambda3=l3, grp=grp, nfolds=nfolds)
      })
      cat("\n")
      s
   })
   cat("crossval done\n")
   m <- which(r == max(r, na.rm=TRUE), arr.ind=TRUE)
   cat("Best/Worst crossval:", max(r), min(r), "; best at",
	 maxL1 * L1mult[m[2]], L3[m[1]], "\n")
   gr <- groupridge3(Xtrain, Ytrain,
	 lambda1=maxL1 * L1mult[m[2]], lambda3=L3[m[1]], grp=grp)
   R2.gr <- R2(as.numeric(Xtest %*% gr), as.numeric(Ytest))

   setwd(oldwd)

   list(
      R2=R2.gr,
      beta=gr
   )
}

nreps <- 3
dir <- "Expr1"
if(!file.exists(dir)) {
   dir.create(dir)
   sapply(1:nreps, makedata, dir=dir, N=100, p=50, K=10)
}

r.gr <- lapply(1:nreps, run.groupridge, dir=dir, r=3, nfolds=3)
r.lasso <- lapply(1:nreps, run.lasso, dir=dir, r=3, nfolds=3)
r.spg <- lapply(1:nreps, run.spg, dir=dir, r=3, nfolds=3)
r.ridge <- lapply(1:nreps, run.ridge, dir=dir, r=3, nfolds=3)

# Measure recovery of non-zeros
recovery <- function(obj, dir)
{
   lapply(seq(along=obj), function(rep) {
      beta <- as.matrix(read.table(sprintf("%s/B_%s.txt", dir, rep)))
      pred <- prediction(predictions=obj[[rep]]$beta, labels=beta != 0)
      roc <- performance(pred, "sens", "spec")
      prc <- performance(pred, "prec", "rec")
      list(roc=roc, prc=prc)
   })
}

rec.gr <- recovery(r.gr, dir)
rec.lasso <- recovery(r.lasso, dir)
rec.spg <- recovery(r.spg, dir)
rec.ridge <- recovery(r.ridge, dir)

stop()

r <- do.call(rbind, r)
t.test(r[,1], r[,2])
colMeans(r)

save.image(r, file="test_groupridge.RData")


stop()


##system.time({
##   g1 <- groupridge(X, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
##      eps=1e-4, verbose=TRUE)
##})
##(err1 <- mean((X %*% g1 - Y)^2))
##system.time({
##   g2 <- groupridge(X, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
##      eps=1e-10, verbose=TRUE)
##})
##(err2 <- mean((X %*% g2 - Y)^2))
##(err3 <- mean((X %*% g3 - Y)^2))
##g4 <- lasso(X, drop(Y), lambda1=1e-2, maxepoch=100)
##(err4 <- mean((X %*% g4 - Y)^2))
##g5 <- lasso.activeset(X, drop(Y), lambda1=1e-2)
##(err5 <- mean((X %*% g5 - Y)^2))
#system.time({
#   g <- glmnet(X, drop(Y), lambda=1e-2)
#})
#b <- as.matrix(coef(g))[-1, ]
##g6 <- groupridge_simple(X, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
##      eps=1e-16,  verbose=TRUE)
##(err6 <- mean((X %*% g6 - Y)^2))
##g7 <- simpler(X, drop(Y), lambda=1e-2)
##(err7 <- mean((X %*% g7 - Y)^2))
##system.time({
##   g8 <- groupridge2(X, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
##      eps=1e-4, verbose=TRUE)
##})
#
#system.time({
#   g2 <- groupridge2(X, Y, lambda1=1e-2, maxiter=1e4, verbose=FALSE)
#})
#system.time({
#   g3 <- groupridge3(X, Y, lambda1=1e-2, maxiter=1e4, verbose=FALSE)
#})
#system.time({
#   l2 <- lasso2(X, drop(Y), lambda1=1e-2, maxiter=1e4, verbose=FALSE)
#})
#system.time({
#   l3 <- lasso3(X, drop(Y), lambda1=1e-2, maxiter=1e4, verbose=FALSE)
#})

#r <- cbind(
#   True=drop(B),
#   glmnet=b,
#   #gr1=drop(g1),
#   #gr2=drop(g2),
#   gr2=drop(g2),
#   gr3=drop(g3),
#   #l1=g4,
#   #l2=g5,
#   #simple=drop(g6),
#   #simpler=drop(g7),
#   #gr8=drop(g8))
#   l2=drop(l2),
#   l3=drop(l3)
#)
#cor(r)
#cor(sign(r))

#Y1 <- matrix(rep(Y, ncol(r)), nrow=N)
#Y2 <- X %*% r
#mse <- apply(Y2 - Y1, 2, function(x) sum(x^2))
#
##Xs <- scale(X0)
##system.time({
##   g31 <- groupridge2(X, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
##      eps=1e-4, verbose=TRUE)
##})
##system.time({
##   g32 <- groupridge2(Xs, Y, lambda1=1e-2, lambda2=0, lambda3=0, g=grp,
##      eps=1e-4, verbose=TRUE)
##})


#################################################################################
## Test groupridge in multi-task setting, lasso only
#N <- 100
#p <- 1000
#K <- 10
#
## scale to zero mean and unit norm (not unit variance)
#X0 <- matrix(rnorm(N * p), N, p)
#X1 <- sweep(X0, 2, colMeans(X0))
#X <- sweep(X1, 2, sqrt(diag(crossprod(X1))), FUN="/")
##B <- matrix(rnorm(p * K), p, K)
##B[sample(p, p - 10),] <- 0
#b <- rnorm(p) * sample(0:1, p, TRUE, prob=c(0.9, 0.1))
#B <- sapply(1:K, function(k) b)
#
#Y <- X %*% B
#
#g2 <- groupridge3(X, Y[,1], lambda1=1e-2, maxiter=1e6)
#g3 <- groupridge3(X, Y, lambda1=1e-2, maxiter=1e6)
#l3 <- lasso3(X, Y[,1], lambda1=1e-2, maxiter=1e6)
#cor(cbind(g2, g3[,1], l3))
#cor(g3)
#
#stop()

################################################################################
# Test groupridge in multi-task setting, lasso + ridge + group ridge
# B is same for all tasks
N <- 100
p <- 1000
K <- 10

# scale to zero mean and unit norm (not unit variance)
X0 <- matrix(rnorm(N * p), N, p)
X1 <- sweep(X0, 2, colMeans(X0))
X <- sweep(X1, 2, sqrt(diag(crossprod(X1))), FUN="/")
#B <- matrix(rnorm(p * K), p, K)
#B[sample(p, p - 10),] <- 0
b <- rnorm(p) * sample(0:1, p, TRUE, prob=c(0.9, 0.1))
B <- sapply(1:K, function(k) b)
Y <- X %*% B

lambda2 <- 1e-2

# Test multi-task, ridge only 
r <- qr.solve(crossprod(X) + diag(p) * lambda2, crossprod(X, Y))
g1 <- groupridge3(X, Y, lambda1=0, lambda2=lambda2, maxiter=1e6)
diag(cor(r, g1))
mean(( X %*% r - Y)^2)
mean(( X %*% g1 - Y)^2)

# Test mixed lasso/ridge
lambda1 <- 1e-2
g2 <- groupridge3(X, Y, lambda1=lambda1, lambda2=1e-3, maxiter=1e6)
g3 <- groupridge3(X, Y, lambda1=lambda1, lambda2=1e0, maxiter=1e6)
g4 <- groupridge3(X, Y, lambda1=lambda1, lambda2=1e2, maxiter=1e6)
g4 <- groupridge3(X, Y, lambda1=lambda1, lambda2=1e4, maxiter=1e6)
cor(cbind(g2, g3, g4))

# Test group ridge, shouldn't make a difference here since B
# is same for all tasks
g2 <- groupridge3(X, Y,
      lambda1=lambda1, lambda2=1e-3, lambda3=1e-3, maxiter=1e6)
g3 <- groupridge3(X, Y,
      lambda1=lambda1, lambda2=1e-3, lambda3=1e0, maxiter=1e6)
g4 <- groupridge3(X, Y,
      lambda1=lambda1, lambda2=1e-3, lambda3=1e3, maxiter=1e6)
diag(cor(cbind(g2, g3, g4)))

################################################################################
# Test groupridge in multi-task setting, lasso + ridge + group ridge
# B weights differ between tasks, to see effect of group ridge
N <- 100
p <- 200
K <- 10

# scale to zero mean and unit norm (not unit variance)
X0 <- matrix(rnorm(N * p), N, p)
X1 <- sweep(X0, 2, colMeans(X0))
X <- sweep(X1, 2, sqrt(diag(crossprod(X1))), FUN="/")
B <- matrix(rnorm(p * K), p, K)
B <- sapply(1:K, function(k) {
   rnorm(p) * sample(0:1, p, TRUE, prob=c(0.9, 0.1))
})
Y <- X %*% B
grp <- rep(1, K)

g2 <- groupridge3(X, Y,
      lambda1=0, lambda2=0, lambda3=1e-3, maxiter=1e6, grp=grp)
g3 <- groupridge3(X, Y,
      lambda1=0, lambda2=0, lambda3=1e0, maxiter=1e6, grp=grp)
g4 <- groupridge3(X, Y,
      lambda1=0, lambda2=0, lambda3=1e2, maxiter=1e6, grp=grp)
cor(cbind(g2[,1], g3[,1], g4[,1]))
sum((g2[,1] - g2[,2])^2)
sum((g3[,1] - g3[,2])^2)
sum((g4[,1] - g4[,2])^2)

