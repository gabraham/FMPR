# Test groupridge

#library(MASS)
library(glmnet)
library(gdata)
library(ROCR)
library(ggplot2)

#dyn.load("groupridge.so")

options(error=dump.frames)

s <- sample(1e6L, 1)
set.seed(s)

source("eval.R")
source("methods.R")

#groupridge <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g,
#      maxiter=1e5, eps=1e-6, verbose=FALSE)
#{
#   p <- ncol(X)
#   K <- ncol(Y)
#
#   if(length(lambda1) == 1)
#      lambda1 <- rep(lambda1, K)
#   
#   if(length(lambda2) == 1)
#      lambda2 <- rep(lambda2, K)
#
#   if(length(lambda3) == 1)
#      lambda3 <- rep(lambda3, K)
#
#   r <- .C("groupridge", as.numeric(X), as.numeric(Y), 
#      numeric(p * K), nrow(X), ncol(X), ncol(Y),
#      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
#      as.integer(g), as.integer(maxiter), as.double(eps),
#      as.integer(verbose)
#   )
#   matrix(r[[3]], p, K)
#}
#
#groupridge_simple <- function(X, Y, lambda1=0, lambda2=0, lambda3=0,
#      g=NULL, maxiter=1e3, eps=1e-6, verbose=FALSE)
#{
#   p <- ncol(X)
#   K <- ncol(Y)
#
#   if(length(lambda1) == 1)
#      lambda1 <- rep(lambda1, K)
#   
#   if(length(lambda2) == 1)
#      lambda2 <- rep(lambda2, K)
#
#   if(length(lambda3) == 1)
#      lambda3 <- rep(lambda3, K)
#
#   if(is.null(g))
#      g <- 1:K
#
#   r <- .C("groupridge_simple", as.numeric(X), as.numeric(Y), 
#      numeric(p * K), nrow(X), ncol(X), ncol(Y),
#      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
#      as.integer(g), as.integer(maxiter), as.double(eps),
#      as.integer(verbose)
#   )
#   matrix(r[[3]], p, K)
#}
#
#
#groupridge2 <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g=NULL,
#      maxiter=1e5, eps=1e-6, verbose=TRUE)
#{
#   p <- ncol(X)
#   K <- ncol(Y)
#
#   if(length(lambda1) == 1)
#      lambda1 <- rep(lambda1, K)
#   
#   if(length(lambda2) == 1)
#      lambda2 <- rep(lambda2, K)
#
#   if(length(lambda3) == 1)
#      lambda3 <- rep(lambda3, K)
#      
#   if(is.null(g))
#      g <- 1:K
#
#   r <- .C("groupridge2", as.numeric(X), as.numeric(Y), 
#      numeric(p * K), nrow(X), ncol(X), ncol(Y),
#      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
#      as.integer(g), as.integer(maxiter),
#      as.double(eps), as.integer(verbose))
#   matrix(r[[3]], p, K)
#}
#
## lambda1: scalar or K-vector
## lambda2: scalar or K-vector
## lambda3: scalar 
#groupridge3 <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, grp=NULL,
#      maxiter=1e5, eps=1e-6, verbose=FALSE)
#{
#   p <- ncol(X)
#   Y <- cbind(Y)
#   K <- ncol(Y)
#
#   if(length(lambda1) == 1)
#      lambda1 <- rep(lambda1, K)
#   
#   if(length(lambda2) == 1)
#      lambda2 <- rep(lambda2, K)
#
#   if(is.null(grp))
#      grp <- 1:K
#
#   r <- .C("groupridge3", as.numeric(X), as.numeric(Y), 
#      numeric(p * K), nrow(X), ncol(X), K,
#      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
#      as.integer(grp), as.integer(maxiter),
#      as.double(eps), as.integer(verbose), integer(1))
#   matrix(r[[3]], p, K)
#}
#
#lasso2 <- function(X, y, lambda1=0,
#      maxiter=1e5, eps=1e-6, verbose=TRUE)
#{
#   p <- ncol(X)
#
#   if(length(lambda1) == 1)
#      lambda1 <- rep(lambda1, K)
#   
#   r <- .C("lasso2", as.numeric(X), as.numeric(y), numeric(p), 
#      nrow(X), ncol(X), as.numeric(lambda1),
#      as.integer(maxiter), as.double(eps),
#      as.integer(verbose)
#   )
#   matrix(r[[3]], p)
#}
#
#lasso3 <- function(X, y, lambda1=0,
#      maxiter=1e5, eps=1e-6, verbose=FALSE)
#{
#   p <- ncol(X)
#
#   r <- .C("lasso3", as.numeric(X), as.numeric(y), numeric(p), 
#      nrow(X), ncol(X), as.numeric(lambda1),
#      as.integer(maxiter), as.double(eps),
#      as.integer(verbose)
#   )
#   matrix(r[[3]], p)
#}
#
#ridge <- function(X, Y, lambda2=0)
#{
#   XX <- crossprod(X)
#   XY <- crossprod(X, Y)
#   p <- ncol(X)
#   qr.solve(XX + diag(p) * lambda2, XY)
#}
#
## a very simple version of lasso, for sanity checking
#simpler <- function(x, y, lambda, maxiter=2000)
#{
#   N <- length(y)
#   p <- ncol(x)
#   beta <- numeric(p)
#   v <- diag(crossprod(x))
#   for(iter in 1:maxiter)
#   {
#      for(j in 1:p)
#      {
#	 d1 <- crossprod(x[,j], x %*% beta - y)
#	 d2 <- v[j]
#	 d <- beta[j] - d1 / d2 - lambda * sign(beta[j])
#
#	 # The solution is zero when:
#	 # 1. beta[j] starts as zero, and Newton step is 
#	 #    less than lambda
#	 # or
#	 # 2. beta[j] doesn't start as zero, but changes sign
#	 if(beta[j] == 0 && abs(d1 / d2) <= lambda) {
#	    beta[j] <- 0
#	 } else if(d * beta[j] < 0) {
#	    beta[j] <- 0
#	 } else {
#	    beta[j] <- d
#	 }
#      }
#   }
#   beta
#}
#
## zero mean, unit norm (not unit variance)
#standardise <- function(x)
#{
#   x1 <- sweep(x, 2, colMeans(x))
#   v <- apply(x1, 2, function(z) sqrt(sum(z^2)))
#   sweep(x1, 2, v, FUN="/")
#}
#
##source("lasso.R")
#
#maxlambda1 <- function(X, Y)
#{
#   Y <- cbind(Y)
#   r <- .C("maxlambda1", as.numeric(X), as.numeric(Y),
#      numeric(ncol(Y)), nrow(X), ncol(X), ncol(Y))
#   r[[3]]
#}
#
#blockX <- function(X, p, C)
#{
#   N <- nrow(X)
#   
#   Xblock <- matrix(0, N * C, p * C) 
#   s1 <- seq(1, N * C, N)
#   s2 <- seq(1, p * C, p)
#
#   for(i in 1:C)
#   {
#      j <- s1[i]:(s1[i] + N - 1)
#      k <- s2[i]:(s2[i] + p - 1)
#      Xblock[j, k] <- X
#   }
#
#   Xblock
#}
#
#R2 <- function(pr, y) 1 - sum((pr - y)^2) / sum((y - mean(y))^2)

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


##################################################################################
### Test groupridge in multi-task setting, lasso only
#N <- 100
#p <- 50
#K <- 5
##
### scale to zero mean and unit norm (not unit variance)
##X0 <- matrix(rnorm(N * p), N, p)
##X1 <- sweep(X0, 2, colMeans(X0))
##X <- sweep(X1, 2, sqrt(diag(crossprod(X1))), FUN="/")
###B <- matrix(rnorm(p * K), p, K)
###B[sample(p, p - 10),] <- 0
##b <- rnorm(p) * sample(0:1, p, TRUE, prob=c(0.9, 0.1))
##B <- sapply(1:K, function(k) b)
##
##Y <- X %*% B
##
#X <- standardise(matrix(rnorm(N * p), N, p))
#Y <- scale(matrix(rnorm(N * K), N, K))
#g4 <- groupridge4(X, Y, lambda1=1e-2, maxiter=1e5, verbose=FALSE)
#stop()
##g3 <- groupridge3(X, Y, lambda1=1e-2, maxiter=1e6)
##l3 <- lasso3(X, Y[,1], lambda1=1e-2, maxiter=1e6)
##cor(cbind(g2, g3[,1], l3))
##cor(g3)
##
#stop()
#
#################################################################################
## Test groupridge in multi-task setting, lasso + ridge + group ridge
## B is same for all tasks
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
#Y <- X %*% B
#
#lambda2 <- 1e-2
#
## Test multi-task, ridge only 
#r <- qr.solve(crossprod(X) + diag(p) * lambda2, crossprod(X, Y))
#g1 <- groupridge3(X, Y, lambda1=0, lambda2=lambda2, maxiter=1e6)
#diag(cor(r, g1))
#mean(( X %*% r - Y)^2)
#mean(( X %*% g1 - Y)^2)
#
## Test mixed lasso/ridge
#lambda1 <- 1e-2
#g2 <- groupridge3(X, Y, lambda1=lambda1, lambda2=1e-3, maxiter=1e6)
#g3 <- groupridge3(X, Y, lambda1=lambda1, lambda2=1e0, maxiter=1e6)
#g4 <- groupridge3(X, Y, lambda1=lambda1, lambda2=1e2, maxiter=1e6)
#g4 <- groupridge3(X, Y, lambda1=lambda1, lambda2=1e4, maxiter=1e6)
#cor(cbind(g2, g3, g4))
#
## Test group ridge, shouldn't make a difference here since B
## is same for all tasks
#g2 <- groupridge3(X, Y,
#      lambda1=lambda1, lambda2=1e-3, lambda3=1e-3, maxiter=1e6)
#g3 <- groupridge3(X, Y,
#      lambda1=lambda1, lambda2=1e-3, lambda3=1e0, maxiter=1e6)
#g4 <- groupridge3(X, Y,
#      lambda1=lambda1, lambda2=1e-3, lambda3=1e3, maxiter=1e6)
#diag(cor(cbind(g2, g3, g4)))
#

################################################################################
# Test groupridge in multi-task setting, lasso + ridge + group ridge
# B weights differ between tasks, to see effect of group ridge
run <- function()
{
   N <- 100
   p <- 200
   K <- 10
   
   # scale to zero mean and unit norm (not unit variance)
   X0 <- matrix(rnorm(N * p), N, p)
   X <- standardise(X0)
   B <- 0
   while(all(B == 0)) {
      b <- rnorm(p) * sample(0:1, p, TRUE, prob=c(0.7, 0.3))
      B <- sapply(1:K, function(k) {
         if(k < 5) {
            b
         } else {
            rev(b)
         }
      })
   }
   Y <- scale(X %*% B + rnorm(N * K, 0, 0.1))
   grp <- rep(1, K)
   
   Xb <- blockX(X, p, K)
   
   G1 <- matrix(1, K, K)
   
   Rc2 <- cor(Y)
   G2 <- sign(Rc2) * (abs(Rc2) > 0.3)
   
   r1 <- optim.groupridge(X=X, Y=Y, G=G2, nfolds=10, maxiter=1e4)
   r2 <- optim.lasso(X=Xb, Y=as.numeric(Y), nfolds=10, maxiter=1e4)
   
   g1 <- groupridge4(X=X, Y=Y, lambda1=r1$opt[1] * 0.1, lambda2=r1$opt[2],
         lambda3=r1$opt[3], maxiter=1e4, G=G2)
   g2 <- matrix(lasso3(X=Xb, y=as.numeric(Y), lambda1=r2$opt[1]), p, K)
   
   auc(as.numeric(B != 0), as.numeric(abs(g1)))
   auc(as.numeric(B != 0), as.numeric(abs(g2)))
   
   pred <- prediction(
      predictions=cbind(as.numeric(abs(g1)), as.numeric(abs(g2))),
      labels=cbind(as.numeric(B != 0), as.numeric(B != 0))
   )
   perf.roc <- performance(pred, "sens", "spec")
   perf.prc <- performance(pred, "prec", "rec")
   
   list(roc=perf.roc, prc=perf.roc)
}

res <- lapply(1:1, function(i) run())

meth <- c("groupridge", "lasso")

resd <- lapply(res, function(r) {
   spec <- do.call(cbind, r$roc@x.values)
   sens <- do.call(cbind, r$roc@y.values)
   prec <- do.call(cbind, r$prc@x.values)
   rec <- do.call(cbind, r$prc@y.values)
   data.frame(
      Specificity=as.numeric(spec),
      Sensitivity=as.numeric(sens),
      Precision=as.numeric(prec),
      Recall=as.numeric(rec),
      Method=rep(meth, each=nrow(spec))
   )
})
resd2 <- do.call(rbind, resd)

gg1 <- ggplot(resd2, aes(x=Specificity, y=Sensitivity, colour=Method))
gg1 <- gg1 + geom_line(size=2)

gg2 <- ggplot(resd2, aes(x=Recall, y=Precision, colour=Method))
gg2 <- gg2 + geom_line(size=2)

pdf("groupridge_block.pdf", width=12)
print(gg1)
print(gg2)
dev.off()

