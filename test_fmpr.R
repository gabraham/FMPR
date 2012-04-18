options(error=dump.frames)

library(FMPR)

#set.seed(18281)

source("FMPR/R/methods.R")

#N <- 100
#p <- 200
#K <- 10

#X <- scale(matrix(rnorm(N * p), N, p))
#b <- 0
#while(all(b == 0)) {
#   b <- rnorm(p) * sample(0:1, size=p, replace=TRUE, prob=c(0.8, 0.2))
#}
#B0 <- matrix(b, p, K)
##B <- B0 * rnorm(p * K, 0, 0.3)
#B <- B0
#Y <- scale(X %*% B + rnorm(N * K, 0, 2))
X <- scale(as.matrix(read.table("Expr13/Xtrain_1.txt")))
Y <- scale(as.matrix(read.table("Expr13/Ytrain_1.txt")))
B <- as.matrix(read.table("Expr13/B_1.txt"))

l <- max(maxlambda1(X, Y))

ngrid <- 50
maxiter <- 1e4

lambda1 <- l * 2^seq(-10, 0, length=ngrid)
lambda2 <- 0
lambda3 <- 1e4
#lambda1 <- 10^seq(-3, 6, length=ngrid)
#lambda2 <- 10^seq(-3, 6, length=ngrid)
#lambda3 <- 10^seq(-3, 6, length=ngrid)
#lambda <- 10^seq(-3, 6, length=ngrid)
gamma <- 10^seq(-3, 6, length=ngrid)

#lambda3 <- numeric(ngrid)
#gamma <- numeric(ngrid)

nfolds <- 5

R <- cor(Y)
G <- graph.sqr(R)
C <- gennetwork(Y, cortype=2, corthresh=0)

r1 <- optim.fmpr(X=X, Y=Y, lambda1=lambda1, lambda2=0, lambda3=lambda3,
      sparse=FALSE)
r2 <- optim.spg(X=X, Y=Y, lambda=lambda1, gamma=gamma, cortype=2, corthresh=0)

f1 <- fmpr(X, Y, G=G, lambda1=r1$opt[1], lambda2=r1$opt[2], lambda3=r1$opt[3],
   sparse=FALSE, maxiter=maxiter, eps=1e-6, verbose=FALSE, simplify=TRUE)
f2 <- spg(X, Y, C=C, lambda=r2$opt[1], gamma=r2$opt[2], tol=1e-6)[[1]][[1]]

#pr1 <- prediction(
#   predictions=lapply(f1, function(x) abs(x[[1]][[1]])),
#   labels=lapply(1:ngrid, function(i) as.numeric(B != 0))
#)
#pr2 <- prediction(
#   predictions=f2,
#   labels=lapply(1:ngrid, function(i) as.numeric(B != 0))
#)
#pr3 <- prediction(
#   predictions=lapply(f3, function(x) abs(x[[1]])),
#   labels=lapply(1:ngrid, function(i) as.numeric(B != 0))
#)

pr1 <- prediction(predictions=as.numeric(abs(f1)), labels=as.numeric(B != 0))
pr2 <- prediction(predictions=as.numeric(abs(f2)), labels=as.numeric(B != 0))

par(mfrow=c(1, 2))
plot(roc1 <- performance(pr1, "sens", "spec"), col=1)
plot(roc2 <- performance(pr2, "sens", "spec"), col=3, add=TRUE)
plot(prc1 <- performance(pr1, "prec", "rec"), col=1)
plot(prc2 <- performance(pr2, "prec", "rec"), col=3, add=TRUE)

#r1 <- optim.fmpr(X=X, Y=Y,
#   lambda1=lambda1,
#   lambda2=0,
#   lambda3=lambda3,
#   maxiter=maxiter
#)
#
#r2 <- optim.spg(X=X, Y=Y,
#   lambda=lambda,
#   gamma=gamma,
#   maxiter=maxiter
#)
#
#r3 <- optim.ridge(X, Y, lambda2=lambda2)
#r4 <- optim.lasso(blockX(X, p, K), as.numeric(Y), lambda1=lambda1)
#
#
#f1 <- fmpr(X, Y, G=G, lambda1=r1$opt[1], lambda2=r1$opt[2], lambda3=r1$opt[3],
#      simplify=TRUE, sparse=FALSE, maxiter=maxiter)
#f2 <- spg(X, Y, C=C, lambda=r2$opt[1], gamma=r2$opt[2], simplify=TRUE,
#      maxiter=maxiter)
#f3 <- ridge(X, Y, lambda2=r3$opt[1])[[1]]
#f4 <- lasso(blockX(X, p, K), as.numeric(Y), lambda1=r4$opt[1])
#
#pred1 <- prediction(
#   predictions=as.numeric(abs(f1)), labels=as.numeric(B != 0)
#)
#pred2 <- prediction(
#   predictions=as.numeric(abs(f2)), labels=as.numeric(B != 0)
#)
#pred3 <- prediction(
#   predictions=as.numeric(abs(f3)), labels=as.numeric(B != 0)
#)
#pred4 <- prediction(
#   predictions=as.numeric(abs(f4)), labels=as.numeric(B != 0)
#)
#
#prc1 <- performance(pred1, "prec", "rec")
#prc2 <- performance(pred2, "prec", "rec")
#prc3 <- performance(pred3, "prec", "rec")
#prc4 <- performance(pred4, "prec", "rec")
#
#roc1 <- performance(pred1, "sens", "spec")
#roc2 <- performance(pred2, "sens", "spec")
#roc3 <- performance(pred3, "sens", "spec")
#roc4 <- performance(pred4, "sens", "spec")
#
#par(mfrow=c(1, 2))
#
#plot(roc1, col=1)
#plot(roc2, add=TRUE, col=2)
#plot(roc3, add=TRUE, col=3)
#plot(roc4, add=TRUE, col=4)
#abline(1, -1, lty=2)
#
#plot(prc1, col=1, ylim=c(0, 1))
#plot(prc2, col=2, add=TRUE)
#plot(prc3, col=3, add=TRUE)
#plot(prc4, col=4, add=TRUE)
#
#
