# Test groupridge

options(error=dump.frames)

library(FMPR)

set.seed(1898)

N <- 100
p <- 10
K <- 3


B <- getB(p=p, K=K, w=0.5, type="same")
d <- makedata(N=N, K=K, B=B, p=p, save=FALSE, sigma=1, rep=1)

R <- cor(d$Ytrain)
diag(R) <- 0
Rthresh <- 0.1

Gt <- sign(R) * (abs(R) > Rthresh)

Xtrain <- scale(d$Xtrain)
Ytrain <- scale(d$Ytrain)

write.table(Xtrain, file="Xtrain.txt", col.names=FALSE, row.names=FALSE)
write.table(Ytrain, file="Ytrain.txt", col.names=FALSE, row.names=FALSE)

#l <- max(maxlambda1(Xtrain, Ytrain))
#
#lambda1 <- l * seq(1, 1e-3, length=ngrid)
#lambda2 <- 10^seq(-3, 6, length=ngrid)
#lambda3 <- 10^seq(-3, 6, length=ngrid)
#
#lambda <- 10^seq(-3, 6, length=ngrid)
#gamma <- 10^seq(-3, 6, length=ngrid)

C <- gennetwork(Ytrain, threshold=Rthresh, weight.fun=sqr)
z.s.t <- spg(Xtrain, Ytrain, C, lambda=1.371, gamma=1.21, maxiter=1e4,
      verbose=TRUE)

#r.s.t <- optim.spg(X=Xtrain, Y=Ytrain, nfolds=nfolds,
#      lambda=lambda, gamma=gamma, C=Gt, verbose=TRUE)

#
#z.t <- fmpr(X=Xtrain, Y=Ytrain, G=Gt,
#   lambda1=1e-3, lambda2=1e-3, lambda3=1e-3,
#   type="threshold", simplify=TRUE)
#
#z.w <- fmpr(X=Xtrain, Y=Ytrain, G=Gt,
#   lambda1=1e-3, lambda2=1e-3, lambda3=1e-3,
#   type="weighted", simplify=TRUE)
#
#sum((z.w-z.t)^2)
#
#
#m.s.t <- spg(X=Xtrain, Y=Ytrain,
#   lambda=r.s.t$opt[1], gamma=r.s.t$opt[2], C=Gt, simplify=TRUE)
#
#r.t <- optim.fmpr(X=Xtrain, Y=Ytrain, nfolds=nfolds,
#      lambda1=lambda1, lambda2=lambda2, lambda3=lambda3,
#      G=Gt, type="threshold")
#
#r.w1 <- optim.fmpr(X=Xtrain, Y=Ytrain, nfolds=nfolds,
#      lambda1=lambda1, lambda2=lambda2, lambda3=lambda3,
#      G=Gw1, type="weighted")
#
#r.w2 <- optim.fmpr(X=Xtrain, Y=Ytrain, nfolds=nfolds,
#      lambda1=lambda1, lambda2=lambda2, lambda3=lambda3,
#      G=Gw2, type="weighted")
#
#m.t <- fmpr(X=Xtrain, Y=Ytrain,
#   lambda1=r.t$opt[1],
#   lambda2=r.t$opt[2],
#   lambda3=r.t$opt[3],
#   G=Gt, type="threshold", simplify=TRUE)
#
#m.w1 <- fmpr(X=Xtrain, Y=Ytrain,
#   lambda1=r.w1$opt[1],
#   lambda2=r.w1$opt[2],
#   lambda3=r.w1$opt[3],
#   G=Gw1, type="weighted", simplify=TRUE)
#
#m.w2 <- fmpr(X=Xtrain, Y=Ytrain,
#   lambda1=r.w2$opt[1],
#   lambda2=r.w2$opt[2],
#   lambda3=r.w2$opt[3],
#   G=Gw2, type="weighted", simplify=TRUE)
#
#
#P.t <- Xtest %*% m.t
#P.w1 <- Xtest %*% m.w1
#P.w2 <- Xtest %*% m.w2
#P.s.t <- Xtest %*% m.s.t
#
#R2 <- function(p, y)
#{
#   1 - sum((p - y)^2) / sum((y - mean(y))^2)
#}
#
#R2.t <- R2(as.numeric(P.t), as.numeric(Ytest))
#R2.w1 <- R2(as.numeric(P.w1), as.numeric(Ytest))
#R2.w2 <- R2(as.numeric(P.w2), as.numeric(Ytest))
#R2.s.t <- R2(as.numeric(P.s.t), as.numeric(Ytest))
#
#c(R2.t, R2.w1, R2.w2, R2.s.t)
#
