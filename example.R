### R code from vignette source 'fmpr.rnw'

library(FMPR)
library(doMC)
registerDoMC(cores=4)



N <- 100
p <- 100
K <- 10 
w <- 0.1
B <- getB(p=p, K=K, w=w, type="same", sparsity=0.8)
print(table(sign(B)))

X <- scale(matrix(rnorm(N * p), N, p))
Y <- scale(X %*% B + rnorm(N * K))


G <- graph.sqr(cor(Y))
C <- gennetwork(Y, corthresh=0, cortype=2)


l <- max(maxlambda1(X, Y))
ngrid <- 25
lambda <- 2^seq(-0.01, -10, length=ngrid) * l
gamma <- c(0, 10^seq(-5, 5, length=ngrid))
nfolds <- 5


# L2 fusion penalty
system.time({
   opt.f <- optim.fmpr(X=X, Y=Y, graph.fun=graph.sqr,
      lambda=lambda, gamma=gamma, nfolds=nfolds)
})

# L1 fusion penalty
system.time({
   opt.s <- optim.spg(X=X, Y=Y, cortype=2, corthresh=0,
      lambda=lambda, gamma=gamma, nfolds=nfolds)
})


f <- fmpr(X=X, Y=Y, G=G, lambda=opt.f$opt["lambda"],
   gamma=opt.f$opt["gamma"], simplify=TRUE)
s <- spg(X=X, Y=Y, C=C, lambda=opt.s$opt["lambda"],
   gamma=opt.s$opt["gamma"], simplify=TRUE)

# Lasso, special case of FMPR
w <- which(opt.f$R2[, 1] == max(opt.f$R2[, 1]))
l <- fmpr(X=X, Y=Y, lambda=lambda[w], gamma=0, simplify=TRUE)


measures <- list(ROC=c("sens", "spec"), PRC=c("prec", "rec"))
res <- lapply(list(FMPR=f, GFlasso=s, Lasso=l), function(B2) {
   lapply(measures, function(m) {
      performance(prediction(
         labels=as.numeric(B != 0),
         predictions=as.numeric(abs(B2))
         ), m[1], m[2])
   })
})


pdf("example.pdf", width=14)

par(mfrow=c(1, 2), mar=c(4, 4, 2, 2) + 0.1)

plot(res$FMPR$ROC, col=1, main="ROC", cex=2, lwd=3,
   ylim=c(0, 1), xlim=c(0, 1))
plot(res$GFlasso$ROC, col=2, add=TRUE, lty=2, lwd=3)
plot(res$Lasso$ROC, col=3, add=TRUE, lty=3, lwd=3)
abline(1, -1, lty=4, col=4, lwd=3)
legend(x=0, y=0.2,
   legend=c("FMPR-w2", "GFlasso-w2", "Lasso", "Random"),
   col=1:4, lwd=4, lty=1:4)

plot(res$FMPR$PRC, col=1, main="PRC", cex=2, lwd=3,
   ylim=c(0, 1), xlim=c(0, 1))
plot(res$GFlasso$PRC, col=2, add=TRUE, lty=2, lwd=3)
plot(res$Lasso$PRC, col=3, add=TRUE, lty=3, lwd=3)
abline(h=mean(B != 0), lty=4, col=4, lwd=3)



X2 <- scale(matrix(rnorm(N * p), N, p))
Y2 <- scale(X2 %*% B + rnorm(N * K))


sapply(list(FMPR=f, GFlasso=s, Lasso=l), function(B2) {
   pr <- X2 %*% B2
   R2(pr, Y2)
})


dev.off()

