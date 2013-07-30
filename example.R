### R code from vignette source 'fmpr.rnw'

library(FMPR)
if(require(doMC)) {
   registerDoMC(cores=3)
}

N <- 100
p <- 100
K <- 10
w <- 0.1
B <- getB(p=p, K=K, w=w, type="same", sparsity=0.8)
print(table(sign(B)))

X <- scale(matrix(rnorm(N * p), N, p))
Y <- scale(X %*% B + rnorm(N * K, 0, 1))
C <- gennetwork(Y, cortype=2)

l <- max(maxlambda1(X, Y))
ngrid <- 20
lambda <- 2^seq(-0.01, -8, length=ngrid) * l
gamma <- c(0, 2^seq(-10, 10, length=ngrid))
nfolds <- 10
maxiter <- 1e5

system.time({
   opt.f <- optim.fmpr(X=X, Y=Y, cortype=2,
      lambda=lambda, gamma=gamma, nfolds=nfolds)
})

# L1 fusion penalty
system.time({
   opt.s1 <- optim.spg(X=X, Y=Y, cortype=2,
      lambda=lambda, gamma=gamma, nfolds=nfolds, type="l1",
      maxiter=maxiter)
})

f <- fmpr(X=X, Y=Y, C=C, lambda=opt.f$opt["lambda"],
   gamma=opt.f$opt["gamma"], simplify=TRUE)

s1 <- spg(X=X, Y=Y, C=C, lambda=opt.s1$opt["lambda"],
   gamma=opt.s1$opt["gamma"], simplify=TRUE, type="l1")

# Lasso, special case of FMPR
w <- which(opt.f$R2[, 1, 1] == max(opt.f$R2[, 1, 1]), arr.ind=TRUE)
l <- fmpr(X=X, Y=Y, lambda=lambda[w], gamma=0, simplify=TRUE)

measures <- list(ROC=c("sens", "spec"), PRC=c("prec", "rec"))
mods <- list(FMPR=f, GFlasso=s1, Lasso=l)
res <- lapply(mods, function(B2) {
   lapply(measures, function(m) {
      performance(
	 prediction(
	    labels=as.numeric(B != 0),
	    predictions=as.numeric(abs(B2))
	 ), m[1], m[2]
      )
   })
})


pdf("example.pdf", width=14)

par(mfrow=c(1, 2), mar=c(4, 4, 2, 2) + 0.1)

plot(NULL, ylim=c(0, 0.4), xlim=c(0.95, 1), main="ROC",
   xlab="Specificity", ylab="Sensitivity")
for(i in seq(along=res)) {
   plot(res[[i]]$ROC, col=i, lty=i, add=TRUE, lwd=3)
}

abline(1, -1, lty=4, col=5, lwd=3)

plot(NULL, ylim=c(0, 1), xlim=c(0, 1), main="PRC",
   xlab="Recall", ylab="Precision")
for(i in seq(along=res)) {
   plot(res[[i]]$PRC, col=i, lty=i, add=TRUE, lwd=3)
}

legend(x=0.7, y=0.9, legend=c(names(mods), "Random"),
   col=1:5, lwd=4, lty=1:5)
abline(h=mean(B != 0), lty=4, col=5, lwd=3)

dev.off()

X2 <- scale(matrix(rnorm(N * p), N, p))
Y2 <- scale(X2 %*% B + rnorm(N * K))

r2 <- sapply(mods, function(B2) {
   pr <- X2 %*% B2
   R2(pr, Y2)
})

sort(r2)

