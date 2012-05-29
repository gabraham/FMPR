
library(FMPR)

nreps <- 50

# Time with increasing N
Ns <- c(100, 250, 500, 750, 1000, 2000, 3000, 4000, 5000)
names(Ns) <- Ns
p <- 400
K <- 10

res1 <- lapply(Ns, function(N) {
   cat("N:", N, "\n")
   sapply(1:nreps, function(rep) {
      B <- getB(p=p, K=K, w=0.5, type="same")
      d1 <- makedata(N=N, K=K, B=B, p=p, save=FALSE, sigma=2, rep=1)
      G <- graph.sqr(cor(d1$Ytrain))
      C <- gennetwork(d1$Ytrain, corthresh=0, cortype=2)
      l <- max(maxlambda1(d1$Xtrain, d1$Ytrain))
      s.fmpr <- system.time({
	 fmpr(X=d1$Xtrain, Y=d1$Ytrain, G=G, lambda=l * 0.9, gamma=1)
      })
      s.spg <- system.time({
	 spg(X=d1$Xtrain, Y=d1$Ytrain, C=C, lambda=l * 0.9,
	    gamma=1)
      })
      cat("rep", rep, "FMPR:", s.fmpr, "SPG:", s.spg, "\n")
      c(FMPR=s.fmpr[3], SPG=s.spg[3])
   })
})

res1m <- sapply(res1, apply, 1, mean)

# Time with increasing p
N <- 100
ps <- c(100, 200, 300, 400, 500, 1000, 1500, 2000)
names(ps) <- ps
K <- 10

res2 <- lapply(ps, function(p) {
   cat("p:", p, "\n")
   sapply(1:nreps, function(rep) {
      B <- getB(p=p, K=K, w=0.5, type="same")
      d1 <- makedata(N=N, K=K, B=B, p=p, save=FALSE, sigma=2, rep=1)
      G <- graph.sqr(cor(d1$Ytrain))
      C <- gennetwork(d1$Ytrain, corthresh=0, cortype=2)
      l <- max(maxlambda1(d1$Xtrain, d1$Ytrain))
      s.fmpr <- system.time({
	 fmpr(X=d1$Xtrain, Y=d1$Ytrain, G=G, lambda=l * 0.9, gamma=1)
      })
      s.spg <- system.time({
	 spg(X=d1$Xtrain, Y=d1$Ytrain, C=C, lambda=l * 0.9, gamma=1)
      })
      cat("rep:", rep, "FMPR:", s.fmpr, "SPG:", s.spg, "\n")
      c(FMPR=s.fmpr[3], SPG=s.spg[3])
   })
})

res2m <- sapply(res2, apply, 1, mean)

# Time with increasing K
N <- 100
p <- 100
Ks <- c(2, 5, 10, 25, 50, 75, 100, 150, 200)

res3 <- lapply(Ks, function(K) {
   cat("K:", K, "\n")
   sapply(1:nreps, function(rep) {
      B <- getB(p=p, K=K, w=0.5, type="same")
      d1 <- makedata(N=N, K=K, B=B, p=p, save=FALSE, sigma=2, rep=1)
      G <- graph.sqr(cor(d1$Ytrain))
      C <- gennetwork(d1$Ytrain, corthresh=0, cortype=2)
      l <- max(maxlambda1(d1$Xtrain, d1$Ytrain))
      s.fmpr <- system.time({
	 fmpr(X=d1$Xtrain, Y=d1$Ytrain, G=G, lambda=l * 0.9, gamma=1)
      })
      s.spg <- system.time({
	 spg(X=d1$Xtrain, Y=d1$Ytrain, C=C, lambda=l * 0.9, gamma=1)
      })
      cat("rep:", rep, "FMPR:", s.fmpr, "SPG:", s.spg, "\n")
      c(FMPR=s.fmpr[3], SPG=s.spg[3])
   })
})

res3m <- sapply(res3, apply, 1, mean)

save(Ns, ps, Ks, res1, res1m, res2, res2m, res3m, res3,
      file="timing.RData")

pdf("FMPR_scaling_samples.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(Ns, t(res1m), xlab="Samples N", ylab="Time (sec)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
legend(x=0, y=max(res1m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2,
      lty=1:2)
dev.off()

pdf("FMPR_scaling_samples_scaled.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(Ns, scale(t(res1m)), xlab="Samples N", ylab="Time (scaled)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
#legend(x=0, y=max(res1m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2)
dev.off()

pdf("FMPR_scaling_variables.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(ps, t(res2m), xlab="Variables p", ylab="Time (sec)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
legend(x=50, y=max(res2m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2,
      lty=1:2)
dev.off()

pdf("FMPR_scaling_variables_scaled.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(ps, scale(t(res2m)), xlab="Variables p", ylab="Time (scaled)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
#legend(x=50, y=max(res2m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2)
dev.off()

pdf("FMPR_scaling_tasks.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(Ks, t(res3m), xlab="Tasks K", ylab="Time (sec)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
legend(x=0, y=max(res3m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2,
      lty=1:2)
dev.off()

pdf("FMPR_scaling_tasks_scaled.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(Ks, scale(t(res3m)), xlab="Tasks K", ylab="Time (scaled)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
#legend(x=0, y=max(res3m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2)
dev.off()

