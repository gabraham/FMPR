
library(FMPR)

options(error=dump.frames)

source("makematlab.R")

nreps <- 50

matlab.path <- "/usr/local/MATLAB/R2012a/bin/matlab"

load("../HAPGEN/sim_9/chr10.RData")
type <- "geno"
stat <- mean

cortype <- 2
corthresh <- 0
tol <- 1e-6
mu <- 1e-4
gamma <- 1e-3
verbose <- FALSE
maxiter <- 1e3

run <- function(dd)
{
   lapply(1:nrow(dd), function(i) {
      N <- dd$N[i]
      p <- dd$p[i]
      K <- dd$K[i]
      
      cat("N:", N, "p:", p, "K:", K, "\n")

      sapply(1:nreps, function(rep) {
         cat("rep", rep, "\n")
         B <- getB(p=p, K=K, w=0.5, type="same")
         d1 <- makedata(rep=rep, N=N, K=K, B=B, p=p, save=FALSE, sigma=2,
   	 type=type, maketest=FALSE)
         C <- gennetwork(d1$Ytrain, corthresh=corthresh, cortype=cortype)
         l <- max(abs(crossprod(d1$Xtrain, d1$Ytrain)))

         s.fmpr <- system.time({
	    f <- fmpr(X=d1$Xtrain, Y=d1$Ytrain, C=C, lambda=l * 0.9 / N,
	       gamma=gamma, verbose=verbose, maxiter=maxiter, eps=tol)
         })
   
         s.spg <- system.time({
   	 s <- spg(X=d1$Xtrain, Y=d1$Ytrain, C=C, lambda=l * 0.9,
   	    gamma=gamma, verbose=verbose, maxiter=maxiter,
   	    tol=tol, mu=mu, simplify=TRUE)
         })
   
         write.table(d1$Xtrain, file=sprintf("Xtrain_%d.txt", rep),
   	 row.names=FALSE, col.names=FALSE, quote=FALSE)
         write.table(d1$Ytrain, file=sprintf("Ytrain_%d.txt", rep),
   	 row.names=FALSE, col.names=FALSE, quote=FALSE)
   
	 m <- makematlab(spg.path="~/Software/SPG_Multi_Graph",
	    xfile=sprintf("Xtrain_%d.txt", rep),
	    yfile=sprintf("Ytrain_%d.txt", rep),
	    cortype=cortype, corthresh=corthresh,
	    maxiter=maxiter, tol=tol, mu=mu,
	    lambda=l * 0.9, gamma=gamma
         )
         cat(m, file="run.m")
         mt <- system(
	    sprintf("%s -nodisplay -nojvm -r run", matlab.path),
	    intern=TRUE
         )
   
         B.mat <- matrix(scan("B.txt", numeric()), nrow=p, byrow=TRUE)
         mse <- mean((B.mat - s)^2)
         sign.err <- mean(sign(B.mat) != sign(s))
         cat("mse:", mse, "sign error:", sign.err, "\n")
   
         s.spg.mat <- as.numeric(
	    strsplit(grep("^time:", mt, value=TRUE), ":")[[1]][[2]]
         )
   	 
         c(
	    FMPR=as.numeric(s.fmpr[3]), 
	    SPG=as.numeric(s.spg[3]),
	    SPG.mat=s.spg.mat
         )
      })
   })
}

# Time with increasing N
Ns <- c(100, 250, 500, 750, 1000, 1500, 2000, 2500)
names(Ns) <- Ns
p <- 500
K <- 10
dN <- data.frame(N=Ns, p=p, K=K)
res1 <- run(dN)
res1m <- sapply(res1, apply, 1, stat)

# Time with increasing p
N <- 100
ps <- c(100, 200, 300, 400, 500, 1000, 1500, 2000)
names(ps) <- ps
K <- 10
dp <- data.frame(N=N, p=ps, K=K)
res2 <- run(dp)
res2m <- sapply(res2, apply, 1, stat)

# Time with increasing K
N <- 1000
p <- 500
Ks <- c(2, 5, 10, 25, 50, 75, 100)
names(Ks) <- Ks
dK <- data.frame(N=N, p=p, K=Ks)
res3 <- run(dK)
res3m <- sapply(res3, apply, 1, stat)

save(Ns, ps, Ks, res1, res1m, res2, res2m, res3m, res3, stat,
      file="timing.RData")

