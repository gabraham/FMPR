# Test groupridge

library(ROCR)

options(error=dump.frames)

s <- sample(1e6L, 1)
set.seed(s)

source("methods.R")
source("eval.R")

makedata <- function(rep, dir=".", N=100, p=50, K=5, B)
{
   cat("rep", rep, "\n")

   Xtrain <- standardise(matrix(sample(0:2, N * p, replace=TRUE), N, p))
   Xtest <- standardise(matrix(sample(0:2, N * p, replace=TRUE), N, p))
   grp <- rep(1, K)
   XtrainB <- blockX(Xtrain, p, K)
   XtestB <- blockX(Xtest, p, K)

   XtrainBs <- standardise(XtrainB)
   XtestBs <- standardise(XtestB)
   
   #b <- 1 * sample(0:1, p, TRUE, prob=c(0.8, 0.2))
   #B <- sapply(1:K, function(k) b)
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
      R2=R2.spg,
      R2lasso=R2.spg.l,
      beta=B.spg,
      betalasso=B.spg.l
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
   R2.gl <- R2(as.numeric(P.gl), as.numeric(ytest))
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

# p: no. variables per task
# K: no. tasks
# w: weight per variable, if using type="same"
# sparsity: [0,1], degree of sparsity per task
getB <- function(p, K, w, sparsity=0.8, type=NULL)
{
   if(type == "same") {
      b <- 0
      B <- 0
      while(all(B == 0)) {
	 b <- w * sample(0:1, p, TRUE, prob=c(sparsity, 1 - sparsity))
	 B <- sapply(1:K, function(k) b)
      }
   }

   B
}

# Evaluate methods over each setup
run <- function(setup, grid=3, nfolds=3, nreps=3)
{
   dir <- setup$dir

   if(!file.exists(dir)) {
      cat("Simulating data in", dir, "\n")
      dir.create(dir)
      m <- sapply(1:nreps, makedata, dir=dir, N=setup$N, p=setup$p, K=setup$K,
	    B=setup$B)
      cat("Simulation done\n")
   }
   
   cat("Running inference\n")
   r.gr <- lapply(1:nreps, run.groupridge, dir=dir, r=grid, nfolds=nfolds)
   r.lasso <- lapply(1:nreps, run.lasso, dir=dir, r=grid, nfolds=nfolds)
   r.spg <- lapply(1:nreps, run.spg, dir=dir, r=grid, nfolds=nfolds)
   r.ridge <- lapply(1:nreps, run.ridge, dir=dir, r=grid, nfolds=nfolds)
   cat("Inference done\n")
   
   # Measure recovery of non-zeros
   recovery <- function(obj, dir)
   {
      lapply(seq(along=obj), function(rep) {
         beta <- as.matrix(read.table(sprintf("%s/B_%s.txt", dir, rep)))
         pred <- prediction(predictions=abs(obj[[rep]]$beta), labels=beta != 0)
         roc <- performance(pred, "sens", "spec")
         prc <- performance(pred, "prec", "rec")
         list(roc=roc, prc=prc)
      })
   }
   
   R2.gr <- sapply(r.gr, function(x) x$R2)
   R2.lasso <- sapply(r.lasso, function(x) x$R2)
   R2.ridge <- sapply(r.ridge, function(x) x$R2)
   R2.spg <- sapply(r.spg, function(x) x$R2[1])
   
   R2.all <- cbind(GR=R2.gr, lasso=R2.lasso, ridge=R2.ridge, SPG=R2.spg)
   
   rec.gr <- recovery(r.gr, dir)
   rec.lasso <- recovery(r.lasso, dir)
   rec.ridge <- recovery(r.ridge, dir)
   rec.spg <- recovery(r.spg, dir)

   list(
      recovery=list(gr=rec.gr, lasso=rec.lasso,
	    ridge=rec.ridge, spg=rec.spg),
      R2=R2.all
   )
}

setup <- list(
   
   # Change sample size, all else fixed
   list(dir=c("Expr1"), N=50, p=50, K=10,
	 B=getB(p=10, K=3, w=0.5, type="same")),
   list(dir=c("Expr2"), N=100, p=50, K=10,
	 B=getB(p=10, K=3, w=0.5, type="same"))
   list(dir=c("Expr2"), N=150, p=50, K=10,
	 B=getB(p=10, K=3, w=0.5, type="same"))
)

res <- lapply(setup, run, nreps=20)
save(setup, res, file="results.RData")

