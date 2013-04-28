
# Wrappers for running the different models on the simulation data

# Make the simulated data.
#
# All inputs X and outputs Y are standardised to zero-mean and unit-variance
makedata <- function(rep, dir=".", N=100, p=50, K=5, B, sigma=0.01,
      save=TRUE, type=c("real", "artificial"), maketest=TRUE)
{
   type <- match.arg(type)

   Xtest <- Ytest <- Xtest2 <- Ytest2 <- noiseTest <- NULL

   if(type == "artificial") {
      cat("makedata: using artificial data\n")
      Xtrain <- matrix(sample(0:2, N * p, replace=TRUE), N, p)
      if(maketest) {
	 Xtest <- matrix(sample(0:2, N * p, replace=TRUE), N, p)
      }
   } else if(!exists("geno", mode="numeric") || is.null(geno)) {
      stop("type = real but geno is missing or NULL")
   } else if(N > nrow(geno)) {
      stop("N too big: geno only has ", nrows(geno), " rows")
   } else {
      cat("makedata: using real data\n")
      # randomly split samples into training and testing, and from each set
      # select N samples
      wtrain <- sample(nrow(geno), N)
      # select a contiguous block of SNPs of size p
      snp.start <- sample(ncol(geno) - p, 1)
      snps <- snp.start + 1:p - 1

      if(maketest) {
	 sample.spl <- sample(c(TRUE, FALSE), nrow(geno), replace=TRUE)
	 wtrain <- sample(which(sample.spl), N)
	 wtest <- sample(which(!sample.spl), N)
	 Xtest <- geno[wtest, snps]
      }
      Xtrain <- geno[wtrain, snps]
   }

   noiseTrain <- rnorm(nrow(Xtrain) * K, 0, sigma)
   Ytrain <- Xtrain %*% B + noiseTrain
   Xtrain2 <- scalefix(Xtrain)
   Ytrain2 <- scalefix(Ytrain)

   if(maketest) {
      noiseTest <- rnorm(nrow(Xtrain) * K, 0, sigma)
      Ytest <- Xtest %*% B + noiseTest
      Xtest2 <- scalefix(Xtest)
      Ytest2 <- scalefix(Ytest)
   }

   if(save) {
      write.table(B, file=sprintf("%s/B_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)

      write.table(Xtrain2, file=sprintf("%s/Xtrain_%s.txt", dir, rep),
            col.names=FALSE, row.names=FALSE, quote=FALSE)
      write.table(Xtrain, file=sprintf("%s/Xtrain_orig_%s.txt", dir, rep),
            col.names=FALSE, row.names=FALSE, quote=FALSE)

      write.table(Ytrain2, file=sprintf("%s/Ytrain_%s.txt", dir, rep),
            col.names=FALSE, row.names=FALSE, quote=FALSE)
      write.table(Ytrain, file=sprintf("%s/Ytrain_orig_%s.txt", dir, rep),
            col.names=FALSE, row.names=FALSE, quote=FALSE)

      if(maketest) {
	 write.table(Xtest2, file=sprintf("%s/Xtest_%s.txt", dir, rep),
      	       col.names=FALSE, row.names=FALSE, quote=FALSE)
      	 write.table(Xtest, file=sprintf("%s/Xtest_orig_%s.txt", dir, rep),
      	       col.names=FALSE, row.names=FALSE, quote=FALSE)

      	 write.table(Ytest2, file=sprintf("%s/Ytest_%s.txt", dir, rep),
      	       col.names=FALSE, row.names=FALSE, quote=FALSE)
      	 write.table(Ytest, file=sprintf("%s/Ytest_orig_%s.txt", dir, rep),
      	       col.names=FALSE, row.names=FALSE, quote=FALSE)
      } 
   } else {
      list(
	 Xtrain=Xtrain2,
	 XtrainOrig=Xtrain,
	 Xtest=Xtest2,
	 XtestOrig=Xtest,
	 Ytrain=Ytrain2,
	 YtrainOrig=Ytrain,
	 Ytest=Ytest2,
	 YtestOrig=Ytest,
	 noiseTrain=noiseTrain,
	 noiseTest=noiseTest
      )
   }
}

run.spg <- function(rep, dir=".", nfolds=10, grid=25,
   lambdar=2^seq(-10, 0, length=grid),
   gamma=10^seq(-3, 6, length=grid),
   corthresh=0, cortype=1, type="l1", divbyN=FALSE)
{
   oldwd <- getwd()
   setwd(dir)

   Xtrain <- as.matrix(read.table(sprintf("Xtrain_%s.txt", rep),
	    header=FALSE))
   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep), header=FALSE))
   Ytrain <- as.matrix(read.table(sprintf("Ytrain_%s.txt", rep), header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep), header=FALSE))
   
   C <- gennetwork(Ytrain, corthresh, cortype)
   if(length(C) == 0)
   {
      R <- cor(Ytrain)
      m <- median(abs(R[upper.tri(R)]))
      cat("run.spg: increasing Rthresh to", m, "\n")
      C <- gennetwork(Ytrain, corthresh, cortype)
   }

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   cat("optim.spg start rep", rep, "\n")
   l <- max(abs(crossprod(Xtrain, Ytrain)))
   if(divbyN) {
      l <- l / N
   }
   lambda <- l * lambdar
   opt <- optim.spg(X=Xtrain, Y=Ytrain,
	 nfolds=nfolds, type=type,
	 corthresh=corthresh, cortype=cortype,
	 lambda=lambda, gamma=gamma, divbyN=divbyN,
	 maxiter=1e5)
   cat("optim.spg end\n")
   g <- spg(X=Xtrain, Y=Ytrain, C=C, type=type,
	 lambda=opt$opt["lambda"], gamma=opt$opt["gamma"],
	 simplify=TRUE, divbyN=divbyN)

   P <- Xtest %*% g
   res <- R2(as.numeric(P), as.numeric(Ytest))

   cat("finished rep", rep, "R2 SPG", res, "\n\n")

   setwd(oldwd)

   list(
      R2=res,
      beta=g
   )
}

run.ridge <- function(rep, dir=".", nfolds=10, grid=25,
      lambda2=10^seq(-3, 6, length=grid))
{
   oldwd <- getwd()
   setwd(dir)

   Xtrain <- as.matrix(read.table(sprintf("Xtrain_%s.txt", rep),
	    header=FALSE))
   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep), header=FALSE))
   Ytrain <- as.matrix(read.table(sprintf("Ytrain_%s.txt", rep),
	    header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep), header=FALSE))

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   ytest <- as.numeric(Ytest)
   
   opt <- optim.ridge(X=Xtrain, Y=Ytrain, nfolds=nfolds, lambda=lambda2)
   g <- ridge(Xtrain, Ytrain, lambda=opt$opt["lambda"])[[1]]
   P <- Xtest %*% g
   res <- R2(as.numeric(P), ytest)
   cat("rep", rep, "R2 ridge:", res, "\n")

   setwd(oldwd)

   list(
      R2=res,
      beta=matrix(g, p, K)
   )
}

run.fmpr <- function(rep, dir=".", nfolds=10, grid=25,
   corthresh=0, cortype=2,
   lambdar=2^seq(-10, 0, length=grid),
   gamma=c(0, 10^seq(-3, 6, length=grid)),
   lambda2=c(0, 10^seq(-3, 6, length=grid)))
{
   oldwd <- getwd()
   setwd(dir)

   cat("run.fmpr rep", rep, "\n")

   Xtrain <- as.matrix(read.table(sprintf("Xtrain_%s.txt", rep),
      header=FALSE))
   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep),
      header=FALSE))
   Ytrain <- as.matrix(read.table(sprintf("Ytrain_%s.txt", rep),
      header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep),
      header=FALSE))

   C <- gennetwork(Ytrain, cortype=cortype, corthresh=corthresh)
   
   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   l <- max(abs(crossprod(Xtrain, Ytrain))) / N
   lambda <- l * lambdar

   cat("optim.fmpr start\n")
   opt <- optim.fmpr(X=Xtrain, Y=Ytrain,
	 nfolds=nfolds,
	 cortype=cortype, corthresh=corthresh,
	 lambda=lambda, gamma=gamma, lambda2=lambda2, maxiter=1e6)

   cat("optim.fmpr end\n")
   g <- fmpr(X=Xtrain, Y=Ytrain,
	 lambda=opt$opt["lambda"],
	 gamma=opt$opt["gamma"],
	 lambda2=opt$opt["lambda2"],
	 maxiter=1e5, C=C, simplify=TRUE)

   P <- Xtest %*% g
   res <- R2(as.numeric(P), as.numeric(Ytest))

   cat("rep", rep, "R2 fmpr", res, "\n")

   setwd(oldwd)

   list(
      R2=res,
      beta=g
   )
}

# p: no. variables per task
# K: no. tasks
# w: weight per variable, if using type="same"
# sparsity: [0,1], degree of sparsity per task
#
# type: 
#       same: same weights, same sparsity across all tasks
#       sparsity: different weights, same sparsity across all tasks
#       random: different weights and different sparsity across the tasks
#       mixed: same absolute weights with different sign, same sparsity
#       cluster: some tasks are related with same weight, some aren't
#       clustersparse: subsets of SNPs and subsets of tasks are related
getB <- function(p, K, w=0.1, sparsity=0.8, type=NULL, ...)
{
   if(type == "same") {
      b <- 0
      B <- 0
      while(all(B == 0)) {
	 b <- w * sample(0:1, p, TRUE, prob=c(sparsity, 1 - sparsity))
	 B <- sapply(1:K, function(k) b)
      }
   } else if(type == "sparsity") {
      B <- 0
      while(all(B == 0)) {
	 B <- matrix(rnorm(p * K, ...), p, K)
	 s <- sample(0:1, p, TRUE, prob=c(sparsity, 1 - sparsity))
	 B[s == 0, ] <- 0
      }
   } else if(type == "random") {
      B <- matrix(rnorm(p * K, ...), p, K)
      B <- B * sample(0:1, p * K, TRUE, prob=c(sparsity, 1 - sparsity))
   } else if(type == "mixed") { 
      b <- 0
      B <- 0
      while(all(B == 0)) {
	 b <- w * sample(0:1, p, TRUE, prob=c(sparsity, 1 - sparsity))
	 B <- sapply(1:K, function(k) sample(c(-1, 1), 1) * b)
      }
   } else if(type == "cluster") {
      b <- 0
      B <- 0
      while(all(B == 0)) {
	 b <- w * sample(0:1, p, TRUE, prob=c(sparsity, 1 - sparsity))
	 B <- sapply(1:K, function(k) sample(c(-1, 0, 1), 1) * b)
      }
   } else if(type == "clustersparse") {
      B <- 0
      while(all(B == 0)) {
	 B <- matrix(rnorm(p * K, ...), p, K)
	 s1 <- sample(0:1, p, TRUE, prob=c(sparsity, 1 - sparsity))
	 s2 <- sample(0:1, K, TRUE, prob=c(sparsity, 1 - sparsity))
	 B[s1 == 0, ] <- 0
	 B[, s2 == 0] <- 0
      }

   }

   B
}

# Measure recovery of non-zeros
recovery <- function(obj, dir, cleanROCR=TRUE)
{
   beta <- sapply(seq(along=obj), function(rep) {
      m <- as.matrix(read.table(sprintf("%s/B_%s.txt", dir, rep)))
      as.numeric(m)
   })

   betahat <- sapply(obj, function(r) as.numeric(abs(r$beta)))

   pred <- prediction(
         predictions=betahat,
         labels=beta != 0
   )
   roc <- performance(pred, "sens", "spec")
   prc <- performance(pred, "prec", "rec")
   auc <- performance(pred, "auc")

   if(cleanROCR) {
      roc <- clean.rocr(roc)
      prc <- clean.rocr(prc)
      auc <- clean.rocr(auc)
   }
   list(roc=roc, prc=prc, auc=auc)
}

# Evaluate methods over each setup
run <- function(setup, grid=3, nfolds=3, nreps=3, cleanROCR=TRUE)
{
   dir <- setup$dir

   if(!file.exists(dir)) {
      cat("Simulating data in", dir, "\n")
      dir.create(dir)
   }
   m <- sapply(1:nreps, makedata, dir=dir,
	 N=setup$N, p=setup$p, K=setup$K, B=setup$B, sigma=setup$sigma,
	 type=setup$type)
   cat("Simulation done\n")

   par.all <- list(dir=dir, grid=grid, nfolds=nfolds)
  
   # proportional to the largest lambda that makes all betas zero
   lambdar <- 2^seq(-10, 0, length=grid)
   gamma <- c(0, 10^seq(-3, 6, length=grid))
   lambda2 <- 10^seq(-3, 6, length=grid)

   param <- list(
      "FMPR-w1"=list(func=run.fmpr, cortype=1,
	 lambdar=lambdar, gamma=gamma, lambda2=0),
      "FMPR-w2"=list(func=run.fmpr, cortype=2,
	 lambdar=lambdar, gamma=gamma, lambda2=0),
      "GFlasso-l1-w1"=list(func=run.spg, cortype=1, lambdar=lambdar,
         gamma=gamma, type="l1"),
      "GFlasso-l1-w2"=list(func=run.spg, cortype=2, lambdar=lambdar,
         gamma=gamma, type="l1"),
      Lasso=list(func=run.fmpr, lambdar=lambdar, gamma=0, lambda2=0),
      Ridge=list(func=run.ridge, lambda2=lambda2),
      Elnet=list(func=run.fmpr,
	 lambdar=lambdar, lambda2=lambda2, gamma=0)
   )

   res <- lapply(seq(along=param), function(i) {
      lapply(1:nreps, function(rep) {
         do.call(param[[i]][[1]],
	    c(rep=rep, par.all, param[[i]][-1]))
      })
   })

   names(res) <- names(param)

   res.R2 <- sapply(res, sapply, function(x) x$R2)
   res.recovery <- lapply(res, recovery, dir=dir, cleanROCR=cleanROCR)
   ex <- list(
      dir=dir,
      weights=lapply(res, function(x) x$beta),
      recovery=res.recovery,
      R2=res.R2
   )
   class(ex) <- "exper"
   ex
}

