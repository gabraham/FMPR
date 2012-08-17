
# Wrappers for running the different models on the simulation data

# Make the simulated data.
#
# All inputs X and outputs Y are standardised to zero-mean and unit-variance
makedata <- function(rep, dir=".", N=100, p=50, K=5, B, sigma=0.01,
      save=TRUE)
{
   Xtrain <- scale(matrix(rnorm(N * p), N, p))
   Xtest <- scale(matrix(rnorm(N * p), N, p))

   noiseTrain <- rnorm(N * K, 0, sigma)
   noiseTest <- rnorm(N * K, 0, sigma)
   Ytrain <- scale(Xtrain %*% B + noiseTrain)
   Ytest <- scale(Xtest %*% B + noiseTest)

   if(save) {
      write.table(B, file=sprintf("%s/B_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)

      write.table(Xtrain, file=sprintf("%s/Xtrain_%s.txt", dir, rep),
            col.names=FALSE, row.names=FALSE, quote=FALSE)
      write.table(Ytrain, file=sprintf("%s/Ytrain_%s.txt", dir, rep),
            col.names=FALSE, row.names=FALSE, quote=FALSE)

      write.table(Xtest, file=sprintf("%s/Xtest_%s.txt", dir, rep),
            col.names=FALSE, row.names=FALSE, quote=FALSE)
      write.table(Ytest, file=sprintf("%s/Ytest_%s.txt", dir, rep),
            col.names=FALSE, row.names=FALSE, quote=FALSE)
   } else {
      list(
	 Xtrain=Xtrain,
	 Xtest=Xtest,
	 Ytrain=Ytrain,
	 Ytest=Ytest,
	 noiseTrain=noiseTrain,
	 noiseTest=noiseTest
      )
   }
}

run.spg <- function(rep, dir=".", nfolds=10, grid=25,
   lambdar=2^seq(-10, 0, length=grid),
   gamma=10^seq(-3, 6, length=grid),
   corthresh=0.5, cortype=1)
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
   l <- max(maxlambda1(Xtrain, Ytrain))
   lambda <- l * lambdar
   opt <- optim.spg(X=Xtrain, Y=Ytrain,
	 nfolds=nfolds,
	 corthresh=corthresh, cortype=cortype,
	 lambda=lambda, gamma=gamma,
	 maxiter=1e5)
   cat("optim.spg end\n")
   g <- spg(X=Xtrain, Y=Ytrain, C=C,
	 lambda=opt$opt["lambda"], gamma=opt$opt["gamma"], simplify=TRUE)

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
      lambda=10^seq(-3, 6, length=grid))
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
   
   opt <- optim.ridge(X=Xtrain, Y=Ytrain, nfolds=nfolds, lambda=lambda)
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
   graph.thresh=0.5, graph.fun=graph.sqr,
   lambdar=2^seq(-10, 0, length=grid),
   gamma=c(0, 10^seq(-3, 6, length=grid)),
   lambda2=c(0, 10^seq(-3, 6, length=grid)), type="l2", huber_mu=NA)
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
   R <- cor(Ytrain)
   G <- graph.fun(R, graph.thresh)
   
   if(all(G == 0))
   {
      thresh <- median(abs(R[upper.tri(R)]))
      cat("G is all zero in run.fmpr, using threshold of", thresh, "\n")
      G <- graph.fun(R, thresh)
   }

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   l <- max(maxlambda1(Xtrain, Ytrain))
   lambda <- l * lambdar

   cat("optim.fmpr start\n")
   opt <- optim.fmpr(X=Xtrain, Y=Ytrain,
	 nfolds=nfolds,
	 graph.fun=graph.fun, graph.thresh=graph.thresh,
	 lambda=lambda, gamma=gamma, lambda2=lambda2, maxiter=1e5, type=type,
	 huber_mu=huber_mu)
   cat("optim.fmpr end\n")
   g <- fmpr(X=Xtrain, Y=Ytrain,
	 lambda=opt$opt["lambda"],
	 gamma=opt$opt["gamma"],
	 lambda2=opt$opt["lambda2"],
	 maxiter=1e5, G=G, simplify=TRUE, type=type, huber_mu=huber_mu)

   P <- Xtest %*% g
   res <- R2(as.numeric(P), as.numeric(Ytest))

   cat("rep", rep, "R2 fmpr", res, "\n")

   setwd(oldwd)

   list(
      R2=res,
      beta=g
   )
}

run.fmpr2 <- function(rep, dir=".", nfolds=10, grid=25,
   corthresh=0.5, cortype=2,
   lambdar=2^seq(-10, 0, length=grid),
   gamma=c(0, 10^seq(-3, 6, length=grid)),
   lambda2=c(0, 10^seq(-3, 6, length=grid)), type="l2", huber_mu=NA)
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
   #R <- cor(Ytrain)
   #G <- graph.fun(R, graph.thresh)

   C <- gennetwork(Ytrain, cortype=cortype, corthresh=corthresh)
   
   #if(all(G == 0))
   #{
   #   thresh <- median(abs(R[upper.tri(R)]))
   #   cat("G is all zero in run.fmpr, using threshold of", thresh, "\n")
   #   G <- graph.fun(R, thresh)
   #}

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   l <- max(maxlambda1(Xtrain, Ytrain))
   lambda <- l * lambdar

   cat("optim.fmpr start\n")
   opt <- optim.fmpr2(X=Xtrain, Y=Ytrain,
	 nfolds=nfolds,
	 cortype=cortype, corthresh=corthresh,
	 lambda=lambda, gamma=gamma, lambda2=lambda2, maxiter=1e5, type=type,
	 huber_mu=huber_mu)
   cat("optim.fmpr end\n")
   g <- fmpr2(X=Xtrain, Y=Ytrain,
	 lambda=opt$opt["lambda"],
	 gamma=opt$opt["gamma"],
	 lambda2=opt$opt["lambda2"],
	 maxiter=1e5, C=C, simplify=TRUE, type=type, huber_mu=huber_mu)

   P <- Xtest %*% g
   res <- R2(as.numeric(P), as.numeric(Ytest))

   cat("rep", rep, "R2 fmpr2", res, "\n")

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
#       same: same weights, same sparsity
#       sparsity: different weights, same sparsity
#       random: different weights and different sparsity
#       mixed: same absolute weights with different sign, same sparsity
#       cluster; some tasks are related with same weight, some aren't
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
recovery <- function(obj, rep, dir, cleanROCR)
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
	 N=setup$N, p=setup$p, K=setup$K, B=setup$B, sigma=setup$sigma)
   cat("Simulation done\n")

   cat("Running methods\n")
   r.fmpr2.w1 <- lapply(1:nreps, run.fmpr2, dir=dir, grid=grid, nfolds=nfolds,
   	 cortype=1, corthresh=0)

   r.fmpr2.w2 <- lapply(1:nreps, run.fmpr, dir=dir, grid=grid, nfolds=nfolds,
   	 cortype=2, corthresh=0)

   r.fmpr.w1 <- lapply(1:nreps, run.fmpr, dir=dir, grid=grid, nfolds=nfolds,
   	 graph.fun=graph.abs, graph.thresh=0)

   r.fmpr.w2 <- lapply(1:nreps, run.fmpr, dir=dir, grid=grid, nfolds=nfolds,
   	 graph.fun=graph.sqr, graph.thresh=0)

   r.fmpr.h1.w1 <- lapply(1:nreps, run.fmpr, dir=dir, grid=grid, nfolds=nfolds,
   	 graph.fun=graph.abs, graph.thresh=0, type="huber", huber_mu=1e-3)

   r.fmpr.h1.w2 <- lapply(1:nreps, run.fmpr, dir=dir, grid=grid, nfolds=nfolds,
   	 graph.fun=graph.sqr, graph.thresh=0, type="huber", huber_mu=1e-3)

   r.fmpr.h2.w1 <- lapply(1:nreps, run.fmpr, dir=dir, grid=grid, nfolds=nfolds,
   	 graph.fun=graph.abs, graph.thresh=0, type="huber", huber_mu=1)

   r.fmpr.h2.w2 <- lapply(1:nreps, run.fmpr, dir=dir, grid=grid, nfolds=nfolds,
   	 graph.fun=graph.sqr, graph.thresh=0, type="huber", huber_mu=1)

   r.spg.w1 <- lapply(1:nreps, run.spg, dir=dir, grid=grid, nfolds=nfolds,
   	 corthresh=0, cortype=1)

   r.spg.w2 <- lapply(1:nreps, run.spg, dir=dir, grid=grid, nfolds=nfolds,
   	 corthresh=0, cortype=2)

   r.lasso <- lapply(1:nreps, run.fmpr, dir=dir,
	 grid=grid, nfolds=nfolds, gamma=0, lambda2=0)

   r.ridge <- lapply(1:nreps, run.ridge, dir=dir, grid=grid, nfolds=nfolds)

   r.elnet <- lapply(1:nreps, run.fmpr, dir=dir,
   	 grid=grid, nfolds=nfolds, gamma=0)
   cat("Done\n")
     
   R2.fmpr.w1 <- sapply(r.fmpr.w1, function(x) x$R2)
   R2.fmpr.w2 <- sapply(r.fmpr.w2, function(x) x$R2)
   R2.fmpr2.w1 <- sapply(r.fmpr2.w1, function(x) x$R2)
   R2.fmpr2.w2 <- sapply(r.fmpr2.w2, function(x) x$R2)
   R2.fmpr.h1.w1 <- sapply(r.fmpr.h1.w1, function(x) x$R2)
   R2.fmpr.h1.w2 <- sapply(r.fmpr.h1.w2, function(x) x$R2)
   R2.fmpr.h2.w1 <- sapply(r.fmpr.h2.w1, function(x) x$R2)
   R2.fmpr.h2.w2 <- sapply(r.fmpr.h2.w2, function(x) x$R2)
   R2.spg.w1 <- sapply(r.spg.w1, function(x) x$R2)
   R2.spg.w2 <- sapply(r.spg.w2, function(x) x$R2)
   R2.lasso <- sapply(r.lasso, function(x) x$R2)
   R2.ridge <- sapply(r.ridge, function(x) x$R2)
   R2.elnet <- sapply(r.elnet, function(x) x$R2)
  
   R2.all <- cbind(
      FMPRw1=R2.fmpr.w1,
      FMPRw2=R2.fmpr.w2,
      FMPR2w1=R2.fmpr2.w1,
      FMPR2w2=R2.fmpr2.w2,
      FMPRh1w1=R2.fmpr.h1.w1,
      FMPRh1w2=R2.fmpr.h1.w2,
      FMPRh2w1=R2.fmpr.h2.w1,
      FMPRh2w2=R2.fmpr.h2.w2,
      SPGw1=R2.spg.w1,
      SPGw2=R2.spg.w2,
      Lasso=R2.lasso,
      Ridge=R2.ridge,
      ElasticNet=R2.elnet
   )

   rec.fmpr.w1 <- recovery(r.fmpr.w1, rep, dir, cleanROCR)
   rec.fmpr.w2 <- recovery(r.fmpr.w2, rep, dir, cleanROCR)
   rec.fmpr2.w1 <- recovery(r.fmpr2.w1, rep, dir, cleanROCR)
   rec.fmpr2.w2 <- recovery(r.fmpr2.w2, rep, dir, cleanROCR)
   rec.fmpr.h1.w1 <- recovery(r.fmpr.h1.w1, rep, dir, cleanROCR)
   rec.fmpr.h1.w2 <- recovery(r.fmpr.h2.w2, rep, dir, cleanROCR)
   rec.fmpr.h2.w1 <- recovery(r.fmpr.h2.w1, rep, dir, cleanROCR)
   rec.fmpr.h2.w2 <- recovery(r.fmpr.h2.w2, rep, dir, cleanROCR)
   rec.spg.w1 <- recovery(r.spg.w1, rep, dir, cleanROCR)
   rec.spg.w2 <- recovery(r.spg.w2, rep, dir, cleanROCR)
   rec.lasso <- recovery(r.lasso, rep, dir, cleanROCR)
   rec.ridge <- recovery(r.ridge, rep, dir, cleanROCR)
   rec.elnet <- recovery(r.elnet, rep, dir, cleanROCR)

   ex <- list(
      dir=dir,
      weights=list(
	 fmpr.w1=r.fmpr.w1,
	 fmpr.w2=r.fmpr.w2,
	 fmpr2.w1=r.fmpr2.w1,
	 fmpr2.w2=r.fmpr2.w2,
	 fmpr.h1.w1=r.fmpr.h1.w1,
	 fmpr.h1.w2=r.fmpr.h1.w2,
	 fmpr.h2.w1=r.fmpr.h2.w1,
	 fmpr.h2.w2=r.fmpr.h2.w2,
	 lasso=r.lasso,
	 ridge=r.ridge,
	 elnet=r.elnet,
	 spg.w1=r.spg.w1,
	 spg.w2=r.spg.w2
      ),
      recovery=list(
	 fmpr.w1=rec.fmpr.w1,
	 fmpr.w2=rec.fmpr.w2,
	 fmpr2.w1=rec.fmpr2.w1,
	 fmpr2.w2=rec.fmpr2.w2,
	 fmpr.h1.w1=rec.fmpr.h1.w1,
	 fmpr.h1.w2=rec.fmpr.h1.w2,
	 fmpr.h2.w1=rec.fmpr.h2.w1,
	 fmpr.h2.w2=rec.fmpr.h2.w2,
	 lasso=rec.lasso,
	 ridge=rec.ridge,
	 elnet=rec.elnet,
	 spg.w1=rec.spg.w1,
	 spg.w2=rec.spg.w2
      ),
      R2=R2.all
   )
   class(ex) <- "exper"
   ex
}

