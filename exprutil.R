
# Wrappers for running the different models on the simulation data

makedata <- function(rep, dir=".", N=100, p=50, K=5, B, sigma=0.01,
      save=TRUE)
{
   cat("rep", rep, "\n")

   Xtrain <- matrix(rnorm(N * p), N, p)
   Xtest <- matrix(rnorm(N * p), N, p)

   noiseTrain <- rnorm(N * K, 0, sigma)
   noiseTest <- rnorm(N * K, 0, sigma)
   Ytrain <- Xtrain %*% B + noiseTrain
   Ytest <- Xtest %*% B + noiseTest

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

   N <- nrow(Xtrain)
   p <- ncol(Xtrain)
   K <- ncol(Ytrain)

   XtrainB <- scale(blockX(Xtrain, p, K))
   XtestB <- scale(blockX(Xtest, p, K))

   ytrain <- as.numeric(center(Ytrain))
   ytest <- as.numeric(center(Ytest))

   l <- maxlambda1(XtrainB, ytrain)
   opt <- optim.lasso(X=XtrainB, Y=ytrain, nfolds=nfolds,
	 L1=seq(l, l * 1e-3, length=r), verbose=FALSE)
   g <- lasso3(X=XtrainB, y=ytrain, lambda1=opt$opt[1])
   P <- XtestB %*% g
   res <- R2(as.numeric(P), ytest)
   cat("rep", rep, "R2 lasso:", res, "\n")

   #if(all(g == 0))
    #  browser()

   setwd(oldwd)

   list(
      R2=res,
      beta=matrix(g, p, K)
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

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   XtrainB <- scale(blockX(Xtrain, p, K))
   XtestB <- scale(blockX(Xtest, p, K)) 

   ytrain <- as.numeric(center(Ytrain))
   ytest <- as.numeric(center(Ytest))
   
   r <- optim.ridge(X=XtrainB, Y=ytrain, nfolds=nfolds,
	 L2=10^seq(-3, 5, length=r))
   g <- ridge(XtrainB, ytrain, lambda2=r$opt[1])
   P <- XtestB %*% g
   res <- R2(as.numeric(P), ytest)
   cat("rep", rep, "R2 ridge:", res, "\n")

   setwd(oldwd)

   list(
      R2=res,
      beta=matrix(g, p, K)
   )
}

run.groupridge <- function(rep, dir=".", nfolds=10, r=25, Rthresh=0.5)
{
   oldwd <- getwd()
   setwd(dir)

   Xtrain <- as.matrix(read.table(sprintf("Xtrain_%s.txt", rep),
	    header=FALSE))
   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep), header=FALSE))
   Ytrain <- as.matrix(read.table(sprintf("Ytrain_%s.txt", rep), header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep), header=FALSE))
   R <- cor(Ytrain)
   diag(R) <- 0
   G <- sign(R) * (abs(R) > Rthresh)
   if(all(G == 0))
      warning("G is all zero in run.groupridge, will not consider groups")

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   Xtrain <- scale(Xtrain)
   Ytrain <- center(Ytrain)

   cat("optim.groupridge start\n")
   l <- max(maxlambda1(Xtrain, Ytrain))
   r <- optim.groupridge(X=Xtrain, Y=Ytrain, G=G,
	 nfolds=nfolds,
	 L1=seq(l, l * 1e-3, length=r),
	 #L2=seq(0, 10, length=r),
	 #L3=seq(0, 10, length=r),
	 L2=10^seq(-3, 5, length=r),
	 L3=10^seq(-3, 5, length=r),
	 maxiter=1e3)
   cat("optim.groupridge end\n")
   g <- groupridge4(X=Xtrain, Y=Ytrain,
	 lambda1=r$opt[1], lambda2=r$opt[2], lambda3=r$opt[3],
	 maxiter=1e4, G=G)
   g <- g[[1]][[1]][[1]]

   P <- scale(Xtest) %*% g
   res <- R2(as.numeric(P), as.numeric(center(Ytest)))

   cat("rep", rep, "R2 groupridge:", res, "\n")

   setwd(oldwd)

   list(
      R2=res,
      beta=g
   )
}

run.elnet.fmpr <- function(rep, dir=".", nfolds=10, r=25, Rthresh=0.5)
{
   oldwd <- getwd()
   setwd(dir)

   Xtrain <- as.matrix(read.table(sprintf("Xtrain_%s.txt", rep),
	    header=FALSE))
   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep), header=FALSE))
   Ytrain <- as.matrix(read.table(sprintf("Ytrain_%s.txt", rep), header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep), header=FALSE))
  
   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   Xtrain <- scale(Xtrain)
   Ytrain <- center(Ytrain)

   G0 <- matrix(0, K, K)
   
   cat("optim.elnet.fmpr start\n")

   # We can use groupridge to do elastic net because we just ignore the task
   # relatedness and set lambda3 to zero.
   l <- max(maxlambda1(Xtrain, Ytrain))
   r <- optim.groupridge(X=Xtrain, Y=Ytrain, G=G0,
	 nfolds=nfolds,
	 L1=seq(l, l * 1e-3, length=r),
	 #L2=seq(0, 10, length=r),
	 L2=10^seq(-3, 5, length=r),
	 L3=0,
	 maxiter=1e3)
   cat("optim.elnet.fmpr end\n")
   g <- groupridge4(X=Xtrain, Y=Ytrain,
	 lambda1=r$opt[1], lambda2=r$opt[2], lambda3=0,
	 maxiter=1e4, G=G0)
   g <- g[[1]][[1]][[1]]

   P <- scale(Xtest) %*% g
   res <- R2(as.numeric(P), as.numeric(center(Ytest)))

   cat("rep", rep, "R2 elasticnet.fmpr:", res, "\n")

   setwd(oldwd)

   list(
      R2=res,
      beta=g
   )
}

run.elnet.glmnet <- function(rep, dir=".", nfolds=10, r=25, Rthresh=0.5)
{
   oldwd <- getwd()
   setwd(dir)

   Xtrain <- as.matrix(read.table(sprintf("Xtrain_%s.txt", rep),
	    header=FALSE))
   Xtest <- as.matrix(read.table(sprintf("Xtest_%s.txt", rep), header=FALSE))
   Ytrain <- as.matrix(read.table(sprintf("Ytrain_%s.txt", rep), header=FALSE))
   Ytest <- as.matrix(read.table(sprintf("Ytest_%s.txt", rep), header=FALSE))

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)

   XtrainB <- scale(blockX(Xtrain, p, K))
   XtestB <- scale(blockX(Xtest, p, K))

   ytrain <- center(as.numeric(Ytrain))
   ytest <- center(as.numeric(Ytest))

   cat("optim.elnet.glmnet start\n")
   l <- maxlambda1(XtrainB, ytrain)
   r <- optim.elnet.glmnet(X=XtrainB, Y=ytrain, nfolds=nfolds,
	 L1=seq(l, l * 1e-3, length=r),
	 alpha=seq(0, 1, length=r))
   cat("optim.elnet.glmnet end\n")

   g <- glmnet(x=XtrainB, y=ytrain, lambda=r$opt[1], alpha=r$opt[2])
   g <- as.matrix(coef(g))[-1, ] # intercept should be almost zero

   P <- XtestB %*% g
   res <- R2(as.numeric(P), ytest)

   cat("rep", rep, "R2 elasticnet.glmnet:", res, "\n")

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
# type: same: same weights, same sparsity
#       sparsity: different weights, same sparsity
#       random: different weights and different sparsity
getB <- function(p, K, w, sparsity=0.8, type=NULL)
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
	 B <- matrix(rnorm(p * K), p, K)
	 s <- sample(0:1, p, TRUE, prob=c(sparsity, 1 - sparsity))
	 B[s == 0, ] <- 0
      }
   } else if(type == "random") {
      B <- matrix(rnorm(p * K), p, K)
      B <- B * sample(0:1, p * K, TRUE, prob=c(sparsity, 1 - sparsity))
   }

   B
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
   
   cat("Running inference\n")
   r.lasso <- lapply(1:nreps, run.lasso, dir=dir, r=grid, nfolds=nfolds)
   r.ridge <- lapply(1:nreps, run.ridge, dir=dir, r=grid, nfolds=nfolds)
   r.elnet.glmnet <- lapply(1:nreps, run.elnet.glmnet, dir=dir,
	 r=grid, nfolds=nfolds, Rthresh=setup$Rthresh)
   r.elnet.fmpr <- lapply(1:nreps, run.elnet.fmpr, dir=dir,
	 r=grid, nfolds=nfolds, Rthresh=setup$Rthresh)
   r.gr <- lapply(1:nreps, run.groupridge, dir=dir, r=grid, nfolds=nfolds,
   	 Rthresh=setup$Rthresh)
   #r.spg <- lapply(1:nreps, run.spg, dir=dir, r=grid, nfolds=nfolds)
   cat("Inference done\n")

   # Measure recovery of non-zeros
   recovery <- function(obj, dir)
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
   
   R2.gr <- sapply(r.gr, function(x) x$R2)
   R2.lasso <- sapply(r.lasso, function(x) x$R2)
   R2.ridge <- sapply(r.ridge, function(x) x$R2)
   R2.elnet.fmpr <- sapply(r.elnet.fmpr, function(x) x$R2)
   R2.elnet.glmnet <- sapply(r.elnet.glmnet, function(x) x$R2)
   #R2.spg <- sapply(r.spg, function(x) x$R2[1])
   
   R2.all <- cbind(
      GR=R2.gr,
      lasso=R2.lasso,
      ridge=R2.ridge,
      ElNetFMPR=R2.elnet.fmpr,
      ElNetGlmnet=R2.elnet.glmnet
      #SPG=R2.spg
   )

   rec.gr <- recovery(r.gr, dir)
   rec.lasso <- recovery(r.lasso, dir)
   rec.ridge <- recovery(r.ridge, dir)
   rec.elnet.fmpr <- recovery(r.elnet.fmpr, dir)
   rec.elnet.glmnet <- recovery(r.elnet.glmnet, dir)
   #rec.spg <- recovery(r.spg, dir)

   list(
      weights=list(
	 gr=r.gr,
	 lasso=r.lasso,
	 ridge=r.ridge,
	 elnet.fmpr=r.elnet.fmpr,
	 elnet.glmnet=r.elnet.glmnet
      ),
      recovery=list(
	 gr=rec.gr,
	 lasso=rec.lasso,
	 ridge=rec.ridge,
	 elnet.fmpr=rec.elnet.fmpr,
	 elnet.glmnet=rec.elnet.glmnet
	 #, spg=rec.spg),
      ),
      R2=R2.all
   )
}

# Remove ROC/PRC replications that had non-sensical results
clean.rocr <- function(obj)
{
   n <- length(obj@x.values)
   len <- sapply(obj@x.values, length)
   w <- len > 2
   obj@x.values <- obj@x.values[w]
   obj@y.values <- obj@y.values[w]
   obj@alpha.values <- obj@alpha.values[w]
   obj
}

