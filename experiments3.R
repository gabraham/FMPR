# Code for simulating multi-task data
#
# Gad Abraham, 2011 (c)
#

v <- strsplit(commandArgs(TRUE), "=")
for(m in v) {
   eval(parse(text=sprintf("%s=\"%s\"", m[1], m[2])))
}

if(!exists("id") || id == "" || is.na(id)) {
   stop("id not specified")
}

id <- as.integer(id)

library(ROCR)

options(error=dump.frames)

s <- sample(1e6L, 1)
set.seed(s)

source("methods.R")
source("eval.R")

makedata <- function(rep, dir=".", N=100, p=50, K=5, B, sigma=0.01)
{
   cat("rep", rep, "\n")

   Xtrain <- standardise(matrix(sample(0:2, N * p, replace=TRUE), N, p))
   Xtest <- standardise(matrix(sample(0:2, N * p, replace=TRUE), N, p))
   grp <- rep(1, K)
   XtrainB <- blockX(Xtrain, p, K)
   XtestB <- blockX(Xtest, p, K)

   XtrainBs <- standardise(XtrainB)
   XtestBs <- standardise(XtestB)
   
   write.table(B, file=sprintf("%s/B_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   
   Ytrain <- scale(Xtrain %*% B + rnorm(N * K, 0, sigma), scale=FALSE)
   Ytest <- scale(Xtest %*% B + rnorm(N * K, 0, sigma), scale=FALSE)
   ytrain <- as.numeric(Ytrain)
   ytest <- as.numeric(Ytest)
   #Ynet <- outer(grp, grp, function(x, y) as.integer(x == y))
   #lowerTriangle(Ynet) <- 0  
   #diag(Ynet) <- 0

   write.table(Xtrain, file=sprintf("%s/Xtrain_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(Ytrain, file=sprintf("%s/Ytrain_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   #write.table(Ynet, file=sprintf("%s/Ynetwork_%s.txt", dir, rep),
   #      col.names=FALSE, row.names=FALSE, quote=FALSE)
   write.table(XtrainBs, file=sprintf("%s/XtrainB_%s.txt", dir, rep),
         col.names=FALSE, row.names=FALSE, quote=FALSE)
   #write.table(cbind(grp), file=sprintf("%s/grp_%s.txt", dir, rep),
   #	 col.names=FALSE, row.names=FALSE, quote=FALSE)

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

   r <- optim.lasso(X=XtrainB, Y=ytrain, nfolds=10, grid=r, maxiter=1e5)
   g <- lasso3(X=XtrainB, y=ytrain, lambda1=r$opt[1])
   P <- XtestB %*% g
   res <- R2(as.numeric(P), as.numeric(ytest))
   cat("R2 lasso:", res, "\n")

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
   XtrainB <- as.matrix(read.table(sprintf("XtrainB_%s.txt", rep),
	    header=FALSE))
   XtestB <- as.matrix(read.table(sprintf("XtestB_%s.txt", rep),
	    header=FALSE))

   ytrain <- as.numeric(Ytrain)
   ytest <- as.numeric(Ytest)

   N <- nrow(Xtrain)
   K <- ncol(Ytrain)
   p <- ncol(Xtrain)
   
   r <- optim.ridge(X=Xtrain, Y=Ytrain, nfolds=10, grid=r)
   g <- ridge(XtrainB, ytrain, lambda2=r$opt[1])
   P <- XtestB %*% g
   res <- R2(as.numeric(P), as.numeric(ytest))
   cat("R2 ridge:", res, "\n")

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

   r <- optim.groupridge(X=Xtrain, Y=Ytrain, G=G, nfolds=10, grid=r, maxiter=1e5)
   g <- groupridge4(X=Xtrain, Y=Ytrain,
      lambda1=r$opt[1], lambda2=r$opt[2], lambda3=r$opt[3], maxiter=1e4, G=G)
   P <- Xtest %*% g
   res <- R2(as.numeric(P), as.numeric(Ytest))

   cat("R2 groupridge:", res, "\n")

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
   }
   m <- sapply(1:nreps, makedata, dir=dir,
	 N=setup$N, p=setup$p, K=setup$K, B=setup$B, sigma=setup$sigma)
   cat("Simulation done\n")
   
   cat("Running inference\n")
   r.gr <- lapply(1:nreps, run.groupridge, dir=dir, r=grid, nfolds=nfolds)
   r.lasso <- lapply(1:nreps, run.lasso, dir=dir, r=grid, nfolds=nfolds)
   #r.spg <- lapply(1:nreps, run.spg, dir=dir, r=grid, nfolds=nfolds)
   r.ridge <- lapply(1:nreps, run.ridge, dir=dir, r=grid, nfolds=nfolds)
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
      list(roc=roc, prc=prc)
   }
   
   R2.gr <- sapply(r.gr, function(x) x$R2)

   R2.lasso <- sapply(r.lasso, function(x) x$R2)
   R2.ridge <- sapply(r.ridge, function(x) x$R2)
   #R2.spg <- sapply(r.spg, function(x) x$R2[1])
   
   R2.all <- cbind(GR=R2.gr, lasso=R2.lasso, ridge=R2.ridge)#, SPG=R2.spg)
   
   rec.gr <- recovery(r.gr, dir)

   rec.lasso <- recovery(r.lasso, dir)
   rec.ridge <- recovery(r.ridge, dir)
   #rec.spg <- recovery(r.spg, dir)

   list(
      res=list(
	 gr=r.gr,
	 lasso=r.lasso,
	 ridge=r.ridge
      ),
      recovery=list(
	 gr=rec.gr,
	 lasso=rec.lasso,
	 ridge=rec.ridge#, spg=rec.spg),
      ),
      R2=R2.all
   )
}

################################################################################
# Configure the experiments here
setup <- list(
   # Different sample size
   Expr1=list(dir=c("Expr1"), N=50, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.5, type="same")),
   Expr2=list(dir=c("Expr2"), N=100, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.5, type="same")),
   Expr3=list(dir=c("Expr3"), N=200, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.5, type="same")),
   Expr4=list(dir=c("Expr4"), N=500, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.5, type="same")),
   # Different noise levels
   Expr5=list(dir=c("Expr5"), N=200, p=50, K=10, sigma=0.01,
	 B=getB(p=50, K=10, w=0.5, type="same")),
   Expr6=list(dir=c("Expr6"), N=200, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.5, type="same")),
   Expr7=list(dir=c("Expr7"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=51, K=10, w=0.5, type="same")),
   Expr8=list(dir=c("Expr8"), N=200, p=50, K=10, sigma=2,
	 B=getB(p=50, K=10, w=0.5, type="same")),
   # Different number of tasks
   Expr9=list(dir=c("Expr9"), N=200, p=50, K=2, sigma=0.5,
	 B=getB(p=50, K=2, w=0.5, type="same")),
   Expr10=list(dir=c("Expr10"), N=200, p=50, K=5, sigma=0.5,
	 B=getB(p=50, K=5, w=0.5, type="same")),
   Expr11=list(dir=c("Expr11"), N=200, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.5, type="same")),
   Expr12=list(dir=c("Expr12"), N=200, p=50, K=20, sigma=0.5,
	 B=getB(p=50, K=20, w=0.5, type="same")),
   # Different weights
   Expr13=list(dir=c("Expr13"), N=200, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.1, type="same")),
   Expr14=list(dir=c("Expr14"), N=200, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.5, type="same")),
   Expr15=list(dir=c("Expr15"), N=200, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.8, type="same")),
   Expr16=list(dir=c("Expr16"), N=200, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=1, type="same"))
)

res <- lapply(setup[id], run, nreps=30, grid=20, nfolds=10)
save(setup, res, id, file=sprintf("results_%s.RData", id))

pdf(sprint("Expr_%s.pdf", id), width=12)
par(mfrow=c(1, 2))
plot(res[[1]]$recovery$gr$roc, avg="threshold", col=1, main="ROC")
plot(res[[1]]$recovery$lasso$roc, avg="threshold", add=TRUE, col=2)
plot(res[[1]]$recovery$ridge$roc, avg="threshold", add=TRUE, col=3)
plot(res[[1]]$recovery$gr$prc, avg="threshold", col=1, main="Precision-Recall")
plot(res[[1]]$recovery$lasso$prc, avg="threshold", add=TRUE, col=2)
plot(res[[1]]$recovery$ridge$prc, avg="threshold", add=TRUE, col=3)
dev.off()

