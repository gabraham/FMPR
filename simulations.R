# Code for simulating multi-task data

v <- strsplit(commandArgs(TRUE), "=")
for(m in v) {
   eval(parse(text=sprintf("%s=\"%s\"", m[1], m[2])))
}

if(!exists("idv") || idv == "" || is.na(idv)) {
   stop("idv not specified")
}

idv <- as.integer(idv)

library(ROCR)
library(glmnet)
library(ggplot2)
library(FMPR)
library(doMC)
registerDoMC(cores=2)

options(error=dump.frames)

seed <- sample(1e6L, 1)
set.seed(seed)

################################################################################
# Configure the experiments here
setup <- list(
   # The reference setup
   Expr1=list(dir=c("Expr1"), N=100, p=100, K=10, sigma=1,
	 B=getB(p=100, K=10, w=0.1, type="same")),
   
   # Different sample size
   Expr2=list(dir=c("Expr2"), N=50, p=100, K=10, sigma=1,
	 B=getB(p=100, K=10, w=0.1, type="same")),
   Expr3=list(dir=c("Expr3"), N=200, p=100, K=10, sigma=1,
	 B=getB(p=100, K=10, w=0.1, type="same")),

   # Different noise levels
   Expr4=list(dir=c("Expr4"), N=100, p=100, K=10, sigma=0.5,
	 B=getB(p=100, K=10, w=0.1, type="same")),
   Expr5=list(dir=c("Expr5"), N=100, p=100, K=10, sigma=2,
	 B=getB(p=100, K=10, w=0.1, type="same")),

   # Different number of tasks
   Expr6=list(dir=c("Expr6"), N=100, p=100, K=5, sigma=1,
	 B=getB(p=100, K=5, w=0.1, type="same")),
   Expr7=list(dir=c("Expr7"), N=100, p=100, K=20, sigma=1,
	 B=getB(p=100, K=20, w=0.1, type="same")),

   # Different weights
   Expr8=list(dir=c("Expr8"), N=100, p=100, K=10, sigma=1,
	 B=getB(p=100, K=10, w=0.05, type="same")),
   Expr9=list(dir=c("Expr9"), N=100, p=100, K=10, sigma=1,
	 B=getB(p=100, K=10, w=0.5, type="same")),

   # High dimensions, increasing p
   Expr10=list(dir=c("Expr10"), N=100, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.1, type="same")),
   Expr11=list(dir=c("Expr11"), N=100, p=200, K=10, sigma=1,
	 B=getB(p=200, K=10, w=0.1, type="same")),

   # Different weights across tasks, same sparsity
   Expr12=list(dir=c("Expr12"), N=100, p=100, K=10, sigma=1,
	 B=getB(p=100, K=10, w=0.1, type="sparsity")),

   # Unrelated tasks
   Expr13=list(dir=c("Expr13"), N=100, p=100, K=10, sigma=1,
	 B=getB(p=100, K=10, w=0.5, type="random")),

   # Same sparsity, different weight, mean=0
   Expr14=list(dir=c("Expr14"), N=100, p=100, K=10, sigma=1,
         B=getB(p=100, K=10, w=0.1, type="sparsity", mean=0, sd=0.5)),
   Expr15=list(dir=c("Expr15"), N=100, p=100, K=10, sigma=1,
         B=getB(p=100, K=10, w=0.1, type="sparsity", mean=0, sd=2)),
   Expr16=list(dir=c("Expr16"), N=100, p=100, K=10, sigma=1,
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0, sd=2)),

   # Same sparsity, different weight, mean=0.5
   Expr17=list(dir=c("Expr17"), N=100, p=100, K=10, sigma=1,
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.05)),
   Expr18=list(dir=c("Expr18"), N=100, p=100, K=10, sigma=1,
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.5)),
   Expr19=list(dir=c("Expr19"), N=100, p=100, K=10, sigma=1,
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=2)),

   Expr20=list(dir=c("Expr20"), N=100, p=100, K=10, sigma=1,
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.05)),

   # Negative correlations
   Expr21=list(dir=c("Expr21"), N=100, p=100, K=10, sigma=1,
	 B=getB(p=100, K=10, w=0.1, type="mixed")),

   # Some tasks are related (clusters), some aren't
   Expr22=list(dir=c("Expr22"), N=100, p=100, K=10, sigma=0.5,
	 B=getB(p=100, K=10, w=0.1, type="cluster"))

)

nreps <- 10
grid <- 25
nfolds <- 5

system.time({
   res <- lapply(setup[idv], run, nreps=nreps, grid=grid, nfolds=nfolds)
})
save(setup, res, idv, nreps, grid, nfolds,
   file=sprintf("%s/results_%s.RData", setup[[idv]]$dir, idv))

################################################################################

lapply(res, plot)

