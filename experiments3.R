# Code for simulating multi-task data
#
# Gad Abraham, 2012 (c)
#

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

options(error=dump.frames)

seed <- sample(1e6L, 1)
#seed <- 705168
set.seed(seed)

source("methods.R")
source("eval.R")
source("exprutil.R")

################################################################################
# Configure the experiments here
setup <- list(
   # Different sample size
   Expr1=list(dir=c("Expr1"), N=50, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr2=list(dir=c("Expr2"), N=100, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr3=list(dir=c("Expr3"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr4=list(dir=c("Expr4"), N=500, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),

   # Different noise levels
   Expr5=list(dir=c("Expr5"), N=200, p=50, K=10, sigma=0.01,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr6=list(dir=c("Expr6"), N=200, p=50, K=10, sigma=0.5,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr7=list(dir=c("Expr7"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr8=list(dir=c("Expr8"), N=200, p=50, K=10, sigma=2,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),

   # Different number of tasks
   Expr9=list(dir=c("Expr9"), N=200, p=50, K=2, sigma=1,
	 B=getB(p=50, K=2, w=0.5, type="same"), Rthresh=0.3),
   Expr10=list(dir=c("Expr10"), N=200, p=50, K=5, sigma=1,
	 B=getB(p=50, K=5, w=0.5, type="same"), Rthresh=0.3),
   Expr11=list(dir=c("Expr11"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr12=list(dir=c("Expr12"), N=200, p=50, K=20, sigma=1,
	 B=getB(p=50, K=20, w=0.5, type="same"), Rthresh=0.3),

   # Different weights
   Expr13=list(dir=c("Expr13"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.1, type="same"), Rthresh=0.3),
   Expr14=list(dir=c("Expr14"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr15=list(dir=c("Expr15"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.8, type="same"), Rthresh=0.3),
   Expr16=list(dir=c("Expr16"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=1, type="same"), Rthresh=0.3),

   # Different correlation thresholds
   Expr17=list(dir=c("Expr17"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.1),
   Expr18=list(dir=c("Expr18"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr19=list(dir=c("Expr19"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.7),
   Expr20=list(dir=c("Expr20"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.9),

   Expr21=list(dir=c("Expr21"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="same"), Rthresh=0.3),

   # Different weights across tasks, same sparsity
   Expr22=list(dir=c("Expr22"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="sparsity"), Rthresh=0.3),

   # Unrelated tasks
   Expr23=list(dir=c("Expr23"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="random"), Rthresh=0.3),

   # Small testing experiment
   Expr24=list(dir=c("Expr24"), N=100, p=30, K=2, sigma=1,
	 B=getB(p=30, K=2, w=0.5, type="same"), Rthresh=0.3),

   # High dimensions, increasing p
   Expr25=list(dir=c("Expr25"), N=100, p=100, K=10, sigma=1,
	 B=getB(p=100, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr26=list(dir=c("Expr26"), N=100, p=200, K=10, sigma=1,
	 B=getB(p=200, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr27=list(dir=c("Expr27"), N=100, p=500, K=10, sigma=1,
	 B=getB(p=500, K=10, w=0.5, type="same"), Rthresh=0.3),
   Expr28=list(dir=c("Expr28"), N=100, p=1000, K=10, sigma=1,
	 B=getB(p=1000, K=10, w=0.5, type="same"), Rthresh=0.3)
)

res <- lapply(setup[idv], run, nreps=50, grid=20, nfolds=5)
save(setup, res, idv, file=sprintf("results_%s.RData", idv))

################################################################################

