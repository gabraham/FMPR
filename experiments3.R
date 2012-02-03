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
library(FMPR)
library(doMC)
registerDoMC(cores=2)

options(error=dump.frames)

seed <- sample(1e6L, 1)
#seed <- 705168
set.seed(seed)


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
	 B=getB(p=50, K=10, w=0.2, type="same"), Rthresh=0.3),
   Expr15=list(dir=c("Expr15"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.3, type="same"), Rthresh=0.3),
   Expr16=list(dir=c("Expr16"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.4, type="same"), Rthresh=0.3),

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
   Expr24=list(dir=c("Expr24"), N=50, p=20, K=2, sigma=1,
	 B=getB(p=20, K=2, w=0.5, type="same", sparsity=0.5), Rthresh=0.3),

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


nreps <- 50
grid <- 10
nfolds <- 5

# Don't call SPG automatically because it's too slow
do.spg <- function()
{
   res.spg <- run.spg2(setup[[idv]], nreps=nreps, grid=grid, nfolds=nfolds)
   save(setup, res.spg, idv, file=sprintf("results_spg_%s.RData", idv))
   load(sprintf("results_%s.RData", idv))
   res[[1]]$weights <- c(res[[1]]$weights, res.spg$weights)
   res[[1]]$recovery <- c(res[[1]]$recovery, res.spg$recovery)
   res[[1]]$R2 <- cbind(res[[1]]$R2, res.spg$R2)
   res
}

system.time({
   res <- lapply(setup[idv], run, nreps=nreps, grid=grid, nfolds=nfolds)
})
save(setup, res, idv, file=sprintf("results_%s.RData", idv))
#res <- do.spg()
#save(setup, res, idv, file=sprintf("results_%s.RData", idv))

################################################################################

source("plotexper.R", echo=TRUE)


