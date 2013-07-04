# Code for simulating multi-task data

v <- strsplit(commandArgs(TRUE), "=")
for(m in v) {
   eval(parse(text=sprintf("%s=\"%s\"", m[1], m[2])))
}

if(!exists("idv") || idv == "" || is.na(idv)) {
   stop("idv not specified")
}

library(ROCR)
library(glmnet)
library(ggplot2)
library(FMPR)
library(doMC)
registerDoMC(cores=5)

options(error=dump.frames)

load("HAPGEN/sim_9/chr10.RData")
#load("HAPGEN/sim_9/chr10_pruned.RData")
if(!exists("geno")) {
   stop("'geno' object does not exist")
}

seed <- sample(1e6L, 1)
set.seed(seed)

################################################################################
# Configure the experiments here
setup <- list(
   # The reference setup
   Expr1=list(dir=c("Expr1"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="same")),
   
   # Different sample size
   Expr2=list(dir=c("Expr2"), N=50, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="same")),
   Expr3=list(dir=c("Expr3"), N=200, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="same")),

   # Different noise levels
   Expr4=list(dir=c("Expr4"), N=100, p=100, K=10, sigma=0.5,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="same")),
   Expr5=list(dir=c("Expr5"), N=100, p=100, K=10, sigma=2,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="same")),

   # Different number of tasks
   Expr6=list(dir=c("Expr6"), N=100, p=100, K=5, sigma=1,
      type="geno", 
      B=getB(p=100, K=5, w=0.1, type="same")),
   Expr7=list(dir=c("Expr7"), N=100, p=100, K=20, sigma=1,
      type="geno", 
      B=getB(p=100, K=20, w=0.1, type="same")),

   # Different weights
   Expr8=list(dir=c("Expr8"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=0.05, type="same")),
   Expr9=list(dir=c("Expr9"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=0.5, type="same")),

   # High dimensions, increasing p
   Expr10=list(dir=c("Expr10"), N=100, p=50, K=10, sigma=1,
      type="geno", 
      B=getB(p=50, K=10, w=0.1, type="same")),
   Expr11=list(dir=c("Expr11"), N=100, p=200, K=10, sigma=1,
      type="geno", 
      B=getB(p=200, K=10, w=0.1, type="same")),

   # Different weights across tasks, same sparsity
   Expr12=list(dir=c("Expr12"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="sparsity")),

   # Unrelated tasks
   Expr13=list(dir=c("Expr13"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="random")),

   # Same sparsity, different weight, mean=0
   Expr14=list(dir=c("Expr14"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="sparsity", mean=0, sd=0.1)),
   Expr15=list(dir=c("Expr15"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="sparsity", mean=0, sd=2)),
   Expr16=list(dir=c("Expr16"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0, sd=2)),

   # Same sparsity, different weight, mean=0.5
   Expr17=list(dir=c("Expr17"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.05)),
   Expr18=list(dir=c("Expr18"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.5)),
   Expr19=list(dir=c("Expr19"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=2)),

   Expr20=list(dir=c("Expr20"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.05)),

   # Negative correlations
   Expr21=list(dir=c("Expr21"), N=100, p=100, K=10, sigma=1, type="geno",
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="mixed")),

   # Some tasks are related (clusters), some aren't
   Expr22=list(dir=c("Expr22"), N=100, p=100, K=10, sigma=0.5,
      type="geno", 
      B=getB(p=100, K=10, w=0.1, type="cluster")),

   Expr23=list(dir=c("Expr23"), N=100, p=100, K=20, sigma=0.5,
      type="geno", 
      B=getB(p=100, K=20, w=0.1, type="clustersparse")),


   ##############################################################################
   # Independent X

   # The reference setup
   Expr1A=list(dir=c("Expr1A"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=getB(p=100, K=10, w=0.1, type="same")),
   
   # Different sample size
   Expr2A=list(dir=c("Expr2A"), N=50, p=100, K=10, sigma=1, type="artificial",
	 B=getB(p=100, K=10, w=0.1, type="same")),
   Expr3A=list(dir=c("Expr3A"), N=200, p=100, K=10, sigma=1, type="artificial",
	 B=getB(p=100, K=10, w=0.1, type="same")),

   # Different noise levels
   Expr4A=list(dir=c("Expr4A"), N=100, p=100, K=10, sigma=0.5, type="artificial",
	 B=getB(p=100, K=10, w=0.1, type="same")),
   Expr5A=list(dir=c("Expr5A"), N=100, p=100, K=10, sigma=2, type="artificial",
	 B=getB(p=100, K=10, w=0.1, type="same")),

   # Different number of tasks
   Expr6A=list(dir=c("Expr6A"), N=100, p=100, K=5, sigma=1, type="artificial",
	 B=getB(p=100, K=5, w=0.1, type="same")),
   Expr7A=list(dir=c("Expr7A"), N=100, p=100, K=20, sigma=1, type="artificial",
	 B=getB(p=100, K=20, w=0.1, type="same")),

   # Different weights
   Expr8A=list(dir=c("Expr8A"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=getB(p=100, K=10, w=0.05, type="same")),
   Expr9A=list(dir=c("Expr9A"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=getB(p=100, K=10, w=0.5, type="same")),

   # High dimensions, increasing p
   Expr10A=list(dir=c("Expr10A"), N=100, p=50, K=10, sigma=1, type="artificial",
	 B=getB(p=50, K=10, w=0.1, type="same")),
   Expr11A=list(dir=c("Expr11A"), N=100, p=200, K=10, sigma=1, type="artificial",
	 B=getB(p=200, K=10, w=0.1, type="same")),

   # Different weights across tasks, same sparsity
   Expr12A=list(dir=c("Expr12A"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=getB(p=100, K=10, w=0.1, type="sparsity")),

   # Unrelated tasks
   Expr13A=list(dir=c("Expr13A"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=getB(p=100, K=10, w=0.1, type="random")),

   # Same sparsity, different weight, mean=0
   Expr14A=list(dir=c("Expr14A"), N=100, p=100, K=10, sigma=1, type="artificial",
         B=getB(p=100, K=10, w=0.1, type="sparsity", mean=0, sd=0.5)),
   Expr15A=list(dir=c("Expr15A"), N=100, p=100, K=10, sigma=1, type="artificial",
         B=getB(p=100, K=10, w=0.1, type="sparsity", mean=0, sd=2)),
   Expr16A=list(dir=c("Expr16A"), N=100, p=100, K=10, sigma=1, type="artificial",
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0, sd=2)),

   # Same sparsity, different weight, mean=0.5
   Expr17A=list(dir=c("Expr17A"), N=100, p=100, K=10, sigma=1, type="artificial",
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.05)),
   Expr18A=list(dir=c("Expr18A"), N=100, p=100, K=10, sigma=1, type="artificial",
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.5)),
   Expr19A=list(dir=c("Expr19A"), N=100, p=100, K=10, sigma=1, type="artificial",
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=2)),

   Expr20A=list(dir=c("Expr20A"), N=100, p=100, K=10, sigma=1, type="artificial",
         B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.05)),

   # Negative correlations
   Expr21A=list(dir=c("Expr21A"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=getB(p=100, K=10, w=0.1, type="mixed")),

   # Some tasks are related (clusters), some aren't
   Expr22A=list(dir=c("Expr22A"), N=100, p=100, K=10, sigma=0.5, type="artificial",
	 B=getB(p=100, K=10, w=0.1, type="cluster")),

   Expr23A=list(dir=c("Expr23A"), N=100, p=100, K=20, sigma=0.5, type="artificial",
	 B=getB(p=100, K=20, w=0.1, type="clustersparse")),

   ExprTest1=list(dir=c("ExprTest1"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=abs(getB(p=100, K=10, sd=sqrt(0.10), type="sparsity"))),

   ExprTest2=list(dir=c("ExprTest2"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=abs(getB(p=100, K=10, sd=sqrt(0.15), type="sparsity"))),

   ExprTest3=list(dir=c("ExprTest3"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=abs(getB(p=100, K=10, sd=sqrt(0.05), type="sparsity"))),

   ExprTest4=list(dir=c("ExprTest4"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=abs(getB(p=100, K=10, sd=sqrt(0.25), type="sparsity"))),

   ExprTest5=list(dir=c("ExprTest5"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=abs(getB(p=100, K=10, sd=sqrt(0.5), type="sparsity"))),

   ExprTest6=list(dir=c("ExprTest6"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=abs(getB(p=100, K=10, sd=1, type="sparsity"))),

   ExprTest7=list(dir=c("ExprTest7"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=abs(getB(p=100, K=10, sd=2, type="sparsity"))),

   ExprTest8=list(dir=c("ExprTest8"), N=100, p=100, K=10, sigma=1, type="artificial",
	 B=abs(getB(p=100, K=10, sd=4, type="sparsity"))),

   Expr17pruned=list(dir=c("Expr17pruned"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=0.05)),

   Expr19pruned=list(dir=c("Expr19pruned"), N=100, p=100, K=10, sigma=1,
      type="geno", 
      B=getB(p=100, K=10, w=NULL, type="sparsity", mean=0.5, sd=2))

)

if(!exists("nreps", mode="integer")) {
   nreps <- 100
}
if(!exists("grid", mode="integer")) {
   grid <- 20
}
if(!exists("nfolds", mode="integer")) {
   nfolds <- 5
}

system.time({
   res <- lapply(setup[idv], run, nreps=nreps, grid=grid, nfolds=nfolds)
})
res.anl <- analyse(res[[1]], dir=setup[[idv]]$dir)
save(setup, res.anl, res, idv, nreps, grid, nfolds, seed,
   file=sprintf("%s/results_%s.RData", setup[[idv]]$dir, idv))

################################################################################

plot.exper(res.anl,
   models.rec=c("FMPR-w1", "FMPR-w2", "GFlasso-w1", "GFlasso-w2", "Lasso",
      "Ridge", "Elnet", "CCA", "PCA", "Univariable"),
   models.R2=c("FMPR-w1", "FMPR-w2", "GFlasso-w1",
      "GFlasso-w2", "Lasso", "Elnet", "Ridge")
)

