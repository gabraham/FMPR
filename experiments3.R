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

s <- sample(1e6L, 1)
set.seed(s)

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
   Expr22=list(dir=c("Expr20"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="sparsity"), Rthresh=0.3),

   # Unrelated tasks
   Expr23=list(dir=c("Expr20"), N=200, p=50, K=10, sigma=1,
	 B=getB(p=50, K=10, w=0.5, type="random"), Rthresh=0.3)
)

res <- lapply(setup[idv], run, nreps=10, grid=20, nfolds=5)
save(setup, res, idv, file=sprintf("results_%s.RData", idv))

################################################################################

pdf(sprintf("Expr%s.pdf", idv), width=12)

par(mfrow=c(1, 2))

plot(res[[1]]$recovery$gr$roc, avg="threshold", col=1, main="Partial ROC",
      xlim=c(0.97, 1), ylim=c(0.97, 1))
plot(res[[1]]$recovery$lasso$roc, avg="threshold", add=TRUE, col=2)
plot(res[[1]]$recovery$ridge$roc, avg="threshold", add=TRUE, col=3, lwd=3)
plot(res[[1]]$recovery$elnet.fmpr$roc, avg="threshold", add=TRUE, col=4, lwd=3)
plot(res[[1]]$recovery$elnet.glmnet$roc, avg="threshold", add=TRUE, col=5, lwd=3)

plot(res[[1]]$recovery$gr$prc, avg="threshold", col=1,
      main="Partial Precision-Recall", xlim=c(0.95, 1), ylim=c(0.95, 1))
plot(res[[1]]$recovery$lasso$prc, avg="threshold", col=2, add=TRUE)
plot(res[[1]]$recovery$ridge$prc, avg="threshold", add=TRUE, col=3)
plot(res[[1]]$recovery$elnet.fmpr$prc, avg="threshold", add=TRUE, col=4)
plot(res[[1]]$recovery$elnet.glmnet$prc, avg="threshold", add=TRUE, col=5)

dev.off()

t.test(res[[1]]$R2[, 1], res[[1]]$R2[,2])

mytheme <- function(base_size=10)
{
   structure(list(
	 axis.text.x=theme_text(size=20),
	 axis.text.y=theme_text(size=20, hjust=1),
	 axis.title.x=theme_text(size=23),
	 axis.title.y=theme_text(size=23),
	 axis.ticks=theme_blank(),
	 plot.title=theme_text(size=30),
	 legend.text=theme_text(size=20),
	 legend.title=theme_text(size=20, hjust=0),
	 legend.key.size=unit(2, "lines"),
	 legend.background=theme_rect(col=0, fill=0),
	 legend.key=theme_blank()
   ), class="options")
}

r2 <- melt(res[[1]]$R2)
colnames(r2) <- c("Replication", "Method", "R2")
g <- ggplot(r2, aes(Method, R2))
g <- g + geom_boxplot()
g <- g + scale_y_continuous(expression(R^2))
g <- g + theme_bw() + mytheme()

pdf(sprintf("Expr%s_R2.pdf", idv), width=12)
print(g)
dev.off()

