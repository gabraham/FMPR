
#library(glmnet)
library(ROCR)
library(ggplot2)
library(MASS)

source("mpgl.R")
source("lasso.R")

options(error=dump.frames)

#d <- read.table("Data/sim.out.raw", header=TRUE, stringsAsFactors=FALSE)
#save(d, file="Data/sim.out.RData")
load("Data/sim.out.RData")

x <- as.matrix(d[, -(1:6)])

v <- apply(x, 2, var)
x <- x[, v >= 0.005]

#w <- seq(1, 2000, by=10)
w <- 1:300 + 1000
x <- x[, w]

# The number of samples in each population is nrow(x)/C
C <- 5
N <- rep(nrow(x) / C, C)
p <- ncol(x)
pop <- rep(1:C, N)
Bpop <- rep(1:C, each=p)
Bp <- rep(1:p, C)

lambda1 <- 1e-1
L <- seq(lambda1, lambda1 / 20, length=5)

getB <- function(C)
{
   #B <- cbind(unlist(lapply(1:C, function(i) {
   #   rnorm(p, i/2, 1)
   #   #runif(p, -5, 5)
   #})))

   B <- matrix(0, p * C, 1)

   # make B group-sparse
   for(j in 1:p) {
      
      # is this SNP non-zero in at least one population?
      z <- sample(c(TRUE, FALSE), size=1, prob=c(0.1, 0.9))
      if(z) {
	 # is the SNP non-zero for all or just some some of the populations?
	 s <- sample(0:1, size=C, prob=c(0.2, 0.8), replace=TRUE)
	 B[Bp == j,] <- rnorm(C, 0, 1) * s
      }
   }
   B
}

B <- getB(C)
while(all(B == 0))
   B <- getB(C)

# how many SNPs are zero across all populations, 1 population, 2 populations,
# etc.
table(sapply(1:p, function(j) sum(B[Bp == j, ] != 0)))

#stop()

#plot(B, col=Bpop)

# Pooling all populations together as one dataset, using the same genotype for
# each population.
X.pooled <- x

## set up block-diagonal X, we use the same genotypes for all populations
X.block <- matrix(0, sum(N), p * C) 
s1 <- seq(1, sum(N), N[1])
s2 <- seq(1, p * C , p)
for(i in 1:C) {
   j <- s1[i]:(s1[i] + N[i] - 1)
   k <- s2[i]:(s2[i] + p - 1)

   #s <- scale(x[j, ])
   #s[is.nan(s)] <- 0
   #X.block[j, k] <- s
   X.block[j, k] <- x[j, ]
}

X.sep <- lapply(1:C, function(i) {
   j <- s1[i]:(s1[i] + N[i] - 1)
   x[j, ]
})


# no intercept since Y is scaled
#Xs <- scale(X)

# zero-mean Y to remove intercept term
Y <- scale(X.block %*% B + rnorm(nrow(X.block) * ncol(B), 0, 1), scale=FALSE)

#g.pooled <- glmnet(X, drop(Y))
#B.lasso.pooled <- lapply(seq(along=g.pooled$lambda), function(i) {
#   abs(as.matrix(coef(g.pooled))[-1, i])
#})

B.mpgl <- lapply(L, function(lambda1) {
   abs(mpgl(X.block, Y, p=p, C=C, lambda1=lambda1))
})

B.lasso.pooled <- lapply(L, function(lambda) {
   abs(lasso.activeset(scale(X.pooled), drop(Y), lambda1=lambda))
})

B.lasso.sep <- lapply(L, function(lambda) {
   r <- lapply(1:C, function(i) {
      j <- s1[i]:(s1[i] + N[i] - 1)
      abs(lasso.activeset(x[j, ], Y[j, ], lambda1=lambda1))
   })
   unlist(r)
})


#B.ridge.pooled <- ginv(Xs) %*% Y

#B.mpgl.mat <- do.call(cbind, B.mpgl)

B.true <- lapply(1:length(L), function(i) {
      drop(B != 0)
})

pred.lasso.pooled <- prediction(
   labels=B.true,
   predictions=lapply(B.lasso.pooled, function(b) {
      do.call("c", lapply(1:C, function(i) b))
   })
)

pred.mpgl <- prediction(
   labels=B.true,
   predictions=B.mpgl
)

pred.lasso.sep <- prediction(
   labels=B.true,
   predictions=B.lasso.sep
)

perf.prc.lasso.pooled <- performance(pred.lasso.pooled, "prec", "rec")
perf.prc.mpgl <- performance(pred.mpgl, "prec", "rec")
perf.prc.lasso.sep <- performance(pred.lasso.sep, "prec", "rec")

d.prc.lasso.pooled <- data.frame(
   Prec=unlist(perf.prc.lasso.pooled@y.values),
   Rec=unlist(perf.prc.lasso.pooled@x.values),
   Method="LassoPool",
   Penalty=rep(L, sapply(perf.prc.lasso.pooled@y.values, length))
)

d.prc.mpgl <- data.frame(
   Prec=unlist(perf.prc.mpgl@y.values),
   Rec=unlist(perf.prc.mpgl@x.values),
   Method="MPGL",
   Penalty=rep(L, sapply(perf.prc.mpgl@y.values, length))
)

d.prc.lasso.sep <- data.frame(
   Prec=unlist(perf.prc.lasso.sep@y.values),
   Rec=unlist(perf.prc.lasso.sep@x.values),
   Method="LassoSep",
   Penalty=rep(L, sapply(perf.prc.lasso.sep@y.values, length))
)

d.prc <- rbind(d.prc.lasso.pooled, d.prc.mpgl, d.prc.lasso.sep)
g.prc <- ggplot(d.prc, 
      aes(x=Rec, y=Prec, shape=factor(Penalty), colour=Method))
g.prc <- g.prc + geom_point(size=4) + geom_line()

print(g.prc)

#perf.auc.lasso.pooled <- performance(pred.lasso.pooled, "sens", "spec")
#perf.auc.mpgl <- performance(pred.mpgl, "sens", "spec")
#
#d.auc.lasso.pooled <- data.frame(
#   Sensitivity=unlist(perf.auc.lasso.pooled@y.values),
#   Specificity=unlist(perf.auc.lasso.pooled@x.values),
#   Method="LassoPool",
#   Penalty=rep(L, sapply(perf.auc.lasso.pooled@y.values, length))
#)
#
#d.auc.mpgl <- data.frame(
#   Sensitivity=unlist(perf.auc.mpgl@y.values),
#   Specificity=unlist(perf.auc.mpgl@x.values),
#   Method="MPGL",
#   Penalty=rep(L, sapply(perf.auc.mpgl@y.values, length))
#)
#
#d.auc <- rbind(d.auc.lasso.pooled, d.auc.mpgl)
#g.auc <- ggplot(d.auc,
#      aes(x=Specificity, y=Sensitivity, shape=factor(Penalty), colour=Method))
#g.auc <- g.auc + geom_point(size=4) + geom_line()
#
#print(g.auc)



