
library(glmnet)
library(ROCR)
library(ggplot2)

source("mpgl.R")

options(error=dump.frames)

#d <- read.table("Data/sim.out.raw", header=TRUE, stringsAsFactors=FALSE)

x <- as.matrix(d[, -(1:6)])

x <- x[, 1:100]

N <- rep(nrow(x) / C, C)
p <- ncol(x)
C <- 10
pop <- rep(1:C, N)
Bpop <- rep(1:C, each=p)
Bp <- rep(1:p, C)

lambda1 <- 1e-2
L <- seq(lambda1, lambda1 / 100, length=20)

getB <- function(C)
{
   B <- cbind(unlist(lapply(1:C, function(i) {
      rnorm(p, i/2, 1)
   })))

   # make B group-sparse
   for(j in 1:p) {
      B[Bp == j,] <- B[Bp == j,] * sample(0:1, size=1, prob=c(0.9, 0.1))
   }
   B
}

B <- getB(C)
while(all(B == 0))
   B <- getB(C)

# set up block-diagonal X
X <- matrix(0, sum(N), p * C) 
s1 <- seq(1, sum(N), N[1])
s2 <- seq(1, p * C , p)
for(i in 1:C) {
   X[s1[i]:(s1[i] + N[i] - 1), s2[i]:(s2[i] + p - 1)] <- scale(rnorm(N[i] * p))
}

# zero-mean Y to remove intercept term
Y <- scale(X %*% B + rnorm(nrow(X) * ncol(B)), scale=FALSE)

g.pooled <- glmnet(X, drop(Y))
B.lasso.pooled <- as.matrix(coef(g.pooled))[-1,]

B.mpgl <- lapply(L, function(lambda1) mpgl(X, Y, p=p, C=C, lambda1=lambda1))
B.mpgl.mat <- do.call(cbind, B.mpgl)


stop("need to separate the models for each penalty")
pred <- prediction(
   labels=list(
      matrix(B != 0, byrow=FALSE, ncol=ncol(B.lasso.pooled),
	 nrow=nrow(B.lasso.pooled)),
      matrix(B != 0, byrow=FALSE, ncol=ncol(B.mpgl.mat),
	 nrow=nrow(B.mpgl.mat))
   ),
   predictions=list(
      LassoPool=abs(B.lasso.pooled),
      LassoSep=abs(B.mpgl.mat)
   )
)

perf.prc <- performance(pred, "prec", "rec")

d.prc <- data.frame(
   Prec=unlist(perf.prc@y.values),
   Rec=unlist(perf.prc@x.values),
   Method=rep(c("LassoPool", "MPGL"), sapply(perf.prc@x.values, length))
)
g.prc <- ggplot(d.prc, aes(x=Rec, y=Prec, colour=Method))
g.prc <- g.prc + geom_point(size=3) + geom_line()

#print(g.roc)
print(g.prc)


