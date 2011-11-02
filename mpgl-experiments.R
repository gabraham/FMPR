
# Formulation of MPGL as group lasso

options(error=dump.frames)

library(ROCR)
library(ggplot2)

#set.seed(343439)

source("mpgl.R")
source("lasso.R")

p <- 200
C <- 10
N <- rep(100, C)
pop <- rep(1:C, N)
Bpop <- rep(1:C, each=p)
Bp <- rep(1:p, C)

lambda1 <- 5e-1

getB <- function(C)
{
   B <- cbind(unlist(lapply(1:C, function(i) {
      rnorm(p, i, 1)
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

# MPGL over the block diagonal data
B.mpgl <- mpgl(X, Y, p=p, C=C, lambda1=lambda1)

# Lasso for each population separately
B.lasso.sep <- unlist(lapply(1:C, function(i) {
   k <- s1[i]:(s1[i] + N[i] - 1)
   m <- s2[i]:(s2[i] + p - 1)
   lasso(X[k, m], Y[k, ], lambda1=lambda1)
}))

X.pool <- do.call(rbind, lapply(1:C, function(i) {
   k <- s1[i]:(s1[i] + N[i] - 1)
   m <- s2[i]:(s2[i] + p - 1)
   X[k, m]
}))

l <- lasso(X.pool, Y, lambda1=lambda1)
B.lasso.pool <- as.numeric(sapply(1:C, function(x) l))

pred <- prediction(
   labels=list(
      B != 0,
      B != 0,
      B != 0),
   predictions=list(
      Lasso=abs(B.lasso.sep),
      LassoPool=abs(B.lasso.pool),
      MPGL=abs(B.mpgl)
   )
)

perf.roc <- performance(pred, "tpr", "fpr")
perf.prc <- performance(pred, "prec", "rec")

d.roc <- data.frame(
   TPR=unlist(perf.roc@y.values),
   FPR=unlist(perf.roc@x.values),
   Method=rep(c("LassoSep", "LassoPool", "MPGL"), sapply(perf.roc@x.values, length))
)
g.roc <- ggplot(d.roc, aes(x=FPR, y=TPR, colour=Method)) + geom_line()

d.prc <- data.frame(
   Prec=unlist(perf.prc@y.values),
   Rec=unlist(perf.prc@x.values),
   Method=rep(c("LassoSep", "LassoPool", "MPGL"), sapply(perf.roc@x.values, length))
)
g.prc <- ggplot(d.prc, aes(x=Rec, y=Prec, colour=Method)) + geom_line()

print(g.roc)
print(g.prc)

