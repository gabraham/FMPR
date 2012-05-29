
library(FMPR)
library(ROCR)
library(doMC)
registerDoMC(cores=2)

N <- 100
p <- 100
K <- 10 
w <- 0.1
B <- getB(p=p, K=K, w=w, type="same", sparsity=0.8)
print(table(sign(B)))

X <- scale(matrix(rnorm(N * p), N, p))
Y <- scale(X %*% B + rnorm(N * K))

G <- graph.sqr(cor(Y))
C <- gennetwork(Y, corthresh=0, cortype=2)
  
l <- max(maxlambda1(X, Y))
ngrid <- 30
lambda1 <- 2^seq(-0.01, -6, length=ngrid) * l
lambda3 <- 10^seq(-5, 5, length=ngrid)
   
g1 <- fmpr(X=X, Y=Y, G=G, lambda1=lambda1, lambda2=0, lambda3=lambda3)
g2 <- spg(X=X, Y=Y, C=C, lambda=lambda1, gamma=lambda3, maxiter=1e5)
g3 <- fmpr(X=X, Y=Y, lambda1=lambda1, lambda2=0, lambda3=0)

g1l <- unlist(unlist(g1, recursive=FALSE), recursive=FALSE)
g2l <- unlist(g2, recursive=FALSE)
g3l <- unlist(unlist(g3, recursive=FALSE), recursive=FALSE)

r1 <- sapply(g1l, function(B2) {
   mean((X %*% B2 - Y)^2)
})

r2 <- sapply(g2l, function(B2) {
   mean((X %*% B2 - Y)^2)
})

r3 <- sapply(g3l, function(B2) {
   mean((X %*% B2 - Y)^2)
})

summary(r1)
summary(r2)
summary(r3)

auc1 <- sapply(g1l, function(B2) {
   performance(prediction(
      labels=as.numeric(B != 0),
      predictions=as.numeric(abs(B2))
      ), "auc"
   )@y.values[[1]]
})

auc2 <- sapply(g2l, function(B2) {
   performance(prediction(
      labels=as.numeric(B != 0),
      predictions=as.numeric(abs(B2))
      ), "auc"
   )@y.values[[1]]
})

auc3 <- sapply(g3l, function(B2) {
   performance(prediction(
      labels=as.numeric(B != 0),
      predictions=as.numeric(abs(B2))
      ), "auc"
   )@y.values[[1]]
})

t.test(auc1, auc2)
t.test(auc1, auc3)

