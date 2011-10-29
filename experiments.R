
options(error=dump.frames)

source("lasso-l1l2.R")

library(glmnet)
library(ROCR)


N <- 200
p <- 2000
X <- scale(matrix(rnorm(N * p), N, p))
k <- 20
#z <- sample(0:1, k * p, TRUE, prob=c(0.99, 0.01))
#B <- matrix(rnorm(k * p) * z, p, k)
B <- matrix(rnorm(k * p), p, k)
for(j in 1:p) {
   B[j,] <- B[j,] * sample(0:1, 1, prob=c(0.9, 0.1))
}
nzB <- apply(B, 1, function(x) sum(x != 0))
Y <- scale(X %*% B + rnorm(N * k, 0, 1))

W1 <- l1l2(X, Y, lambda1=1e-3, lambda2=1e-3)
W2 <- l1l2(X, Y, lambda1=1e-4, lambda2=1e-3)
W3 <- l1l2(X, Y, lambda1=5e-3, lambda2=1e-3)

W.gl <- lapply(1:k, function(j) {
   g <- glmnet(X, Y[,j])
   as.matrix(coef(g))[-1,]
})

nz.gl <- lapply(W.gl, function(w) {
   apply(w, 2, function(x) sum(x != 0))
})


# Evaluate precision/recall
pred.l1l2.1 <- prediction(predictions=abs(W1[,1]), labels=B[,1] !=0) 
perf.l1l2.1 <- performance(pred.l1l2.1, "prec", "rec")

pred.l1l2.2 <- prediction(predictions=abs(W2[,1]), labels=B[,1] !=0) 
perf.l1l2.2 <- performance(pred.l1l2.2, "prec", "rec")

pred.l1l2.3 <- prediction(predictions=abs(W3[,1]), labels=B[,1] !=0) 
perf.l1l2.3 <- performance(pred.l1l2.3, "prec", "rec")

Bl <- lapply(1:k, function(j) matrix(B[,j] != 0, nrow=p, ncol=ncol(W.gl[[j]])))
#pred.gl <- prediction(predictions=lapply(W.gl, abs), labels=Bl)
#perf.gl <- performance(pred.gl, "prec", "rec")

pred.gl <- prediction(predictions=abs(W.gl[[1]]), labels=Bl[[1]])
perf.gl <- performance(pred.gl, "prec", "rec")

plot(perf.gl)
plot(perf.l1l2.1, add=TRUE, col=2, lwd=3)
plot(perf.l1l2.2, add=TRUE, col=3, lwd=3)
plot(perf.l1l2.3, add=TRUE, col=4, lwd=3)

