
library(MASS)


#seed <- sample(1e9, 1)
#set.seed(seed)
set.seed(4398790)

N <- 100
p <- 20
K <- 3

## save to file
#write.table(X, file="X.txt", col.names=FALSE, row.names=FALSE)
#write.table(Y, file="Y.txt", col.names=FALSE, row.names=FALSE)
#
#adj <- outer(g, g, function(X, Y) { as.integer(X == Y) })
#write.table(adj, file="Ynetwork.txt", col.names=FALSE, row.names=FALSE)


# g: a vector of length K partitioned into M values
#groupridge <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g)
#{
#   N <- nrow(X)
#   p <- ncol(X)
#   K <- ncol(Y)
#   B1 <- matrix(0, p, K)
#   colnames(B1) <- g
#   S <- table(g)
#   loss <- 0
#   oldloss <- Inf
#   active <- matrix(1, p, K)
#   oldactive <- matrix(1, p, K)
#
#   v <- diag(crossprod(X))
#
#   for(iter in 1:5000)
#   {
#      cat("iter:", iter, loss, "\n")
#      for(k in 1:K)
#      {
#         for(j in 1:p)
#         {
#            d1 <- crossprod(X[, j], X %*% B1[, k] - Y[, k])
#            d2 <- v[j]
#	    
#	    # soft thresholding
#	    if(lambda1 > 0) 
#	    {
#	       s <- B1[j, k] - d1 / d2
#	       s2 <- sign(s) * max(abs(s) - lambda1, 0)
#
#	       if(s2 == 0)
#	       {
#		  B1[j, k] <- 0
#		  active[j, k] <- 0
#		  break
#	       }
#	    
#	       # l1 shrinkage when weight isn't zero
#	       d1 <- d1 + lambda1 * sign(B1[j, k])
#	    }
#   
#	    # Penalise by the differences from all other variables
#   	    # in the same partition, don't penalise an empty partition
#   	    w <- which(g == g[k])
#   	    d1 <- d1 + lambda2 * sum(B1[j, k] - B1[j, w])
#   	    d2 <- d2 + lambda2 * N
#   
#   	    # standard ridge penalty
#   	    d1 <- d1 + lambda3 * B1[j, k]
#   	    d2 <- d2 + lambda3
#   
#            B1[j, k] <- B1[j, k] - d1 / d2
#         }
#      }
#
#
#      loss <- mean((X %*% B1 - Y)^2)
#      if(all(active == oldactive) && abs(oldloss - loss) < 1e-3)
#	 break
#      oldloss <- loss
#      oldactive <- active
#   }
#   
#   B1
#}

groupridge <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, g)
{
   if(length(lambda1) == 1)
      lambda1 <- rep(lambda1, K)
   
   if(length(lambda2) == 1)
      lambda2 <- rep(lambda2, K)

   if(length(lambda3) == 1)
      lambda3 <- rep(lambda3, K)

   r <- .C("groupridge", as.numeric(X), as.numeric(Y), 
      numeric(p * K), nrow(X), ncol(X), ncol(Y),
      as.numeric(lambda1), as.numeric(lambda2), as.numeric(lambda3),
      as.integer(g))
   matrix(r[[3]], p, K)
}

maxlambda1 <- function(X, Y)
{
   r <- .C("maxlambda1", as.numeric(X), as.numeric(Y),
      numeric(K), nrow(X), ncol(X), ncol(Y))
   r[[3]]
}

blockX <- function(X, p, C)
{
   N <- nrow(X)
   
   Xblock <- matrix(0, N * C, p * C) 
   s1 <- seq(1, N * C, N)
   s2 <- seq(1, p * C, p)

   for(i in 1:C)
   {
      j <- s1[i]:(s1[i] + N - 1)
      k <- s2[i]:(s2[i] + p - 1)
      Xblock[j, k] <- X
   }

   Xblock
}

#L <- 10^seq(-4, 2, length=10)

dyn.load("groupridge.so")
X <- scale(matrix(rnorm(N * p), N, p))
g <- sample(3, size=K, replace=TRUE)
#B <- sapply(g, function(k) rep(k, p))
B <- matrix(sample(0:1, p * K, TRUE), p, K)
Y <- scale(X %*% B + rnorm(N, 0, 1), scale=FALSE)
B0 <- ginv(X) %*% Y

maxL <- maxlambda1(X, Y)
B1 <- groupridge(X, Y, lambda1=0, lambda2=0, lambda3=0, g=g)
mean((B0 - B1)^2)
B1.1 <- groupridge(X, Y, lambda1=maxL, lambda2=0, lambda3=0, g=g)
B1.2 <- groupridge(X, Y, lambda1=maxL * 0.99, lambda2=0, lambda3=0, g=g)
B1.3 <- groupridge(X, Y, lambda1=maxL * 0.9, lambda2=0, lambda3=0, g=g)
B1.4 <- groupridge(X, Y, lambda1=maxL * 0.8, lambda2=0, lambda3=0, g=g)

stop()

#B1.2 <- groupridge(X, Y, lambda1=0.8, lambda2=0, lambda3=0, g=g)
#B1.3 <- groupridge(X, Y, lambda1=0.7, lambda2=0, lambda3=0, g=g)
#B1.4 <- groupridge(X, Y, lambda1=0.6, lambda2=0, lambda3=0, g=g)
#B1.5 <- groupridge(X, Y, lambda1=0.5, lambda2=0, lambda3=0, g=g)
#B1.6 <- groupridge(X, Y, lambda1=0.1, lambda2=0, lambda3=0, g=g)
#B1.7 <- groupridge(X, Y, lambda1=0.01, lambda2=0, lambda3=0, g=g)

stop()


run <- function()
{
   X <- scale(matrix(rnorm(N * p), N, p))
   Xtest <- scale(matrix(rnorm(N * p), N, p))
   
   # put tasks in groups
   g <- sample(3, size=K, replace=TRUE)
   
   # same weights across all variables in one task for convenience
   B <- sapply(g, function(k) rep(k, p))
   
   # shared sparsity pattern for same variable across all tasks
   for(j in 1:p)
      B[j, ] <- B[j, ] * sample(0:1, 1)
   Y <- scale(X %*% B + rnorm(N, 0, 1), scale=FALSE)
   Ytest <- scale(Xtest %*% B + rnorm(N, 0, 1), scale=FALSE)

   B0 <- ginv(X) %*% Y

   B1 <- lapply(L, function(l1) {
      cat("lambda:", l1, "\n")
      lapply(L, function(l2) {
	 groupridge(X, Y, lambda1=0, lambda2=l1, lambda3=l2, g=g)
      })
   })

   B1 <- unlist(B1, recursive=FALSE)

   # standard ridge regression over each task separately
   B2 <- lapply(L, function(l) {
      qr.solve(crossprod(X) + l * diag(p), crossprod(X, Y))
   })
   
   # ridge regression, combining the tasks into one task
   Yv <- as.numeric(Y)
   Xblock <- blockX(X, p, K)
   B3 <- lapply(L, function(l) {
      b <- qr.solve(crossprod(Xblock) + l * diag(ncol(Xblock)),
   	 crossprod(Xblock, Yv))
      matrix(b, p, K)
   })
   
   # training
   loss0 <- colSums((Y - X %*% B0)^2)
   loss1 <- sapply(B1, function(b) colSums((Y - X %*% b)^2))
   loss2 <- sapply(B2, function(b) colSums((Y - X %*% b)^2))
   loss3 <- sapply(B3, function(b) colSums((Y - X %*% b)^2))
   
   R2.0 <- 1 - loss0 / colSums(sweep(Y, 2, colMeans(Y))^2)
   R2.1 <- rbind(apply(rbind(loss1), 2, function(l) {
      1 - l / colSums(sweep(Y, 2, colMeans(Y))^2)
   }))
   R2.2 <- rbind(apply(rbind(loss2), 2, function(l) {
      1 - l / colSums(sweep(Y, 2, colMeans(Y))^2)
   }))
   R2.3 <- rbind(apply(rbind(loss3), 2, function(l) {
      1 - l / colSums(sweep(Y, 2, colMeans(Y))^2)
   }))
   
   # test
   loss0t <- colSums((Ytest - Xtest %*% B0)^2)
   loss1t <- sapply(B1, function(b) colSums((Ytest - Xtest %*% b)^2))
   loss2t <- sapply(B2, function(b) colSums((Ytest - Xtest %*% b)^2))
   loss3t <- sapply(B3, function(b) colSums((Ytest - Xtest %*% b)^2))
   
   R2.0t <- 1 - loss0t / colSums(sweep(Ytest, 2, colMeans(Ytest))^2)
   R2.1t <- rbind(apply(rbind(loss1t), 2, function(l) {
      1 - l / colSums(sweep(Ytest, 2, colMeans(Ytest))^2)
   }))
   R2.2t <- rbind(apply(rbind(loss2t), 2, function(l) {
      1 - l / colSums(sweep(Ytest, 2, colMeans(Ytest))^2)
   }))
   R2.3t <- rbind(apply(rbind(loss3t), 2, function(l) {
      1 - l / colSums(sweep(Ytest, 2, colMeans(Ytest))^2)
   }))

   list(mean(R2.0t), colMeans(R2.1t), colMeans(R2.2t), colMeans(R2.3t))
}

Nrep <- 10
res <- replicate(Nrep, run(), simplify=FALSE)
r0 <- sapply(res, function(x) x[[1]])
r1 <- sapply(res, function(x) x[[2]])
r2 <- sapply(res, function(x) x[[3]])
r3 <- sapply(res, function(x) x[[4]])

r0m <- mean(r0)
r1m <- matrix(rowMeans(r1), length(L))
r2m <- rowMeans(r2)
r3m <- rowMeans(r3)

rowSD <- function(x) apply(x, 1, sd)

r0s <- sd(r0) / sqrt(Nrep)
r1s <- matrix(rowSD(r1) / sqrt(Nrep), length(L))
r2s <- rowSD(r2) / sqrt(Nrep)
r3s <- rowSD(r3) / sqrt(Nrep)

r <- cbind(r1m, r2m)

pdf("groupridge.pdf")
# Average R^2 across the K outputs
matplot(L, r, log="x",
      type="b", main="test", ylab=expression(R^2), lty=1,
      lwd=2,
      pch=c(rep(20, length(L), 21)),
      col=c(rep(1, length(L)), 2))
dev.off()
