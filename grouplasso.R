options(error=dump.frames)

# Convert the input matrix X to a block diagonal matrix, where each block is
# identical to X
# 
# X: N by p matrix
# p: number of variables (same in all groups)
# C: number of groups
#
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

#
# X: block diagonal matrix of inputs
# Y: column vector of outputs
# p: number of variables in each population
# C: number of populations
#
grouplasso <- function(X, Y, p, C, lambda1=0, maxiter=1e4L, eps=1e-4)
{
   B <- cbind(numeric(p * C))
   s <- seq(0, p * C - 1, p)
   oldloss <- 0
   loss <- Inf
   
   for(iter in 1:maxiter)
   {
      cat("iter", iter, "loss:", loss, "\n")

      # operates over the p groups in the C populations
      for(j in 1:p)
      {
	 Err <- X %*% B - Y
	 loss <- mean(Err^2)

	 # qr will only work if input matrix is invertible
	 w <- apply(X[, s+j], 2, function(x) all(x == 0))
	 S <- numeric(length(s))
	 S[!w] <- B[(s + j)[!w], ] - qr.solve(X[, (s+j)[!w]], Err)
         normS <- sqrt(drop(crossprod(S)))
         if(normS > 0) {
	    B[s + j, ] <- pmax(0, 1 - lambda1 / normS) * S
         }
      }
   
      if(abs(loss - oldloss) <= eps) {
	 cat("converged in", iter, "iterations\n")
         break
      }
      oldloss <- loss
   }
   
   B
}

#test.grouplasso <- function()
#{
   N <- 70
   p <- 20
   C <- 5
   X <- matrix(rnorm(N * p), N, p)
   b <- rnorm(p)
   # same B for all groups
   B <- matrix(rep(b, C), p, C)
   Y <- X %*% B + rnorm(N * C)

   Y2 <- cbind(as.numeric(Y))
   Xblock <- blockX(X, p, C)
   B.gl <- grouplasso(Xblock, Y2, p, C) 
   B.gl.2 <- matrix(B.gl, p, C)
   diag(cor(B, B.gl.2))
   sum((B - B.gl.2)^2)

   # unpenalised least squares
   B.qr <- qr.solve(X, Y)
   diag(cor(B, B.qr))
   sum((B - B.qr)^2)

   # ridge regression
   library(MASS)
   B.ridge <- ginv(X) %*% Y
   diag(cor(B, B.ridge))
   sum((B - B.ridge)^2)
#}
