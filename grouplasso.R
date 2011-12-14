options(error=dump.frames)

library(MASS)

seed <- 93184091
set.seed(seed)

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

#
# X: block diagonal matrix of inputs
# Y: column vector of outputs
# p: number of variables in each population
# C: number of populations
# W: p times C matrix of weights for inclusion/exclusion of
#     each input from a group, used to define the groups
# lambda1: inter-group l1 penalty
# lambda2: inter-group l2 penalty
# lambda2g: intra-group l2 penalty
# 
# implements an L1 penalty across groups (to encourage sparsity at the group
# levels) and an L2 penalty inside groups (to shrink correlated outputs
# towards each other)
#
grouplasso.weighted <- function(X, Y, p, C, lambda1=0, lambda2=0, lambda2g=0,
      maxiter=1e4L, eps=1e-4)
{
   B <- cbind(numeric(p * C))
   s <- seq(0, p * C - 1, p)
   oldloss <- 0
   loss <- Inf
   
   for(iter in 1:maxiter)
   {
      cat("iter", iter, "loss:", loss, "B:", drop(crossprod(B)), "\n")

      # operates over the p variables in the C groups
      for(j in 1:p)
      {
	 Err <- X %*% B - Y
	 loss <- mean(Err^2)

	 # qr will only work if input matrix is invertible
	 # so only update variables that are not all zero
	 # TODO: move this check out of the loop
	 z <- apply(X[, s+j, drop=FALSE], 2, function(x) all(x == 0))

	 S <- numeric(length(s))
	 #w <- W[j, !z]
	 x <- X[, (s+j)[!z], drop=FALSE]
	 L <- lambda2g * diag(ncol(x))

	 # we need to solve B for each block, 
	 step <- qr.solve(crossprod(x) + L, crossprod(x, Err))
	 S[!z] <- B[(s + j)[!z], , drop=FALSE] - step
         normS <- sqrt(drop(crossprod(S)))
	 if(j == 1) cat("normS:", normS, "step:", drop(crossprod(step)), "\n")

         if(normS > .Machine$double.eps) {
	    B[s + j, ] <- pmax(0, 1 - lambda1 / normS) * S / (1 + lambda2)
         }
	 if(j == 1) browser()
      }

      #if(iter == 2000) browser()

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
   C <- 3
   X <- (matrix(rnorm(N * p), N, p))
   b <- rnorm(p) #* sample(0:1, p, replace=TRUE)
   # same B for all groups
   B <- matrix(rep(b, C), p, C)
   Y <- X %*% B + rnorm(N * C)

   Y2 <- cbind(as.numeric(Y))
   Xblock <- blockX(X, p, C)

   #B.gl <- grouplasso(Xblock, Y2, p, C) 
   #B.glm <- matrix(B.gl, p, C)
   #diag(cor(B, B.glm))
   #sum((B - B.glm)^2)
   #diag(cor(sign(B), sign(B.glm)))

   ## unpenalised least squares
   B.qr <- qr.solve(X, Y)
   #diag(cor(B, B.qr))
   #sum((B - B.qr)^2)
   #diag(cor(sign(B), sign(B.qr)))

   ## ridge regression separately on each output (yes it is!)
   #library(MASS)
   B.ridge <- ginv(X) %*% Y
   #diag(cor(B, B.ridge))
   #sum((B - B.ridge)^2)
   #diag(cor(sign(B), sign(B.ridge)))

   ## now with l1 penalisation
   #B.gl.p <- grouplasso(Xblock, Y2, p, C, lambda1=0.7) 
   #B.gl.pm <- matrix(B.gl.p, p, C)
   #diag(cor(B, B.gl.pm))
   #sum((B - B.gl.pm)^2)
   #diag(cor(sign(B), sign(B.gl.pm)))

   B.gl.w2 <- grouplasso.weighted(Xblock, Y2, p, C, lambda1=0, lambda2=0,
	 lambda2g=1e4) 
   B.gl.wm2 <- matrix(B.gl.w2, p, C)
   sum((B - B.gl.wm2)^2)
   #diag(cor(sign(B.gl.wm2), sign(B.gl.pm)))
   diag(cor(sign(B.gl.wm2), sign(B)))
   #diag(cor(sign(B.gl.wm2), sign(B.ridge)))

   summary(diag(tcrossprod(B.qr)))
   summary(diag(tcrossprod(B.ridge)))
   summary(diag(tcrossprod(B.gl.wm2)))

   plot(B.ridge, B.gl.wm2)


   ##W <- matrix(1, p, C)
   #B.gl.w <- grouplasso.weighted(Xblock, Y2, p, C, lambda1=0.7, lambda2=1e-1) 
   #B.gl.wm <- matrix(B.gl.w, p, C)
   #sum((B - B.gl.wm)^2)
   #diag(cor(sign(B.gl.wm), sign(B.gl.pm)))

#}
