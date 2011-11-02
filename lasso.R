# plain l1 lasso for single outputs, very naive implementation
#
# Assumes X is standardised and y is zero-mean
lasso <- function(X, y, lambda1=0, maxiter=1e4L, eps=1e-4)
{
   p <- ncol(X)
   beta <- numeric(p)
   oldloss <- Inf
   loss <- Inf

   for(iter in 1:maxiter)
   {
      cat("iter", iter, "loss:", loss, "\n")   
      for(j in 1:p)
      {
	 Err <- X %*% beta - y
	 loss <- mean(Err^2)
	 s <- beta[j] - crossprod(X[,j], Err) / sum(X[,j]^2)
	 #beta[j] <- pmax(0, 1 - lambda1) * s
	 beta[j] <- sign(s) * max(abs(s) - lambda1, 0)
      }

      if(abs(oldloss - loss) < eps) {
	 cat("converged in", iter, "iterations\n")
	 break
      }

      oldloss <- loss
   }
   beta
}

