
# Multi-population group lasso
#
# X: block diagonal matrix of inputs
# Y: column vector of outputs
# p: number of variables in each population
# C: number of populations
#
mpgl <- function(X, Y, p, C, lambda1=0, maxiter=1e4L, eps=1e-4)
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
	 S <- B[s + j, ] - qr.solve(X[, s + j], Err)
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

