
# l1/l2
#
# lambda1: penalises the same weight across tasks
# lambda2: same, but shrinks
# step: # reduce Newton step size even further, don't know why but it works
#
l1l2 <- function(X, Y, lambda1=0, lambda2=0, maxiter=1e3, step=1e-2,
      eps=1e-3)
{
   W <- matrix(0, p, k)
   loss <- Inf
   oldloss <- 0
   
   nz <- numeric(p)
   
   for(iter in 1:maxiter)
   {
      E <- X %*% W - Y
      loss <- mean(E^2)
      cat(iter, loss, "\n")
      grad <- crossprod(X, E) / (N - 1) * step    
      for(j in 1:p)
      {
         S <- W[j,] - grad[j, ]
         normS <- sqrt(drop(crossprod(S)))
         if(normS > 0) {
            W[j,] <- pmax(0, 1 - lambda1 / normS) * S / (1 + lambda2)
         }
         nz[j] <- sum(W[j,] != 0)
      }
      cat(summary(nz), "\n")
      if(abs(loss - oldloss) <= eps) {
         break
      }
      oldloss <- loss
   }

   W
}


