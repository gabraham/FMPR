# specnorm(X) is equivalent to MATLAB eigs(X'X, 1)
specnorm <- function(x)
{
   safe.svd(x, nu=0, nv=1)$d[1]^2
}

soft_threshold <- function(a, b)
{
   sign(a) * pmax(abs(a) - b, 0)
}

hard_threshold <- function(a, b=1)
{
   sign(a) * pmin(abs(a), b)
}

spgL1 <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, tol=1e-4, mu=1e-4)
{
   N <- nrow(X)
   J <- ncol(X)
   K <- ncol(Y)

   loss <- numeric(maxiter)
   B <- matrix(0, J, K)
   W <- B
   theta <- 1

   XX <- crossprod(X)
   XY <- crossprod(X, Y)
   CNorm <- 2 * max(colSums(C^2))

   C <- gamma * C
   L <- specnorm(X) + gamma^2 * CNorm / mu

   for(iter in 1:maxiter)
   {
      A <- hard_threshold(tcrossprod(C, W / mu), 1)
      grad <- XX %*% W - XY + crossprod(A, C)
      V <- W - (1 / L) * grad
      B_new <- soft_threshold(V, lambda / L)
      theta_new <- (sqrt(theta^4 + 4 * theta^2) - theta^2) / 2
      W <- B_new + (1 - theta) / theta * theta_new * (B_new - B)

      loss[iter] <- (sum((Y - X %*% B_new)^2) / 2
	 + sum(abs(tcrossprod(B_new, C)))
	 + lambda * sum(abs(B_new)))

      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < tol)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")

      B <- B_new
      theta <- theta_new
   }

   list(B=B, loss=loss[1:iter])
}

spgL2 <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, tol=1e-4, mu=1e-4)
{
   N <- nrow(X)
   J <- ncol(X)
   K <- ncol(Y)

   loss <- numeric(maxiter)
   B <- matrix(0, J, K)
   W <- B
   theta <- 1

   XX <- crossprod(X)
   XY <- crossprod(X, Y)
   CC <- crossprod(C)

   L <- specnorm(X) + gamma * specnorm(C)

   for(iter in 1:maxiter)
   {
      grad <- XX %*% W - XY + gamma * W %*% CC
      V <- W - (1 / L) * grad
      B_new <- soft_threshold(V, lambda / L)
      theta_new <- (sqrt(theta^4 + 4 * theta^2) - theta^2) / 2
      W <- B_new + (1 - theta) / theta * theta_new * (B_new - B)
      loss[iter] <- (sum((Y - X %*% B_new)^2) / 2
	 + gamma * sum(tcrossprod(B_new, C)^2)
	 + lambda * sum(abs(B_new)))

      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < tol)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")

      B <- B_new
      theta <- theta_new
   }

   list(B=B, loss=loss[1:iter])
}

