# This function is a bit confusing: 
# it computes the squared spectral norm of X^T X, using SVD of X.
#0
# The squared spectral norm is the same as the largest eigenvalue.
#
# maxeigen(X) is equivalent to R's eigen(crossprod(X))$values[1] and
# to MATLAB's eigs(X'X, 1)
maxeigen <- function(X)
{
   safe.svd(X, nu=0, nv=1)$d[1]^2
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
   L <- maxeigen(X) + gamma^2 * CNorm / mu

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

spgL1cd <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, tol=1e-4, mu=1e-4)
{
   N <- nrow(X)
   J <- ncol(X)
   K <- ncol(Y)

   loss <- numeric(maxiter)
   B <- matrix(0, J, K)
   Err <- matrix(0, N, p)

   CNorm <- 2 * max(colSums(C^2))
   C <- gamma * C
   L <- maxeigen(X) + gamma^2 * CNorm / mu

   for(iter in 1:maxiter)
   {
      for(k in 1:K)
      {
	 for(j in 1:J)
      	 {
	    Err[,k] <- X %*% B[,k] - Y[,k]
	    A <- hard_threshold(tcrossprod(C, B / mu), 1)
	    grad <- crossprod(X[,j], Err[,k]) + crossprod(A, C)[j, k]
	    B[j, k] <- soft_threshold(B[j, k] - grad / L, lambda / L)
      	 }
      }

      loss[iter] <- (sum((Y - X %*% B)^2) * 0.5
	 + sum(abs(tcrossprod(B, C)))
	 + lambda * sum(abs(B)))

      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < tol)
	 break
      cat("iter", iter, "loss:", loss[iter], "\n")
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

   L <- maxeigen(X) + gamma * maxeigen(C)

   for(iter in 1:maxiter)
   {
      grad <- XX %*% W - XY + gamma * W %*% CC
      V <- W - (1 / L) * grad
      B_new <- soft_threshold(V, lambda / L)
      theta_new <- (sqrt(theta^4 + 4 * theta^2) - theta^2) / 2
      W <- B_new + (1 - theta) / theta * theta_new * (B_new - B)
      loss[iter] <- (0.5 * sum((Y - X %*% B_new)^2)
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

# We scale the loss by 1/N, and subsequently the Lipschitz constant
spgL1v2 <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, tol=1e-4, mu=1e-4)
{
   N <- nrow(X)
   J <- ncol(X)
   K <- ncol(Y)

   loss <- numeric(maxiter)
   B <- matrix(0, J, K)
   W <- B
   theta <- 1

   XX <- crossprod(X) / N
   XY <- crossprod(X, Y) / N
   CNorm <- 2 * max(colSums(C^2))

   C <- gamma * C
   L <- maxeigen(X) / N + gamma^2 * CNorm / mu

   for(iter in 1:maxiter)
   {
      A <- hard_threshold(tcrossprod(C, W / mu), 1)
      grad <- XX %*% W - XY + crossprod(A, C)
      V <- W - (1 / L) * grad
      B_new <- soft_threshold(V, lambda / L)
      theta_new <- (sqrt(theta^4 + 4 * theta^2) - theta^2) / 2
      W <- B_new + (1 - theta) / theta * theta_new * (B_new - B)

      loss[iter] <- (0.5 * sum((Y - X %*% B_new)^2) / N
	 + sum(abs(tcrossprod(B_new, C))) # C already includes gamma
	 + lambda * sum(abs(B_new)))

      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < tol)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")

      B <- B_new
      theta <- theta_new
   }

   list(B=B, loss=loss[1:iter])
}

# We scale the loss by 1/N, and subsequently the Lipschitz constant
spgL2v2 <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, tol=1e-4, mu=1e-4)
{
   N <- nrow(X)
   J <- ncol(X)
   K <- ncol(Y)

   loss <- numeric(maxiter)
   B <- matrix(0, J, K)
   W <- B
   theta <- 1

   XX <- crossprod(X) / N
   XY <- crossprod(X, Y) / N
   CC <- crossprod(C)

   L <- maxeigen(X) / N + gamma * maxeigen(C)

   for(iter in 1:maxiter)
   {
      grad <- XX %*% W - XY + gamma * W %*% CC
      V <- W - (1 / L) * grad
      B_new <- soft_threshold(V, lambda / L)
      theta_new <- (sqrt(theta^4 + 4 * theta^2) - theta^2) / 2
      W <- B_new + (1 - theta) / theta * theta_new * (B_new - B)

      loss[iter] <- (0.5 * sum((Y - X %*% B_new)^2) / N
	 + 0.5 * gamma * sum(tcrossprod(B_new, C)^2)
	 + lambda * sum(abs(B_new)))

      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < tol)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")

      B <- B_new
      theta <- theta_new
   }

   list(B=B, loss=loss[1:iter])
}

fmprR <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, eps=1e-4)
{
   p <- ncol(X)
   K <- ncol(Y)
   N <- nrow(X)
   B <- matrix(0, p, K)

   loss <- numeric(maxiter)
   Err <- matrix(0, N, K)

   CC <- diag(crossprod(C))
   XX <- diag(crossprod(X) / N)

   for(iter in 1:maxiter)
   {
      for(k in 1:K)
      {
	 for(j in 1:p)
	 {
	    Err[,k] <- X %*% B[,k] - Y[,k]
	    d1 <- crossprod(X[, j], Err[,k]) / N + gamma * B[j, k] * CC[k]
	    #d2 <- XX[j] + gamma * CC[k]
	    d2 <- 1e4
	    B[j, k] <- soft_threshold(B[j, k] - d1 / d2, lambda)
      
	    loss[iter] <- (0.5 * sum((Y - X %*% B)^2) / N
	       + 0.5 * gamma * sum(tcrossprod(B, C)^2)
	       + lambda * sum(abs(B)))
	 }
      }
    
      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < eps)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")
   }

   B
}

fmprR2 <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, eps=1e-4)
{
   p <- ncol(X)
   K <- ncol(Y)
   N <- nrow(X)
   B <- matrix(0, p, K)

   loss <- numeric(maxiter)
   Err <- matrix(0, N, K)
   
   CC <- crossprod(C)
   XX <- diag(crossprod(X) / N)

   for(iter in 1:maxiter)
   {
      for(k in 1:K)
      {
	 for(j in 1:p)
	 {
	    Err[,k] <- X %*% B[,k] - Y[,k]
	    #d1 <- crossprod(X[, j], Err[,k]) / N + gamma * B[j, k] * CC[k]
	    d1 <- crossprod(X[, j], Err[,k]) / N + gamma * (B %*% CC)[j, k]
	    #d2 <- L
	    W <- B
	    W[] <- 0
	    W[j, k] <- 1
	    d2 <- XX[j] + gamma * (W %*% CC)[j, k]
	    B[j, k] <- soft_threshold(B[j, k] - d1 / d2, lambda)
      
	    loss[iter] <- (0.5 * sum((Y - X %*% B)^2) / N
	       + 0.5 * gamma * sum(tcrossprod(B, C)^2)
	       + lambda * sum(abs(B)))
	 }
      }
    
      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < eps)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")
   }

   B
}

# Lipschitz constant for step size
fmprR3 <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, eps=1e-4)
{
   p <- ncol(X)
   K <- ncol(Y)
   N <- nrow(X)
   B <- matrix(0, p, K)

   loss <- numeric(maxiter)
   Err <- matrix(0, N, K)
   
   L <- maxeigen(X) / N + gamma * maxeigen(C)

   CC <- crossprod(C)
   XX <- diag(crossprod(X) / N)

   for(iter in 1:maxiter)
   {
      for(k in 1:K)
      {
	 for(j in 1:p)
	 {
	    Err[,k] <- X %*% B[,k] - Y[,k]
	    d1 <- crossprod(X[, j], Err[,k]) / N + gamma * (B %*% CC)[j, k]
	    B[j, k] <- soft_threshold(B[j, k] - d1 / L, lambda / L)
      
	    loss[iter] <- (0.5 * sum((Y - X %*% B)^2) / N
	       + 0.5 * gamma * sum(tcrossprod(B, C)^2)
	       + lambda * sum(abs(B)))
	 }
      }
    
      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < eps)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")
   }

   B
}

fmprR5 <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, eps=1e-4)
{
   p <- ncol(X)
   K <- ncol(Y)
   N <- nrow(X)
   B <- matrix(0, p, K)

   loss <- numeric(maxiter)
   Err <- matrix(0, N, K)
   
   CC <- crossprod(C)
   XX <- diag(crossprod(X) / N)

   for(iter in 1:maxiter)
   {
      for(k in 1:K)
      {
	 for(j in 1:p)
	 {
	    Err[,k] <- X %*% B[,k] - Y[,k]
	    #d1 <- crossprod(X[, j], Err[,k]) / N + gamma * B[j, k] * CC[k]
	    d1 <- crossprod(X[, j], Err[,k]) / N + gamma * (B %*% CC)[j, k]
	    #d2 <- L
	    W <- B
	    W[] <- 0
	    W[j, k] <- 1
	    d2 <- XX[j] + gamma * (W %*% CC)[j, k]
	    B[j, k] <- soft_threshold(B[j, k] - d1 / d2, lambda)
      
	    loss[iter] <- (0.5 * sum((Y - X %*% B)^2) / N
	       + 0.5 * gamma * sum(tcrossprod(B, C)^2)
	       + lambda * sum(abs(B)))
	 }
      }
    
      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < eps)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")
   }

   B
}

fmprR6 <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e3, eps=1e-4)
{
   p <- ncol(X)
   K <- ncol(Y)
   N <- nrow(X)
   B <- matrix(0, p, K)

   loss <- numeric(maxiter)
   Err <- matrix(0, N, K)
   
   CC <- crossprod(C)
   XX <- diag(crossprod(X) / N)
   L <- maxeigen(X) / N + gamma * maxeigen(C)

   for(iter in 1:maxiter)
   {
      for(k in 1:K)
      {
	 for(j in 1:p)
	 {
	    Err[,k] <- X %*% B[,k] - Y[,k]
	    d1 <- crossprod(X[, j], Err[,k]) / N 
	    #g <- sapply(1:nrow(C), function(e) {
	    #   B[j, ] %*% C[e, ] * C[e, k]
	    #})
	    #d1 <- d1 + gamma * sum(g)
	    d1 <- d1 + gamma * (B %*% CC)[j, k]
	    #d2 <- XX[j] + gamma * sum(C[,k]^2)
	    #B[j, k] <- soft_threshold(B[j, k] - d1 / d2, lambda)
	    B[j, k] <- soft_threshold(B[j, k] - d1 / L, lambda / L)
      
	    loss[iter] <- (0.5 * sum((Y - X %*% B)^2) / N
	       + 0.5 * gamma * sum(tcrossprod(B, C)^2)
	       + lambda * sum(abs(B)))
	 }
      }
    
      if(iter > 1 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < eps)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")
   }

   B
}

gd <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e6, eps=1e-4, mu=1e-4)
{
   p <- ncol(X)
   K <- ncol(Y)
   N <- nrow(X)
   B <- matrix(0, p, K)

   loss <- numeric(maxiter)
   Err <- matrix(0, N, K)

   CC <- crossprod(C)
   XX <- diag(crossprod(X) / N)

   for(iter in 1:maxiter)
   {
      Err <- X %*% B - Y
      d1 <- crossprod(X, Err) / N + gamma * B %*% CC
      B <- soft_threshold(B - d1 * mu, lambda)
      
      loss[iter] <- (0.5 * sum((Y - X %*% B)^2) / N
         + 0.5 * gamma * sum(tcrossprod(B, C)^2)
         + lambda * sum(abs(B)))
    
      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < eps)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")
   }

   B
}

gd2 <- function(X, Y, C, lambda=0, gamma=0, maxiter=1e6, eps=1e-4)
{
   p <- ncol(X)
   K <- ncol(Y)
   N <- nrow(X)
   B <- matrix(0, p, K)

   loss <- numeric(maxiter)
   Err <- matrix(0, N, K)

   CC <- crossprod(C)
   L <- maxeigen(X) / N + gamma * maxeigen(C)

   for(iter in 1:maxiter)
   {
      Err <- X %*% B - Y
      grad <- crossprod(X, Err) / N + gamma * B %*% CC
      B <- soft_threshold(B - grad / L, lambda / L)
      
      loss[iter] <- (0.5 * sum((Y - X %*% B)^2) / N
         + 0.5 * gamma * sum(tcrossprod(B, C)^2)
         + lambda * sum(abs(B)))
    
      if(iter > 10 && abs(loss[iter] - loss[iter - 1]) / abs(loss[iter - 1]) < eps)
	 break

      cat("iter", iter, "loss:", loss[iter], "\n")
   }

   B
}

