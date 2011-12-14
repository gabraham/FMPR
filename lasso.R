# plain l1 lasso for single outputs, very naive implementation, no active set
# convergence
#
# Assumes X is standardised and y is zero-mean
lasso <- function(X, y, lambda1=0, maxepoch=5000, eps=1e-3)
{
   p <- ncol(X)
   beta <- numeric(p)
   oldloss <- Inf
   loss <- Inf
   lp <- X %*% beta
   Err <- lp - y
   loss <- mean(Err^2)

   V <- apply(X, 2, function(x) sum(x^2))

   for(epoch in 1:maxepoch)
   {
      cat("epoch", epoch, "loss:", loss, "\n")   
      for(j in 1:p)
      {
	 grad <- crossprod(X[,j], Err)
	 if(grad != 0) {
	    beta_new <- soft(beta[j] - grad / V[j], lambda1)
	    diff <- beta_new - beta[j]
	    lp <- lp + x[,j] * diff
	    Err <- lp - y
	    beta[j] <- beta_new
	 }
      }
      loss <- mean(Err^2)
   }
   beta
}

soft <- function(b, g)
{
   sign(b) * pmax(abs(b) - g, 0)
}

lasso.activeset <- function(x, y, lambda1=0,
   maxepoch=5000, eps=1e-4)
{
   p <- ncol(x)
   n <- nrow(x)
   beta <- numeric(p)
   lp <- numeric(n)
   beta_old <- beta
   active_new <-  rep(TRUE, p)
   active_old <- rep(FALSE, p)
   allconverged <- 1
   loss <- Inf
   loss_old <- Inf
   converged <- logical(p)

   l <- coef(lm(y ~ x - 1))
   V <- apply(x, 2, function(x) sum(x^2))
   loss_null <- mean((y - mean(y))^2)
   cat("null loss:", loss_null, "\n")
   
   for(epoch in 1:maxepoch)
   {
      cat("epoch:", epoch, "loss:", loss, "\n")

      for(j in 1:p)
      {
	 loss_old <- loss
	 converged[j] <- TRUE
	 if(active_new[j]) 
	 {
	    converged[j] <- FALSE
	    s <- crossprod(x[,j], lp - y) / V[j]
	    beta_new <- soft(beta[j] - s, lambda1)
	    diff <- beta_new - beta[j]
	    lp <- lp + x[,j] * diff
	    beta[j] <- beta_new
	    active_new[j] <- beta[j] != 0
	    loss <- mean((lp - y)^2)
	    converged[j] <- abs(loss - loss_old) < 1e-5 * loss_null
	 }
      }

      if(all(converged))
      {
	 if(allconverged == 1){
	    active_old <- active_new
	    active_new[] <- TRUE
	    allconverged <- 2
	 } else {
	    if(all(active_old == active_new)) {
	       cat("terminating at epoch", epoch, "\n")
	       break
	    } else {
	       allconverged <- 1
	       active_old <- active_new
	       active_new[] <- TRUE
	    }
	 }
      }

      beta_old <- beta
      loss_old <- loss
   }

   cat("\n")
   beta
}

sphere <- function(X, y, lambda1max, lambda)
{
   X <- scale(X) / sqrt(nrow(X)-1)
   z <- drop(abs(crossprod(X, y))) 
   z < lambda1max * (1 - 2 * sqrt(1/lambda1max^2 - 1) * (lambda1max/lambda - 1))
}

# Xiang NIPS 2011 screening
lasso.sphere <- function(X, y, lambda1=0, maxepoch=5000, eps=1e-3)
{
   X <- scale(X)
   p <- ncol(X)
   beta <- numeric(p)
   oldloss <- Inf
   loss <- Inf
   lp <- X %*% beta
   Err <- lp - y
   loss <- mean(Err^2)

   s <- sphere(X, y, lambda1)

   V <- apply(X, 2, function(x) sum(x^2))

   for(epoch in 1:maxepoch)
   {
      cat("epoch", epoch, "loss:", loss, "\n")   
      for(j in 1:p)
      {
	 grad <- crossprod(X[,j], Err)
	 if(grad != 0) {
	    beta_new <- soft(beta[j] - grad / V[j], lambda1)
	    diff <- beta_new - beta[j]
	    lp <- lp + x[,j] * diff
	    Err <- lp - y
	    beta[j] <- beta_new
	 }
      }
      loss <- mean(Err^2)
   }
   beta
}

