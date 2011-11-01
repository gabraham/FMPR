
# Formulation of MPGL as group lasso

p <- 10
C <- 5
N <- rep(100, C)
pop <- rep(1:C, N)
Bpop <- rep(1:C, each=p)

lambda1 <- 0
lambda2 <- 0

B <- cbind(unlist(lapply(1:C, function(i) {
   rnorm(p, i, 0.01)
})))

X <- matrix(0, sum(N), p * C) 
s1 <- seq(1, sum(N), N[1])
s2 <- seq(1, p * C , p)

for(i in 1:C)
{
   X[s1[i]:(s1[i] + N[i] - 1), s2[i]:(s2[i] + p - 1)] <- scale(rnorm(N[i] * p))
}

Y <- scale(X %*% B, scale=FALSE)

Bhat <- cbind(numeric(p * C))
s <- seq(0, p * C - 1, p)

oldloss <- 0
eps <- 1e-4
maxiter <- 10000

for(iter in 1:maxiter)
{
   Err <- X %*% Bhat - Y
   loss <- mean(Err^2)
   cat(iter, loss, "\n")
   grad <- crossprod(X, Err) / (N - 1) / 200 
   for(j in 1:p)
   {
      S <- Bhat[s + j,] - grad[s + j, ]
      normS <- sqrt(drop(crossprod(S)))
      #if(normS > 0) {
      #	 Bhat[s + g, ] <- pmax(0, 1 - lambda1 / normS) * S / (1 + lambda2)
      #}
      Bhat[s + j, ] <- S
   }

   if(abs(loss - oldloss) <= eps) {
      break
   }
   oldloss <- loss
   
}

