
graph.sqr <- function(R, threshold=NA) 
{
   s <- sign(R) * R^2
   diag(s) <- 0
   s
}

graph.ind <- function(R, threshold)
{
   s <- sign(R) * (abs(R) > threshold)
   diag(s) <- 0
   s
}

graph.abs <- function(R, threshold=NA)
{
   diag(R) <- 0
   R
}

fmpr <- function(..., type=c("cold", "warm"))
{
   type <- match.arg(type)
   if(type == "warm") {
      fmpr.warm(...)
   } else {
      fmpr.cold(...)
   }
}

# lambda1: scalar or K-vector
# lambda2: scalar or K-vector
# lambda3: scalar 
fmpr.cold <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, G=NULL,
      maxiter=1e5, eps=1e-6, verbose=FALSE, simplify=FALSE,
      sparse=TRUE, nzmax=nrow(X))
{
   if(length(X) == 0 || length(Y) == 0)
      stop("X and/or Y have zero length")

   p <- ncol(X)
   Y <- cbind(Y)
   K <- ncol(Y)

   if(nrow(X) != nrow(Y))
      stop("dimensions of X and Y don't agree")

   if(!is.null(G)) {
      if(nrow(G) != ncol(G) || nrow(G) != K) {
	 stop("dimensions of G and Y don't agree")
      }
   } else {
      G <- matrix(0, K, K)
      lambda3 <- 0
   }
	       
   B0 <- if(sparse) {
      sparseMatrix(i={}, j={}, dims=c(p, K))
   } else matrix(0, p, K)

   # fit models in increasing order of lambda1, without messing with
   # the original ordering of lambda1 requested by user
   l1ord <- order(lambda1, decreasing=FALSE)

   # We parallelise the lambda2 and lambda3 but not the lambda1. Assuming that
   # the lambda1 are in increasing order, we can stop if there are no active
   # variables for a given lambda2,lambda3, as incresing lambda1 will only
   # result in no active variables again. This allows us to not waste time on
   # models that will have zero active variables anyway.
   Btmp <- foreach(j=seq(along=lambda2)) %:%
	 foreach(k=seq(along=lambda3)) %dopar% {
	    Bjk <- vector("list", length(lambda1))
	    nactive <- numeric(length(lambda1))

	    for(i in seq(along=lambda1))
	    {
	       Bjk[[i]] <- B0
	    }
	       
	    for(i in seq(along=lambda1))
	    {
	       #if(verbose)
	          cat("\t", l1ord[i], j, k, ": ")

	       if(i > 1 && nactive[l1ord[i-1]] == 0)
	       {
		  cat("\nterminating early, no active variables left\n")
		  break
	       }
	       
	       r <- .C("fmpr_weighted", as.numeric(X), as.numeric(Y), 
	          numeric(p * K), nrow(X), ncol(X), K,
       	          as.numeric(rep(lambda1[l1ord[i]], K)),
	          as.numeric(rep(lambda2[j], K)),
	          as.numeric(rep(lambda3[k], K)),
       	          as.numeric(G), as.integer(maxiter),
       	          as.double(eps), as.integer(verbose), integer(1),
	          integer(1), integer(1)
	       )
	       status <- r[[14]]
	       nactive[l1ord[i]] <- r[[16]]
	       numiter <- r[[15]]
	       if(!status) {
	          cat("fmpr failed to converge within ",
	             maxiter, " iterations")
	       } else #if(verbose) {
	          {
		  cat("converged in", numiter, "iterations",
	             "with", nactive[l1ord[i]], "active variables\n\n")
	       }
	          
	       m <- matrix(r[[3]], p, K)

	       Bjk[[l1ord[i]]] <- if(sparse) {
	          w <- which(m != 0, arr.ind=TRUE)
	          sparseMatrix(i=w[,1], w[,2], x=m[w], dims=c(p, K))
	       } else m
	    }
	    Bjk
   }

   # Invert the list so that lambda1 is on the first dimension again
   B <- lapply(seq(along=lambda1), function(i) {
      lapply(seq(along=lambda2), function(j) {
	 lapply(seq(along=lambda3), function(k) {
	    Btmp[[j]][[k]][[i]]
	 })
      })
   })

   if(simplify && length(lambda1) == 1 
      && length(lambda2) == 1 
      && length(lambda3) == 1) {
      B[[1]][[1]][[1]]
   } else {
      B
   }
}

fmpr.warm <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, G=NULL,
      maxiter=1e5, eps=1e-6, verbose=FALSE, simplify=FALSE,
      sparse=FALSE, nzmax=nrow(X))
{
   if(length(X) == 0 || length(Y) == 0)
      stop("X and/or Y have zero length")

   p <- ncol(X)
   Y <- cbind(Y)
   K <- ncol(Y)
   N <- nrow(X)

   if(nrow(X) != nrow(Y))
      stop("dimensions of X and Y don't agree")

   if(!is.null(G)) {
      if(nrow(G) != ncol(G) || nrow(G) != K) {
	 stop("dimensions of G and Y don't agree")
      }
   } else {
      G <- matrix(0, K, K)
      lambda3 <- 0
   }
	       
   #B0 <- if(sparse) {
   #   sparseMatrix(i={}, j={}, dims=c(p, K))
   #} else matrix(0, p, K)
   B0 <- matrix(0, p, K)
   LP0 <- matrix(0, N, K)

   # fit models in increasing order of lambda1, without messing with
   # the original ordering of lambda1 requested by user
   l1ord <- order(lambda1, decreasing=TRUE)

   # We parallelise the lambda2 and lambda3 but not the lambda1. Assuming that
   # the lambda1 are in increasing order, we can stop if there are no active
   # variables for a given lambda2,lambda3, as incresing lambda1 will only
   # result in no active variables again. This allows us to not waste time on
   # models that will have zero active variables anyway.
   Btmp <- foreach(j=seq(along=lambda2)) %:%
	 foreach(k=seq(along=lambda3)) %dopar% {
	    Bjk <- vector("list", length(lambda1))
	    LPjk <- vector("list", length(lambda1))
	    nactive <- numeric(length(lambda1))

	    #for(i in seq(along=lambda1))
	    #{
	    #   Bjk[[i]] <- B0
	    #}
	       
	    for(i in seq(along=lambda1))
	    {
	       #if(verbose)
	          cat("\t", l1ord[i], j, k, ": ")

	       #if(i > 1 && nactive[l1ord[i-1]] == 0)
	       #{
	       #   cat("\nterminating early, no active variables left\n")
	       #   break
	       #}

	       if(i == 1) {
		  B <- B0
		  LP <- LP0
	       } else {
		  B <- Bjk[[l1ord[i-1]]]
		  LP <- LPjk[[l1ord[i-1]]]
	       }
	       
	       r <- .C("fmpr_weighted_warm",
		  as.numeric(X),
		  as.numeric(Y), 
	          as.numeric(B),
		  as.numeric(LP),
		  nrow(X),
		  ncol(X),
		  K,
       	          as.numeric(rep(lambda1[l1ord[i]], K)),
	          as.numeric(rep(lambda2[j], K)),
	          as.numeric(rep(lambda3[k], K)),
       	          as.numeric(G),
		  as.integer(maxiter),
       	          as.double(eps),
		  as.integer(verbose),
		  integer(1),
	          integer(1),
		  integer(1)
	       )
	       status <- r[[15]]
	       numiter <- r[[16]]
	       nactive[l1ord[i]] <- r[[17]]
	       if(!status) {
	          cat("fmpr failed to converge within ",
	             maxiter, " iterations")
	       } else #if(verbose) {
	          {
		  cat("converged in", numiter, "iterations",
	             "with", nactive[l1ord[i]], "active variables\n\n")
	       }
	          
	       #m <- matrix(r[[3]], p, K)
	       #
	       #Bjk[[l1ord[i]]] <- if(sparse) {
	       #   w <- which(m != 0, arr.ind=TRUE)
	       #   sparseMatrix(i=w[,1], w[,2], x=m[w], dims=c(p, K))
	       #} else m
	       Bjk[[l1ord[i]]] <- matrix(r[[3]], p, K)
	       LPjk[[l1ord[i]]] <- matrix(r[[4]], N, K)
	    }
	    Bjk
   }

   # Invert the list so that lambda1 is on the first dimension again
   B <- lapply(seq(along=lambda1), function(i) {
      lapply(seq(along=lambda2), function(j) {
	 lapply(seq(along=lambda3), function(k) {
	    Btmp[[j]][[k]][[i]]
	 })
      })
   })

   if(simplify && length(lambda1) == 1 
      && length(lambda2) == 1 
      && length(lambda3) == 1) {
      B[[1]][[1]][[1]]
   } else {
      B
   }
}

lasso <- function(X, y, lambda1=0,
      maxiter=1e5, eps=1e-6, verbose=FALSE)
{
   p <- ncol(X)

   sapply(lambda1, function(l) {
      r <- .C("lasso", as.numeric(X), as.numeric(y), numeric(p), 
         nrow(X), ncol(X), as.numeric(l),
         as.integer(maxiter), as.double(eps),
         as.integer(verbose)
      )
      matrix(r[[3]], p)
   })
}

# Can handle ncol(Y) > 1
ridge <- function(X, Y, lambda2=0)
{
   Y <- cbind(Y)
   XX <- crossprod(X)
   XY <- crossprod(X, Y)
   K <- ncol(Y)
   p <- ncol(X)
   l <- foreach(l=lambda2) %dopar% {
      b <- try(qr.solve(XX + diag(p) * l, XY), silent=TRUE)
      if(is(b, "try-error")) {
	 matrix(0, p, K)
      } else {
	 b
      }
   }
   l
}

# zero mean, unit norm (not unit variance)
standardise <- function(x)
{
   x1 <- sweep(x, 2, colMeans(x))
   v <- apply(x1, 2, function(z) sqrt(sum(z^2)))
   s <- sweep(x1, 2, v, FUN="/")
   s[is.na(s)] <- 0
   s
}

# Returns the maximal l1 penalty for each task, i.e., the smallest l1 penalty
# that makes all the weights zero.
maxlambda1 <- function(X, Y)
{
   Y <- cbind(Y)
   r <- .C("maxlambda1", as.numeric(X), as.numeric(Y),
      numeric(ncol(Y)), nrow(X), ncol(X), ncol(Y))
   r[[3]]
}

# Converts an N by p matrix into a block-diagonal N*K by p*K matrix
blockX <- function(X, p, K)
{
   N <- nrow(X)
   
   Xblock <- matrix(0, N * K, p * K) 
   s1 <- seq(1, N * K, N)
   s2 <- seq(1, p * K, p)

   for(i in 1:K)
   {
      j <- s1[i]:(s1[i] + N - 1)
      k <- s2[i]:(s2[i] + p - 1)
      Xblock[j, k] <- X
   }

   Xblock
}

# Makes the edges by vertices matrix for SPG
# Same as gennetwork.m
# 
# cortype: 0: thresholding to binary values
#          1: |R|
#          2: R^2
gennetwork <- function(Y, corthresh=0.5, cortype=1)
{
   R <- cor(Y)

   R[abs(R) < corthresh] <- 0

   K <- ncol(Y)
   nV <- K
   UR <- R
   UR[lower.tri(UR)] <- 0
   diag(UR) <- 0
   
   W <- if(cortype == 0) {
      (abs(UR) > corthresh) + 0
   } else if(cortype == 1) {
      abs(UR)
   } else if(cortype == 2) {
      UR^2
   } else {
      stop("unknown cortype:", cortype)
   }

   nzUR <- which(UR != 0)
   E <- which(UR != 0, arr.ind=TRUE)
   Ecoef <- W[nzUR]
   Esign <- sign(R[nzUR])
   nE <- nrow(E)
   C_I <- c(1:nE, 1:nE)
   C_J <- as.numeric(E)
   C_S <- cbind(Ecoef, -Ecoef * Esign)
   C <- matrix(0, nE, nV)
   C[cbind(C_I, C_J)] <- C_S
   C
}

spg <- function(X, Y, C=NULL, lambda=0, gamma=0, tol=1e-6,
      mu=1e-4, maxiter=1e4, simplify=FALSE, verbose=FALSE)
{
   K <- ncol(Y)
   N <- nrow(Y)
   p <- ncol(X)

   # SPG works on this scale
   lambda <- lambda * N

   if(length(C) == 0)
      C <- matrix(0, 1, K)

   XX <- crossprod(X)
   XY <- crossprod(X, Y)
   CNorm <- 2 * max(colSums(C^2))
   C0 <- C
   L0 <- eigen(XX, only.values=TRUE, symmetric=TRUE)$values[1]

   cat("CNorm:", CNorm, "\n")

   B <- foreach(i=seq(along=lambda)) %:% 
      foreach(j=seq(along=gamma)) %dopar% {
	 L <- L0 + gamma[j]^2 * CNorm / mu
	 C <- gamma[j] * C0

	 if(verbose)
	    cat("spg gamma:", gamma[j], "lambda:", lambda[i], "\n")

	 r <- .C("spg_core",
   	    as.numeric(XX), as.numeric(XY),
   	    as.numeric(X), as.numeric(Y),
   	    as.integer(N), as.integer(p), as.integer(K),
   	    numeric(p * K),
   	    as.numeric(C), as.integer(nrow(C)), as.numeric(L), 
   	    as.numeric(gamma[j]), as.numeric(lambda[i]),
   	    as.numeric(tol), as.numeric(mu), as.integer(maxiter),
   	    as.integer(verbose), integer(1), integer(1)
   	 )
   	 niter <- r[[18]]
	 status <- r[[19]]
	 if(!status) {
	    cat("SPG failed to converge within", niter, "iterations\n")
	 }

	 matrix(r[[8]], p, K)
   }
   
   if(simplify
      && length(lambda) == 1 
      && length(gamma) == 1) {
      B[[1]][[1]]
   } else {
      B
   }
}

fmprR <- function(X, Y, G, lambda1, lambda2, lambda3, maxiter=1e3, eps=1e-6)
{
   K <- ncol(Y)
   p <- ncol(X)
   N <- nrow(X)
   B <- matrix(0, p, K)

   v <- diag(crossprod(X))

   loss <- numeric(maxiter)
   for(iter in 1:maxiter)
   {
      for(k in 1:K)
      {
	 for(j in 1:p)
      	 {
      	    Err <- X %*% B[,k] - Y[, k]
      	    d1 <- crossprod(X[,j], Err)
      	    d2 <- v[j]

	    d1 <- d1 + lambda3 * sum(abs(G[k, -k]) 
		  * (B[j, k] - sign(G[k, -k]) * B[j, -k]))
	    d2 <- d2 + lambda3 * sum(abs(G[k, -k]))

	    s <- (B[j, k] - d1 / d2) / (1 + lambda2)
	    if(abs(s) <= lambda1)
	    {
	       B[j, k] <- 0
	    }
	    else
	    {
	       B[j, k] <- s - lambda1 * sign(s)
	    }
      	 }
      }
      loss[iter] <- mean((X %*% B - Y)^2)
      cat("iter", iter, "loss",
	 format(loss[iter], digits=-log10(eps)), "\n")
      if(iter > 1 && abs(loss[iter] - loss[iter-1]) <= eps)
	 break
   }
   B
}

