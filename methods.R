
# Wrappers for calling lasso and groupridge

require(foreach)

dyn.load("~/Code/L1L2/groupridge.so")

require(Matrix)

# lambda1: scalar or K-vector
# lambda2: scalar or K-vector
# lambda3: scalar 
groupridge <- function(X, Y, lambda1=0, lambda2=0, lambda3=0, G=NULL,
      maxiter=1e5, eps=1e-8, type="threshold", verbose=FALSE, simplify=FALSE,
      sparse=TRUE)
{
   p <- ncol(X)
   Y <- cbind(Y)
   K <- ncol(Y)

   #if(length(lambda1) == 1)
   #   lambda1 <- rep(lambda1, K)
   
   #if(length(lambda2) == 1)
   #   lambda2 <- rep(lambda2, K)

   if(is.null(G))
      G <- matrix(0, K, K)

   if(type == "threshold") {
      fn <- "groupridge_threshold"
      g <- as.integer(G)
   } else {
      fn <- "groupridge_weighted"
      g <- as.numeric(G)
   }

   B0 <- if(sparse) {
      sparseMatrix(i={}, j={}, dims=c(p, K))
   } else {
      matrix(0, p, K)
   }
   
   res <- lapply(seq(along=lambda1), function(i) {
      lapply(seq(along=lambda2), function(j) {
	 lapply(seq(along=lambda3), function(k) {
	    B0 
	 })
      })
   })

   #nz <- matrix(0, 

   # Assumes that lambda2, lambda3 are sorted in
   # increasing order, and that lambda1 is in decreasing order
   #for(i in seq(along=lambda1))
   foreach(i=seq(along=lambda1)) %dopar% {
      #for(j in seq(along=lambda2))
      foreach(j=seq(along=lambda2)) %dopar% {
	 # lambda2 is in increasing order. If a given lambda1 value
	 # with the previous lambda2 produces a model with all zeros,
	 # then the next lambda2 will also be all zero since
	 # increasing lambda2 only reduces the gradient.

	 #for(k in seq(along=lambda3))
	 foreach(k=seq(along=lambda3)) %dopar% {
	    # groupridge expects l1/l2/l3 to be a vector of length K,
	    # allowing for a different penalty for each task, but we
	    # don't use this feature here
	    r <- .C(fn, as.numeric(X), as.numeric(Y), 
	       numeric(p * K), nrow(X), ncol(X), K,
       	       as.numeric(rep(lambda1[i], K)),
	       as.numeric(rep(lambda2[j], K)),
	       as.numeric(rep(lambda3[k], K)),
       	       g, as.integer(maxiter),
       	       as.double(eps), as.integer(verbose), integer(1),
	       integer(1)
	    )
	    status <- r[[14]]
	    if(!status) {
	       warning("groupridge failed to converge within ",
	          maxiter, " iterations")
	    } else if(verbose) {
	       cat("converged in", r[[15]], "iterations\n\n")
	    }
	       
	    m <- matrix(r[[3]], p, K)

	    res[[i]][[j]][[k]] <- if(sparse) {
	       w <- which(m != 0, arr.ind=TRUE)
	       sparseMatrix(i=w[,1], w[,2], x=m[w], dims=c(p, K))
	    } else m
	 }
      }
   }

   if(simplify && length(lambda1) == 1 
      && length(lambda2) == 1 
      && length(lambda3) == 1) {
      res[[1]][[1]][[1]]
   } else {
      res
   }
}

# warm restarts
groupridge.warm <- function(X, Y, B=NULL, lambda1=0, lambda2=0, lambda3=0, G=NULL,
      maxiter=1e5, eps=1e-8, verbose=FALSE)
{
   p <- ncol(X)
   Y <- cbind(Y)
   K <- ncol(Y)

   #if(length(lambda1) == 1)
   #   lambda1 <- rep(lambda1, K)
   
   #if(length(lambda2) == 1)
   #   lambda2 <- rep(lambda2, K)

   if(is.null(G))
      G <- matrix(0, K, K)

   #if(is.null(B)) {
   #   LP <- numeric(N * K)
   #} else {
   #   LP <- as.numeric(X %*% B)
   #}

   len1 <- length(lambda1)
   len2 <- length(lambda2)
   len3 <- length(lambda3)

   B <- array(0, c(len1, len2, len3, p, K))
   LP <- array(0, c(len1, len2, len3, N, K))
   
   for(i in 1:len1)
   {
      L <- if(i == 1) {
	 numeric(N * K) 
      } else {
	 LP[i - 1, j, k , ,]
      }
      for(j in 1:len2)
      {
	 for(k in 1:len3)
	 {
	    r <- .C("groupridge5", as.numeric(X), as.numeric(Y), 
	       as.numeric(B[i, j, k, , ]), as.numeric(L),
	       nrow(X), ncol(X), K,
	       as.numeric(rep(lambda1[i], K)), as.numeric(rep(lambda2[j], K)),
	       as.numeric(rep(lambda3[k], K)),
	       as.integer(G), as.integer(maxiter),
	       as.double(eps), as.integer(verbose), integer(1),
	       integer(1))
	    status <- r[[15]]
	    if(!status) {
	       warning("groupridge failed to converge within ",
		  maxiter, " iterations")
	    } else if(verbose) {
	       cat("converged in", r[[16]], "iterations\n")
	    }
	    B[i, j, k, , ] <- matrix(r[[3]], p, K)
	    LP[i, j, k, , ] <- matrix(r[[4]], N, K)
	 }
      }
   }
   B
}

lasso <- function(X, y, lambda1=0,
      maxiter=1e5, eps=1e-8, verbose=FALSE)
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

ridge <- function(X, Y, lambda2=0)
{
   Y <- cbind(Y)
   XX <- crossprod(X)
   XY <- crossprod(X, Y)
   K <- ncol(Y)
   p <- ncol(X)
   sapply(lambda2, function(l) {
      b <- try(qr.solve(XX + diag(p) * l, XY), silent=TRUE)
      if(is(b, "try-error")) {
	 matrix(0, p, K)
      } else {
	 b
      }
   })
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

## Threshold the correlation matrix and convert into groups, finding the graph
## that spans the connected vertices. Each connected subgraph is a group.
#cor2grp <- function(Y, thresh=0.7)
#{
#   R <- abs(cor(Y))
#   diag(R) <- 0
#   w <- which(R >= thresh, arr.ind=TRUE)
#   K <- numeric(length(unique(w)))
#}

