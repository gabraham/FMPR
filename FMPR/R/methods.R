
lasso <- function(X, y, lambda=0,
      maxiter=1e5, eps=1e-4, verbose=FALSE)
{
   p <- ncol(X)

   sapply(lambda, function(l) {
      r <- .C("lasso", as.numeric(X), as.numeric(y), numeric(p), 
         nrow(X), ncol(X), as.numeric(l),
         as.integer(maxiter), as.double(eps),
         as.integer(verbose)
      )
      matrix(r[[3]], p)
   })
}

## Can handle ncol(Y) > 1
#ridge <- function(X, Y, lambda=0)
#{
#   Y <- cbind(Y)
#   XX <- crossprod(X)
#   XY <- crossprod(X, Y)
#   res <- lapply(lambda, function(l){
#      b <- NULL
#      while(is.null(b) || is(b, "try-error")) {
#	 b <- try(qr.solve(XX + diag(p) * l, XY), silent=FALSE)
#	 l <- l * 1.1
#      }
#      b
#   })
#   res
#}

ridge <- function(X, Y, lambda=1e-3)
{
   Y <- cbind(Y)
   s <- try(svd(X))

   # sometimes svd works for X^T but not for X
   if(is(s, "try-error")) {
      s <- try(svd(t(X)))
      if(is(s, "try-error")) {
	 stop(s)
      }
      # transposed X
      U <- s$v
      V <- s$u
   } else {
      U <- s$u
      V <- s$v
   }

   D <- diag(s$d)
   R <- U %*% D

   RR <- crossprod(R)
   RY <- crossprod(R, Y)
   res <- lapply(lambda, function(l){
      b <- NULL
      while(is.null(b) || is(b, "try-error")) {
	 b <- try(
	    V %*% qr.solve(RR + diag(ncol(R)) * l, RY), silent=FALSE
	 )
	 l <- l * 1.1
      }
      b
   })
   res
}

lm.simple <- function(X, Y)
{
   res <- lapply(1:ncol(Y), function(i) {
      lapply(1:ncol(X), function(j) {
	 l <- coef(summary(lm(Y[,i] ~ X[,j])))
	 if(nrow(l) == 2) {
	    c(beta=l[2, 1], pval=l[2, 4])
	 } else {
	    c(beta=NA, pval=NA)
	 }
      })
   })

   B <- sapply(res, sapply, function(x) x[1])
   pval <- sapply(res, sapply, function(x) x[2])

   list(B=B, pval=pval)
}

# CCA of all phenotypes on one predictor at a time, like PLINK.multivariate
cca <- function(X, Y, min.p=1e-16)
{
   res <- apply(X, 2, function(x) {
      if(var(x) <= 1e-9) {
	 return(NA)
      }
      r <- cancor(x, Y)
      capture.output({
	 p <- p.asym(r$cor, N=length(x), p=1, q=ncol(Y), tstat="Wilks")$p.value
      })
      max(p, min.p)
   })

   list(B=NULL, pval=res)
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
blockX <- function(X, K)
{
   N <- nrow(X)
   p <- ncol(X)
   
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
# cortype: 1: |R|
#          2: R^2
#
# full: logical. if TRUE, the full (K-1)K by K  matrix will be created, even
# for zero-weight edges. Otherwise, only the edges with non-zero weight will
# be selected.
gennetwork <- function(Y, cortype=1, cormethod=c("pearson", "partial"))
{
   cormethod <- match.arg(cormethod)
   
   R <- if(cormethod == "pearson") {
      cor(Y)
   } else if(cormethod == "partial") {
      pcor.shrink(Y)
   } else {
      stop("unknown cormethod '", type, "'")
   }

   K <- ncol(Y)
   nV <- K
   UR <- R
   UR[lower.tri(UR)] <- 0
   diag(UR) <- 0
   
   W <- if(cortype == 1) {
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

# svd can fail for perfectly correlated data, so use irlba if it fails.
safe.svd <- function(x, ...)
{
   s <- try(svd(x, ...), silent=FALSE)
   if("try-error" %in% class(s)) {
      cat("using irlba instead\n")
      irlba(x, ...)
   } else {
      s
   }
}

spg <- function(X, Y, C=NULL, lambda=0, gamma=0, tol=1e-4,
   mu=1e-4, maxiter=1e4, simplify=FALSE, verbose=FALSE, type=c("l1", "l2"),
   divbyN=FALSE)
{
   K <- ncol(Y)
   N <- nrow(Y)
   p <- ncol(X)

   if(length(C) == 0)
      C <- matrix(0, 1, K)

   # Don't normalise by N here, we do in the C code
   XX <- crossprod(X)
   XY <- crossprod(X, Y)
   CNorm <- 2 * max(colSums(C^2))
   C0 <- C

   type <- match.arg(type)
   fun <- switch(type, "l1"="spg_core", "l2"="spg_l2_core")

   if(verbose) {
      cat("\nSPG: type=", type, ", lambda=", lambda, "gamma=", gamma,
	 "CNorm=", CNorm, "\n")
   }

   # Lipschitz constant of squared loss
   #
   # The largest eigenvalue of X^T X is upper-bounded by the Frobenius norm
   L0 <- if(p < 1e4L) {
      maxeigen(X)
   } else {
      sum(XX^2)
   }

   if(divbyN) {
      L0 <- L0 / N
   }

   # Lipschitz constant of l2 fusion penalty
   eigCC <- if(type == "l1") {
      NULL
   } else {
      maxeigen(C)
   }

   B <- lapply(seq(along=lambda), function(i) {
      lapply(seq(along=gamma), function(j) {
	 
	 # Compute final Lipschitz constant.
	 # SPG with l1 fusion folds the penalty gamma into the matrix C. The
	 # l2 fusion doesn't, because
	 # \gamma ||B C^T C||_F^2 != ||B Cg^T Cg||_F^2 where Cg = \gamma * C
	 if(type == "l1") {
	    L <- L0 + gamma[j]^2 * CNorm / mu
	    Cj <- gamma[j] * C0
	 } else if(type == "l2") {
	    L <- L0 + gamma[j] * eigCC
	    Cj <- C0
	 }

	 if(verbose)
	    cat("spg L:", L, "gamma:", gamma[j], "lambda:", lambda[i], "\n")

	 r <- .C(fun,
   	    as.numeric(XX),        # 1
	    as.numeric(XY),        # 2
   	    as.numeric(X),         # 3
	    as.numeric(Y),         # 4
   	    as.integer(N),         # 5
	    as.integer(p),         # 6
	    as.integer(K),         # 7
   	    numeric(p * K),        # 8   B
   	    as.numeric(Cj),        # 9   C
	    as.integer(nrow(Cj)),  # 10  nE
	    as.numeric(L),         # 11  L
   	    as.numeric(gamma[j]),  # 12  gamma
	    as.numeric(lambda[i]), # 13  lambda
   	    as.numeric(tol),       # 14  tol
	    as.numeric(mu),        # 15  mu
	    as.integer(maxiter),   # 16  maxiter
   	    as.integer(verbose),   # 17  verbose
	    integer(1),            # 18  niter
	    integer(1),            # 19  status
	    as.integer(divbyN)     # 20  divbyN
   	 )
   	 niter <- r[[18]]
	 status <- r[[19]]
	 if(!status) {
	    warning("SPG failed to converge within ", niter, " iterations,",
	       " lambda=", lambda[i], " gamma=", gamma[j], "\n")
	 }

	 if(verbose)
	    cat("\n")

	 matrix(r[[8]], p, K)
      })
   })
   
   if(simplify
      && length(lambda) == 1 
      && length(gamma) == 1) {
      B[[1]][[1]]
   } else {
      B
   }
}

fmpr <- function(X, Y, lambda=0, lambda2=0, gamma=0, C=NULL,
      maxiter=1e5, eps=1e-6, verbose=FALSE, simplify=FALSE,
      sparse=FALSE, nzmax=NULL, warm=TRUE, divbyN=TRUE, dup=TRUE)
{
   if(length(X) == 0 || length(Y) == 0)
      stop("X and/or Y have zero length")

   p <- ncol(X)
   Y <- cbind(Y)
   K <- ncol(Y)
   N <- nrow(X)

   if(nrow(X) != nrow(Y))
      stop("dimensions of X and Y don't agree")

   B0 <- matrix(0, p, K)
   LP0 <- matrix(0, N, K)

   pairs <- edges <- 0
   if(is.null(C)) {
      if(length(gamma) > 1 && gamma > 0) {
	 stop("gamma is non-zero but C is not set")
      }
      nE <- K * (K - 1) / 2
      C <- 0
      pairs <- 0
      edges <- 0
   } else {
      # mapping of edges to vertices (two vertices per edge), zero-based index
      # assumes that C is full size, i.e., all K(K-1)/2 edges are in it.
      # This is the same as cbind(C_J, ncol=2) from gennetwork()
      if(is(C, "matrix")) {
	 pairs <- t(apply(C, 1, function(r) which(r != 0))) - 1
	 # each kth column represents which edges task k is involved in
	 edges <- matrix(which(C != 0, arr.ind=TRUE)[,1], K - 1) - 1
      } else if(is(C, "list") 
	 && exists("pairs", where=C)
	 && exists("edges", where=C)) {
	 pairs <- C$pairs
	 edges <- C$edges
      } else {
	 stop("C is neither a matrix nor a list containing pairs and edges")
      }
   }

   # fit models in decreasing order of lambda, without messing with
   # the original ordering of lambda requested by user
   l1ord <- order(lambda, decreasing=TRUE)

   # We parallelise the gamma but not the lambda. Assuming that
   # the lambda are in increasing order, we can stop if there are no active
   # variables for a given gamma, as increasing lambda will only
   # result in no active variables again. This allows us to not waste time on
   # models that will have zero active variables anyway.
   Btmp <- lapply(seq(along=gamma), function(j) {
      lapply(seq(along=lambda2), function(m) {
	    Bjk <- vector("list", length(lambda))
	    LPjk <- vector("list", length(lambda))
	    nactive <- numeric(length(lambda))

	    # process sequential along the l1 penalty
	    for(i in seq(along=lambda))
	    {
	       if(verbose) {
		  cat("\t", "lambda[", l1ord[i], "] gamma[", j,
		     "] lambda2[", m, "] : ", sep="")
	       }

	       if(i == 1 || !warm) {
		  B <- B0
		  LP <- LP0
	       } else {
		  LP <- LPjk[[l1ord[i-1]]]
		  B <- Bjk[[l1ord[i-1]]]
	       }

	       r <- .C("fmpr",
		  as.numeric(X),	    # 1: X
		  as.numeric(Y),       	    # 2: Y
	          as.numeric(B),       	    # 3: B
		  as.numeric(LP),      	    # 4: LP
		  nrow(X),	       	    # 5: N
		  ncol(X),	       	    # 6: p
		  K,		       	    # 7: K
       	          rep(lambda[l1ord[i]], K), # 8: lambda
		  rep(lambda2[m], K),       # 9: lambda2
		  gamma[j],	       	    # 10: gamma
       	          as.numeric(C),            # 11: C
		  as.integer(pairs),        # 12: pairs
		  as.integer(edges),        # 13: edges
		  as.integer(maxiter),	    # 14: maxiter
       	          as.double(eps),      	    # 15: eps
		  as.integer(verbose), 	    # 16: verbose
		  integer(1),	       	    # 17: status
	          integer(1),	       	    # 18: iter
		  integer(1),    	    # 19: numactive
		  as.integer(divbyN),       # 20: divbyN
		  DUP=dup
	       )
	       status <- r[[17]]
	       numiter <- r[[18]]
	       nactive[l1ord[i]] <- r[[19]]
	       if(!status) {
	          cat("fmpr failed to converge within ",
	             maxiter, " iterations, lambda=", lambda[l1ord[i]],
		     "gamma=", gamma[j], "\n")
	       } else if(verbose) {
		  cat("fmpr converged in", numiter, "iterations",
	             "with", nactive[l1ord[i]], "active variables\n\n")
	       }
	          
	       # The linear predictor isn't sparse
	       LPjk[[l1ord[i]]] <- matrix(r[[4]], N, K)

	       if(sparse) {
		  Bjk[[l1ord[i]]] <- Matrix(matrix(r[[3]], p, K), sparse=TRUE)
	       } else {
		  Bjk[[l1ord[i]]] <- matrix(r[[3]], p, K)
	       }

	       if(!is.null(nzmax) && nactive[l1ord[i]] > nzmax) {
		  cat("fmpr reached maximum number of non-zero variables",
		     " (", nzmax, "): ", nactive[l1ord[i]], ", stopping\n",
		     sep="")
		  break
	       }

	    }

	    Bjk
      })
   })

   # Invert the list so that lambda is on the first dimension again
   B <- lapply(seq(along=lambda), function(i) {
	    lapply(seq(along=lambda2), function(m) {
	       lapply(seq(along=gamma), function(j) {
		  Btmp[[j]][[m]][[i]]
	       })
	    })
   })

   if(simplify 
      && length(lambda) == 1 
      && length(gamma) == 1 
      && length(lambda2) == 1) {
      B[[1]][[1]][[1]]
   } else {
      B
   }
}

# here, a lambda/lambda2/gamma vector represents applying a different penalty
# to each task
fmpr.single <- function(X, Y, lambda=0, lambda2=0, gamma=0, C=NULL,
      maxiter=1e5, eps=1e-6, verbose=FALSE, sparse=FALSE,
      nzmax=NULL, divbyN=TRUE, dup=TRUE)
{
   if(length(X) == 0 || length(Y) == 0)
      stop("X and/or Y have zero length")

   p <- ncol(X)
   Y <- cbind(Y)
   K <- ncol(Y)
   N <- nrow(X)

   if(length(gamma) > 1) {
      gamma <- gamma[1]
      warning("length(gamma) > 1, using first value only")
   } else if(length(gamma) == 0) {
      gamma <- 0
   }

   if(nrow(X) != nrow(Y))
      stop("dimensions of X and Y don't agree")

   pairs <- edges <- 0
   if(is.null(C)) {
      if(length(gamma) > 1 && gamma > 0) {
	 stop("gamma is non-zero but C is not set")
      }
      nE <- K * (K - 1) / 2
      C <- 0
      pairs <- 0
      edges <- 0
   } else {
      # mapping of edges to vertices (two vertices per edge), zero-based index
      # assumes that C is full size, i.e., all K(K-1)/2 edges are in it.
      # This is the same as cbind(C_J, ncol=2) from gennetwork()
      if(is(C, "matrix")) {
	 pairs <- t(apply(C, 1, function(r) which(r != 0))) - 1
	 # each kth column represents which edges task k is involved in
	 edges <- matrix(which(C != 0, arr.ind=TRUE)[,1], K - 1) - 1
      } else if(is(C, "list") 
	 && exists("pairs", where=C)
	 && exists("edges", where=C)) {
	 pairs <- C$pairs
	 edges <- C$edges
      } else {
	 stop("C is neither a matrix nor a list containing pairs and edges")
      }
   }

   B0 <- matrix(0, p, K)
   LP0 <- matrix(0, N, K)

   r <- .C("fmpr",
      as.numeric(X),	   # 1: X
      as.numeric(Y),       # 2: Y
      as.numeric(B0),      # 3: B
      as.numeric(LP0),     # 4: LP
      nrow(X),	       	   # 5: N
      ncol(X),	       	   # 6: p
      K,		   # 7: K
      as.numeric(lambda),  # 8: lambda
      as.numeric(lambda2), # 9: lambda2
      gamma,	       	   # 10: gamma
      as.numeric(C),       # 11: C
      as.integer(pairs),   # 12: pairs
      as.integer(edges),   # 13: edges
      as.integer(maxiter), # 14: maxiter
      as.double(eps),      # 15: eps
      as.integer(verbose), # 16: verbose
      integer(1),	   # 17: status
      integer(1),	   # 18: iter
      integer(1),    	   # 19: numactive
      as.integer(divbyN),  # 20: divbyN
      DUP=dup
   )
   status <- r[[17]]
   numiter <- r[[18]]
   nactive <- r[[19]]
   if(!status) {
      cat("fmpr failed to converge within ",
         maxiter, " iterations, lambda=", lambda,
         "gamma=", gamma, "\n")
   } else if(verbose) {
      cat("fmpr converged in", numiter, "iterations",
         "with", nactive, "active variables\n\n")
   }
      
   B <- if(sparse) {
      Matrix(matrix(r[[3]], p, K), sparse=TRUE)
   } else {
      matrix(r[[3]], p, K)
   }

   if(!is.null(nzmax) && nactive > nzmax) {
      cat("fmpr reached maximum number of non-zero variables",
         " (", nzmax, "): ", nactive, ", stopping\n",
         sep="")
      break
   }

   B
}

