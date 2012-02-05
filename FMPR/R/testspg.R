
# 
# Tests the C implementation of SPG versus the MATLAB version, using random
# datasets.
# 
# The correct path to the MATLAB SPG code must be set.
#
# Returns statistics:
# Cmse: mean squared error of C matrix
# Csign: mean error in sign 
# Bmse: mean squared error of weights B
#
# The smaller the statistics, the better the agreement.
run.spg.test <- function(N, p, K, spgpath)
{
   X <- scale(matrix(rnorm(N * p), N, p))
   Y <- scale(matrix(rnorm(N * K), N, K))
   
   write.table(X, file="X.txt", col.names=FALSE, row.names=FALSE)
   write.table(Y, file="Y.txt", col.names=FALSE, row.names=FALSE)
   
   gamma <- runif(1, 0, 1e3)
   lambda <- runif(1, 0, 1e3)
   tol <- 1e-4
   maxiter <- 1e4
   mu <- 1e-6
   R <- abs(cor(Y))
   corthresh <- if(K <= 2) {
      runif(1, 0, R[1, 2])
   } else {
      # ensure threshold selects at least one edge, and add small dithering
      # factor to account for loss of precision in translating floating point
      # to ascii and back
      median(R[upper.tri(R)]) + rnorm(1, 0, 1e-3)
   }
   cortype <- 1

   C <- gennetwork(Y, corthresh=corthresh, cortype=cortype)
   g <- spg(X, Y, C, lambda=lambda, gamma=gamma, mu=mu,
      tol=tol, maxiter=maxiter, simplify=TRUE, verbose=TRUE)
   
   s <- paste(
      c(
         paste("addpath '", spgpath, "';"),
         "X = dlmread('X.txt');",
         "Y = dlmread('Y.txt');",
         paste("option.cortype = ", cortype, ";"),
         paste("option.corthreshold = ", corthresh, ";"),
         paste("option.maxiter = ", maxiter, ";"),
         paste("option.tol = ", tol, ";"),
         "option.verbose = true;",
         "option.display_iter = 1;",
         paste("option.mu = ", mu, ";"),
         paste("lambda = ", lambda, ";"),
         paste("gamma = ", gamma, ";"),
         "[C, CNorm, E, Ecoef, Esign, R] = gennetwork(Y, option);",
         paste("[beta, obj, density, iter, time] =",
         "SPG_multi(Y, X, gamma, lambda, C, CNorm, option);"),
         "save('beta.txt', 'beta', '-ascii');",
         "Cf = full(C);",
         "save('C.txt', 'Cf', '-ascii');",
         "quit;"
      ),
      collapse="\n"
   )
   
   cat(s, file="run.m")
   system("matlab -nodisplay -r run")
   
   beta <- matrix(scan("beta.txt"), p, K, byrow=TRUE)
   Cm <- as.matrix(read.table("C.txt", header=FALSE))

   if(nrow(Cm) != nrow(C))
      stop("C matrices don't agree in dimension")

   c(
      Cmse=mean((C - Cm)^2),
      Csign=mean(sign(C) != sign(Cm)),
      Bmse=mean((beta - g)^2)
   )
}

spg.test <- function(nreps=100, spgpath="~/Software/SPG_Multi_Graph")
{
   res <- sapply(1:nreps, function(i) {
      N <- sample(1e3, 1)
      p <- sample(5e2, 1)
      K <- sample(1e2, 1) + 1 # ensures K>=2
      run.spg.test(N, p, K, spgpath)
   })
   apply(res, 1, max)
}

