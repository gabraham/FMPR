makematlab <- function(spg.path, xfile, yfile, cortype, corthresh, maxiter,
   tol, mu, lambda, gamma)
{
   s <- paste(
      c(
         paste("addpath '", spg.path, "';", sep=""),
         paste("X = dlmread('", xfile, "');", sep=""),
         paste("Y = dlmread('", yfile, "');", sep=""),
	 paste("size(X)"),
	 paste("size(Y)"),
         paste("option.cortype = ", cortype, ";"),
         paste("option.corthreshold = ", corthresh, ";"),
         paste("option.maxiter = ", maxiter, ";"),
         paste("option.tol = ", tol, ";"),
         "option.verbose = false;",
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
         "save('B.txt', 'beta', '-ascii');",
	 "fprintf('time:%d\\n', time)",
         "quit;"
      ),
      collapse="\n"
   )
   s
}
