R2 <- function(pr, y) 1 - sum((pr - y)^2) / sum((y - mean(y))^2)

crossval <- function(X, Y, nfolds=5, fun, ...)
{
   N <- nrow(X)
   Y <- cbind(Y)
   folds <- sample(1:nfolds, N, TRUE)
   s <- sapply(1:nfolds, function(fold) {
      g <- fun(X[folds != fold, ], Y[folds != fold, ], ...)
      p <- X[folds == fold, ] %*% g
      R2(as.numeric(p), as.numeric(Y[folds == fold, ]))
   })
   mean(s)
}

