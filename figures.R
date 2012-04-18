library(FMPR)
library(ggplot2)

N <- 100
p <- 30
K <- 10 
w <- 1
#B0 <- getB(p, K, w=w, type="same")
B1 <- getB(p, K, w=NULL, type="sparsity", mean=w, sd=0.5)
B2 <- getB(p, K, w=NULL, type="random", mean=w, sd=0.5)
B0 <- sign(B1) * w

trans <- function(X)
{
   (X - min(X)) / (max(X) - min(X))
}

m0 <- data.frame(melt(trans(B0)), Type="Same sparsity, same weights", B="weights")
m1 <- data.frame(melt(trans(B1)), Type="Same sparsity, different weights",
   B="weights")
m2 <- data.frame(melt(trans(B2)), Type="Unrelated", B="weights")

X <- matrix(rnorm(N * p), N, p)
E <- matrix(rnorm(N * K, 0, 1), N, K)

Y0 <- X %*% B0 + E
Y1 <- X %*% B1 + E
Y2 <- X %*% B2 + E

r0 <- data.frame(melt(cor(Y0)), Type="Same sparsity, same weights",
   B="correlation")
r1 <- data.frame(melt(cor(Y1)), Type="Same sparsity, different weights",
   B="correlation")
r2 <- data.frame(melt(cor(Y2)), Type="Unrelated",
   B="correlation")

m <- rbind(m0, m1, m2, r0, r1, r2)
g <- ggplot(m, aes(X2, X1)) + geom_tile(aes(fill=value)) 
g <- g + facet_grid(B ~ Type, scales="free_y")
g <- g + theme_bw()
g <- g + scale_fill_gradient(name=expression(beta), low="white", high="black")
g <- g + scale_x_continuous("Tasks k")
g <- g + scale_y_continuous("Variables p")
g <- g + opts(legend.position="none")

pdf("sparsity.pdf", width=10)
print(g)
dev.off()


