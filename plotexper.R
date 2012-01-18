
# Plot the experiment results

library(ROCR)
library(ggplot2)

pdf(sprintf("Expr%s.pdf", idv), width=12)

ymin <- 0
ymax <- 1
xmin <- 0
xmax <- 1

par(mfrow=c(1, 2))

plot(clean.rocr(res[[1]]$recovery$gr$roc), avg="threshold", col=1, main="ROC",
      xlim=c(xmin, xmax), ylim=c(ymin, ymax), lwd=3, cex=1.5, cex.axis=1.5,
      cex.lab=1.5)
plot(clean.rocr(res[[1]]$recovery$lasso$roc), avg="threshold", add=TRUE,
      col=2, lwd=3, lty=2)
plot(clean.rocr(res[[1]]$recovery$ridge$roc), avg="threshold", add=TRUE,
      col=3, lwd=3, lty=3)
plot(clean.rocr(res[[1]]$recovery$elnet.fmpr$roc), avg="threshold", add=TRUE,
      col=4, lwd=3, lty=4)
plot(clean.rocr(res[[1]]$recovery$elnet.glmnet$roc), avg="threshold", add=TRUE, col=5,
      lwd=3, lty=5)

legend(xmin, ymin + 0.25,
   c("FMPR", "Lasso", "Ridge", "ElasticNet-FMPR", "ElasticNet-glmnet"),
   col=1:5, lwd=3, lty=1:5)

plot(clean.rocr(res[[1]]$recovery$gr$prc), avg="threshold", col=1,
      main="Precision-Recall", ylim=c(ymin, ymax), xlim=c(xmin, xmax), lwd=3,
      cex=1.5, cex.axis=1.5, cex.lab=1.5)
plot(clean.rocr(res[[1]]$recovery$lasso$prc), avg="threshold", col=2,
      add=TRUE, lwd=3, lty=2)
plot(clean.rocr(res[[1]]$recovery$ridge$prc), avg="threshold", add=TRUE,
      col=3, lwd=3, lty=3)
plot(clean.rocr(res[[1]]$recovery$elnet.fmpr$prc), avg="threshold", add=TRUE,
      col=4, lwd=3, lty=4)
plot(clean.rocr(res[[1]]$recovery$elnet.glmnet$prc), avg="threshold",
      add=TRUE, col=5, lwd=3, lty=5)

legend(xmin, ymin + 0.25,
   c("FMPR", "Lasso", "Ridge", "ElasticNet-FMPR", "ElasticNet-glmnet"),
   col=1:5, lwd=3, lty=1:5)

dev.off()

t.test(res[[1]]$R2[, 1], res[[1]]$R2[,2])

mytheme <- function(base_size=10)
{
   structure(list(
	 axis.text.x=theme_text(size=20),
	 axis.text.y=theme_text(size=20, hjust=1),
	 axis.title.x=theme_text(size=23),
	 axis.title.y=theme_text(size=23),
	 axis.ticks=theme_blank(),
	 plot.title=theme_text(size=30),
	 legend.text=theme_text(size=20),
	 legend.title=theme_text(size=20, hjust=0),
	 legend.key.size=unit(2, "lines"),
	 legend.background=theme_rect(col=0, fill=0),
	 legend.key=theme_blank()
   ), class="options")
}

r2 <- melt(res[[1]]$R2)
colnames(r2) <- c("Replication", "Method", "R2")
g <- ggplot(r2, aes(Method, R2))
g <- g + geom_boxplot()
g <- g + scale_y_continuous(expression(R^2))
g <- g + theme_bw() + mytheme()

pdf(sprintf("Expr%s_R2.pdf", idv), width=12)
print(g)
dev.off()

