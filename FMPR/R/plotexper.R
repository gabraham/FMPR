
# Remove ROC/PRC replications that had non-sensical results
clean.rocr <- function(obj)
{
   if(is.null(obj))
      return(NULL)
   n <- length(obj@x.values)
   if(n == 0)
      return(NULL)
   len <- sapply(obj@x.values, length)
   w <- len > 2
   obj@x.values <- obj@x.values[w]
   obj@y.values <- obj@y.values[w]
   obj@alpha.values <- obj@alpha.values[w]
   obj
}

# Plot the experiment results

plot.rocr <- function(obj, ...)
{
   obj <- clean.rocr(obj)
   l <- list(...)
   if(is.null(l$add)) {
      plot(obj, ...)
   } else if(!is.null(obj)) {
      plot(obj, ...)
   }
}

mytheme <- function(base_size=10)
{
   structure(list(
	 axis.text.x=theme_text(size=24, angle=-30),
	 axis.text.y=theme_text(size=24, hjust=1),
	 axis.title.x=theme_text(size=28),
	 axis.title.y=theme_text(size=28),
	 axis.ticks=theme_blank(),
	 plot.title=theme_text(size=30),
	 legend.text=theme_text(size=20),
	 legend.title=theme_text(size=20, hjust=0),
	 legend.key.size=unit(2, "lines"),
	 legend.background=theme_rect(col=0, fill=0),
	 legend.key=theme_blank()
   ), class="options")
}

plot.exper <- function()
{
   pdf(sprintf("Expr%s.pdf", idv), width=12)
   
   ymin <- 0
   ymax <- 1
   xmin <- 0
   xmax <- 1
   
   par(mfrow=c(1, 2), mar=c(5, 4, 3, 0.1) + 0.1)

   nm <- c("roc", "prc")
   titles <- c("ROC", "PRC")
   xlab <- c("Specificity", "Recall")
   ylab <- c("Sensitivity", "Precision")
   
   for(i in 1:2)
   {
      plot.rocr(res[[1]]$recovery$fmpr.w1[[nm[i]]], avg="threshold", col=1,
            main=titles[i], xlim=c(xmin, xmax), ylim=c(ymin, ymax), lwd=3,
            cex=1.5, cex.axis=1.5, cex.lab=1.5, lty=1,
            xlab=xlab[i], ylab=ylab[i])
      plot.rocr(res[[1]]$recovery$fmpr.w2[[nm[i]]], avg="threshold", add=TRUE,
            col=2, lwd=3, lty=1)
      plot.rocr(res[[1]]$recovery$spg.w1[[nm[i]]], avg="threshold", add=TRUE,
            col=3, lwd=3, lty=2)
      plot.rocr(res[[1]]$recovery$spg.w2[[nm[i]]], avg="threshold", add=TRUE,
            col=4, lwd=3, lty=2)
      plot.rocr(res[[1]]$recovery$lasso[[nm[i]]], avg="threshold", add=TRUE,
            col=5, lwd=3, lty=3)
      plot.rocr(res[[1]]$recovery$ridge[[nm[i]]], avg="threshold", add=TRUE,
            col=6, lwd=3, lty=4)
      plot.rocr(res[[1]]$recovery$elnet.glmnet[[nm[i]]], avg="threshold",
	    add=TRUE, col=7, lwd=3, lty=5)
      
      legend(xmin, ymin + 0.3,
         c("FMPR-w1", "FMPR-w2", "GFlasso-w1", "GFlasso-w2",
	    "Lasso", "Ridge", "ElasticNet"),
         col=1:7, lwd=3, lty=c(1, 1, 2, 2, 3, 4, 5)
      )
   }

   dev.off()
   
   r2 <- melt(res[[1]]$R2)
   m <- as.character(r2[,2])
   #m[m == "GRt"] <- "FMPR-t"
   #m[m == "GRw1"] <- "FMPR-w1"
   #m[m == "GRw2"] <- "FMPR-w2"
   #m[m == "SPGt"] <- "GFlasso-t"
   m[m == "SPGw1"] <- "GFlasso-w1"
   m[m == "SPGw2"] <- "GFlasso-w2"
   r2[, 2] <- factor(m)
   colnames(r2) <- c("Replication", "Method", "R2")
   g <- ggplot(r2, aes(Method, R2))
   g <- g + geom_boxplot()
   g <- g + scale_y_continuous(expression(R^2))
   g <- g + theme_bw() + mytheme()
   
   pdf(sprintf("Expr%s_R2.pdf", idv), width=12)
   print(g)
   dev.off()
   
}
