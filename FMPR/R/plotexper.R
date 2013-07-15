
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
   theme(
	 axis.text.x=element_text(size=24, angle=-90, hjust=0),
	 axis.text.y=element_text(size=24, hjust=1),
	 axis.title.x=element_text(size=28),
	 axis.title.y=element_text(size=28),
	 axis.ticks=element_blank(),
	 plot.title=element_text(size=30),
	 legend.text=element_text(size=20),
	 legend.title=element_text(size=20, hjust=0),
	 legend.key.size=unit(2, "lines"),
	 legend.background=element_rect(colour=0, fill=0),
	 legend.key=element_blank()
   )
}

plot.exper <- function(x,
   lim=list(roc=c(0.95, 1, 0, 1), prc=c(0, 1, 0, 1)), ci=FALSE,
      models.rec=NULL, models.R2=NULL,
      cex=2, suffix="", col=NULL,
      legend.x=NULL, legend.y=NULL, plot.prc=TRUE, plot.legend=TRUE,
      main=NULL, outlier.shape=1, ...)
{
   if(is.null(col)) {
      # remove yellow, make it black
      col <- rainbow(length(x$recovery))
      col[3] <- 1
   }

   if(plot.prc) {
      width <- 12
      sq <- 1:2
   } else {
      sq <- 1
      width <- 6
   }

   pdf(sprintf("%s/%s%s.pdf", x$dir, x$dir, suffix), width=width)
   par(mfrow=c(1, length(sq)), mar=c(4, 4.5, 2, 0.5) + 0.1, pty="s")
   usr <- par("usr")

   nm <- c("roc", "prc")
   titles <- c("ROC", "PRC")
   xlab <- c("Specificity", "Recall")
   ylab <- c("Sensitivity", "Precision")

   if(ci) {
      spread.estimate <- "stderror"
      spread.scale <- 2
   } else {
      spread.estimate <- "none"
      spread.scale <- 1
   }

   if(is.null(legend.x)) {
      legend.x <- usr[2] - 0.0275
   }
   if(is.null(legend.y)) {
      legend.y <- usr[4] + 0.02
   }

   for(i in sq)
   {
      maint <- ifelse(is.null(main), titles[i], main)
      plot(NULL, main=maint, xlim=lim[[i]][1:2], ylim=lim[[i]][3:4],
         cex=cex, cex.axis=cex, cex.lab=cex, xlab=xlab[i], ylab=ylab[i],
	 cex.main=cex)
      
      rr <- if(is.null(models.rec)) {
	 z <- seq(along=x$recovery)
	 names(z) <- names(x$recovery)
	 z
      } else {
	 names(models.rec) <- models.rec
	 models.rec
      }

      for(j in seq(along=rr))
      {
	 plot.rocr(x$recovery[[rr[j]]][[nm[i]]], avg="threshold", add=TRUE,
	    col=col[j], lwd=5, lty=j,
	    spread.estimate=spread.estimate,
	    spread.scale=spread.scale)
      }

      if(plot.legend && i == 1) {
	 legend(x=legend.x, y=legend.y,
      	    legend=names(rr),
      	    col=col, lty=1:length(rr),
	    lwd=6, bty="n", cex=1.6, seg.len=3
      	 )
      } 
   }

   dev.off()
   
   r2 <- melt(x$R2)
   m <- as.character(r2[,2])
   r2[, 2] <- factor(m)
   colnames(r2) <- c("Replication", "Method", "R2")

   if(!is.null(models.R2)) {
      r2 <- droplevels(subset(r2, Method %in% models.R2))
      r2$Method <- factor(as.character(r2$Method), levels=models.R2)
   }

   g <- ggplot(r2, aes(Method, R2))
   g <- g + scale_y_continuous(expression(R^2))
   g <- g + theme_bw() + mytheme()
   g <- g + stat_summary(fun.data="mean_cl_normal", size=1.5)
   g <- g + ggtitle(main)
   
   pdf(sprintf("%s/%s%s_R2.pdf", x$dir, x$dir, suffix))
   print(g)
   dev.off()
}

