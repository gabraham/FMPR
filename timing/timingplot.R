
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

load("timing.RData")

names(res1) <- Ns
names(res2) <- ps
names(res3) <- Ks

pdf("FMPR_timing_samples.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(Ns, t(res1m), xlab="Samples N", ylab="Time (sec)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
legend(x=0, y=max(res1m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2,
      lty=1:2)
dev.off()

pdf("FMPR_timing_samples_scaled.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(Ns, scale(t(res1m)), xlab="Samples N", ylab="Time (scaled)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
#legend(x=0, y=max(res1m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2)
dev.off()

pdf("FMPR_timing_variables.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(ps, t(res2m), xlab="Variables p", ylab="Time (sec)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
legend(x=50, y=max(res2m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2,
      lty=1:2)
dev.off()

pdf("FMPR_timing_variables_scaled.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(ps, scale(t(res2m)), xlab="Variables p", ylab="Time (scaled)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
#legend(x=50, y=max(res2m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2)
dev.off()

pdf("FMPR_timing_tasks.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(Ks, t(res3m), xlab="Tasks K", ylab="Time (sec)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
legend(x=0, y=max(res3m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2,
      lty=1:2)
dev.off()

pdf("FMPR_timing_tasks_scaled.pdf")
par(mar=c(5, 4.5, 4, 2) + 0.1)
matplot(Ks, scale(t(res3m)), xlab="Tasks K", ylab="Time (scaled)",
   cex=2, cex.axis=2, cex.lab=2, pch=20, type="b", lty=1:2, lwd=2)
#legend(x=0, y=max(res3m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2)
dev.off()

pdf("FMPR_timing.pdf", width=14, height=4)
lg  <- ""
par(mfrow=c(1, 3), mar=c(5, 4.5, 4, 2) + 0.1)
matplot(Ns, t(res1m), xlab="Samples N", ylab="Time (sec)", log=lg,
   cex=2, cex.axis=2, cex.lab=2, pch=c(20, 21), type="b", lty=1:2, lwd=2)
legend(x=0, y=max(res1m), legend=c("FMPR", "SPG"), col=1:2, lwd=5, cex=2,
      lty=1:2, pch=c(20, 21))
matplot(ps, t(res2m), xlab="Variables p", ylab="Time (sec)", log=lg,
   cex=2, cex.axis=2, cex.lab=2, pch=c(20, 21), type="b", lty=1:2, lwd=2)
matplot(Ks, t(res3m), xlab="Tasks K", ylab="Time (sec)", log=lg,
   cex=2, cex.axis=2, cex.lab=2, pch=c(20, 21), type="b", lty=1:2, lwd=2)
dev.off()

res1d <- melt(res1)
res2d <- melt(res2)
res3d <- melt(res3)

res1d$Var1 <- factor(as.character(res1d$Var1),
   labels=c("FMPR", "SPG / R", "SPG / MATLAB"))
res2d$Var1 <- factor(as.character(res2d$Var1),
   labels=c("FMPR", "SPG / R", "SPG / MATLAB"))
res3d$Var1 <- factor(as.character(res3d$Var1),
   labels=c("FMPR", "SPG / R", "SPG / MATLAB"))

res1d$L1 <- as.integer(res1d$L1)
res2d$L1 <- as.integer(res2d$L1)
res3d$L1 <- as.integer(res3d$L1)

m1 <- melt(with(res1d, tapply(value, list(Var1, L1), mean)))
m2 <- melt(with(res2d, tapply(value, list(Var1, L1), mean)))
m3 <- melt(with(res3d, tapply(value, list(Var1, L1), mean)))

pos <- ""
sz <- 0.8

g1 <- ggplot(res1d, aes(x=L1, y=value, colour=Var1, linetype=Var1))
g1 <- g1 + stat_summary(fun.data="mean_cl_normal", mult=2, size=sz,
   position=pos)
g1 <- g1 + theme_bw()
g1 <- g1 + scale_x_continuous("Samples N")
g1 <- g1 + scale_y_continuous("Time (sec)")
g1 <- g1 + theme(
   legend.position=c(0.28, 0.86),
   legend.title=element_blank(),
   legend.background=element_blank(),
   legend.key.width=unit(2.5, "lines")
)
g1 <- g1 + scale_colour_manual(values=c("black", "red", "blue"))
g1 <- g1 + geom_line(aes(x=Var2, y=value, colour=Var1, linetype=Var1), data=m1)

#pdf("FMPR_timing_samples_boxplot.pdf", width=8)
#print(g1)
#dev.off()

g2 <- ggplot(res2d, aes(x=L1, y=value, colour=Var1))
g2 <- g2 + stat_summary(fun.data="mean_cl_normal", mult=2, size=sz,
   position=pos)
g2 <- g2 + theme_bw()
g2 <- g2 + scale_x_continuous("Variables p")
g2 <- g2 + scale_y_continuous("Time (sec)")
g2 <- g2 + theme(legend.position="none")
g2 <- g2 + scale_colour_manual(values=c("black", "red", "blue"))
g2 <- g2 + geom_line(aes(x=Var2, y=value, colour=Var1, linetype=Var1), data=m2)

#pdf("FMPR_timing_variables_boxplot.pdf", width=10)
#print(g2)
#dev.off()

g3 <- ggplot(res3d, aes(x=L1, y=value, colour=Var1, linetype=Var1))
g3 <- g3 + stat_summary(fun.data="mean_cl_normal", mult=2, size=sz,
   position=pos)
g3 <- g3 + theme_bw()
g3 <- g3 + scale_x_continuous("Tasks K")
g3 <- g3 + scale_y_continuous("Time (sec)")
g3 <- g3 + theme(legend.position="none")
g3 <- g3 + scale_colour_manual(values=c("black", "red", "blue"))
g3 <- g3 + geom_line(aes(x=Var2, y=value, colour=Var1, linetype=Var1), data=m3)

#pdf("FMPR_timing_tasks_boxplot.pdf", width=8, height=6)
#print(g3)
#dev.off()

pdf("FMPR_timing_boxplot.pdf", height=3.3, width=10)
grid.arrange(g1, g2, g3, ncol=3, nrow=1)
dev.off()


