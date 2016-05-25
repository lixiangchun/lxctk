

plot.depth <- function(df, 
	left.jitter.amount=0.15, left.bp.notch=T, left.bp.ylim=c(0,300), left.bp.outline=F, left.bp.main="CCDS", left.bp.cex.main=0.8, left.bp.axis.label='CCDS',
	right.jitter.amount=0.3, right.bp.notch=T, right.bp.ylim=c(0.7,1), right.bp.outline=F, right.bp.main="Cov", right.bp.cex.main=0.8, right.bp.axis.labels=NA,
	left.bp.point.col=rgb(0.9, 0.5, 0.5, alpha=0.8),
	right.bp.point.col=rgb(0.9, 0.4, 0.5, alpha=0.8),
	pdffile=NA)
{

df$V1 <- NULL
df$V2 <- NULL
depth <- df$V7
df$V7 <- NULL

y = jitter(rep(1, length(depth)), amount=left.jitter.amount)

x <- stack(df)

if (!is.na(pdffile)) { pdf(pdffile) }
par(las=1)
layout(matrix(c(1,2), byrow=T, nrow=1), widths=c(30,50))
boxplot(depth, ylab="Depth", notch=left.bp.notch, ylim=left.bp.ylim, axes=F, outline=left.bp.outline, xaxs='i', yaxs='i', main=left.bp.main, cex.main=left.bp.cex.main)
points(y, depth, pch=20, col=left.bp.point.col, xlab="", ylab="")
axis(side=2)
axis(side=1, at=c(0.5, 1, 1.5), labels=c("",left.bp.axis.label,""))

text(x=1.4, y=median(depth), format(median(depth), digits=2), col="blue")


par(mar=c(5, 2, 4, 2)+0.1)
boxplot(values ~ ind, data=x, axes=F, ylim=right.bp.ylim, outline=right.bp.outline, xaxs='i', yaxs='i', notch=right.bp.notch, main=right.bp.main, cex.main=right.bp.cex.main)
y = jitter(rep(1, length(df$V3)), amount=right.jitter.amount)
points(y, df$V3, pch=20, col=right.bp.point.col, cex=0.3, xlab="", ylab="")
y = jitter(rep(2, length(df$V4)), amount=right.jitter.amount)
points(y, df$V4, pch=20, col=right.bp.point.col, cex=0.4, xlab="", ylab="")
y = jitter(rep(3, length(df$V5)), amount=right.jitter.amount)
points(y, df$V5, pch=20, col=right.bp.point.col, cex=0.6, xlab="", ylab="")
y = jitter(rep(4, length(df$V6)), amount=right.jitter.amount)
points(y, df$V6, pch=20, col=right.bp.point.col, cex=0.8, xlab="", ylab="")

if (any(is.na(right.bp.axis.labels))) {
	right.bp.axis.labels <- c(">=1X", ">=4X", ">=10X", ">=20X")
}
axis(side=1, at=1:4, labels=right.bp.axis.labels)
axis(side=2)
mtext(text="Fraction", side=2, las=3, padj=-5)

# set xpd = TRUE to allow operation out of canvas
#par(xpd=TRUE)
text(x=1.6, y=median(df$V3), format(median(df$V3), digits=2), col='blue', xpd=TRUE)
text(x=2.6, y=median(df$V4), format(median(df$V4), digits=2), col='blue', xpd=TRUE)
text(x=3.6, y=median(df$V5), format(median(df$V5), digits=2), col='blue', xpd=TRUE)
text(x=4.6, y=median(df$V6), format(median(df$V6), digits=2), col='blue', xpd=TRUE)

if (!is.na(pdffile)) { dev.off() }

}
