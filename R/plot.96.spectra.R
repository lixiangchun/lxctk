

plot.96.spectra <- function(x, col=brewer.pal(6, "Paired"), main="spectra.96", cex.main=1, context.col=col, draw.legend=TRUE, bottom=FALSE)
{
	if (any(is.na(col))) {
		col = brewer.pal(6, 'Dark2')
		tmp = col[3]
		col[3] = col[4]
		col[4] = tmp
	} else {
		col = col
	}
	if (!is.numeric(x)) {
		x = as.numeric(x)
	}
	if (any(x > 1)) {
		x = x / sum(x)
	}

	if (!bottom) {
		par(mar=c(1,0.5,2,0.5))
	} else {
		par(mar=c(4,0.5,2,0.5))
	}
	#rect(0, 0, 96, 0.105, lty="dashed")
	#p = barplot(x, border=NA, ylim=c(0,0.2), space=0, col=rep(col, each=16), axes=F, main=main)
	p = barplot(x, border=NA, ylim=c(0,0.2), space=0, col=NA, axes=F, main="")
	segments(0,0, 0, 0.13, col="gray", lwd=1.2)
	segments(0,0, 96,0, col="gray", lwd=1.2)
	segments(0, 0.13/2, 96, 0.13/2, col="gray", lty="dashed", lwd=0.6)

	y <- seq(0, 0.13, length.out=3)	
	text(x=c(0,0,0), y=y, labels=as.character(y * 100), xpd=T, adj=1.2)

	par(new=T)
	p = barplot(x, border=NA, ylim=c(0,0.2), space=0, col=rep(col, each=16), axes=F, main=main, cex.main=cex.main)
	#axis(side=2, at=c(0, 0.05, 0.1), labels=c("0","5%","10%"), las=2, hadj=0.7)
	#abline(v=0, col='grey20')
	#abline(h=0.126, col="grey20")
	texts = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
	texts.x = c(8.5, 24.5, 40.5, 56.5, 72.5, 88.5)
	if (draw.legend && !bottom) {
		rect.yleft = 0.15
		YLIMMAX = 0.16
		rect(p[1]-0.5, rect.yleft, p[16], YLIMMAX, col=col[1], border=NA)
		rect(p[17], rect.yleft, p[32], YLIMMAX, col=col[2], border=NA)
		rect(p[33], rect.yleft, p[48], YLIMMAX, col=col[3], border=NA)
		rect(p[49], rect.yleft, p[64], YLIMMAX, col=col[4], border=NA)
		rect(p[65], rect.yleft, p[80], YLIMMAX, col=col[5], border=NA)
		rect(p[81], rect.yleft, p[96]+0.5, YLIMMAX, col=col[6], border=NA)
		texts.y = rep(YLIMMAX + 0.015, 6)
		text(texts.x, texts.y, texts, cex=1, xpd=TRUE, col=context.col, font=2)
	}
	if (draw.legend && bottom) {
		rect.yleft = -0.01
		YLIMMAX = -0.038
		rect(p[1]-0.5, rect.yleft, p[16], YLIMMAX, col=col[1], border=NA, xpd=T)
		rect(p[17], rect.yleft, p[32], YLIMMAX, col=col[2], border=NA, xpd=T)
		rect(p[33], rect.yleft, p[48], YLIMMAX, col=col[3], border=NA, xpd=T)
		rect(p[49], rect.yleft, p[64], YLIMMAX, col=col[4], border=NA, xpd=T)
		rect(p[65], rect.yleft, p[80], YLIMMAX, col=col[5], border=NA, xpd=T)
		rect(p[81], rect.yleft, p[96]+0.5, YLIMMAX, col=col[6], border=NA, xpd=T)
		texts.y = rep(-0.055, 6)
		text(texts.x, texts.y, texts, cex=1.5, xpd=TRUE, col=context.col, font=2)
	}
}

#d = read.table("sampleid_mapping")
#sampleid=as.character(d$V2)
#names(sampleid) = as.character(d$V1)

#d=read.table("spectra.txt")
#rownames(d) = sampleid[rownames(d)]

#pdf("spectra.pdf", width=9, height=7)
#par(mfrow=c(6,2), family="sans")

#x = d["1-1",]
#plot.96.spectrum(x, main="1-1", bottom=F)
#x = d["1-2",]
#plot.96.spectrum(x, main="1-2", bottom=F)

#x = d["2-1",]
#plot.96.spectrum(x, main="2-1", draw.legend=T)
#x = d["2-2",]
#plot.96.spectrum(x, main="2-2", draw.legend=T)

#x = d['3-1',]
#plot.96.spectrum(x, main="3-1", draw.legend=T)
#x = d["3-2",]
#plot.96.spectrum(x, main="3-2", draw.legend=T)

#x = d['4-1',]
#plot.96.spectrum(x, main="4-1", bottom=F)
#x = d["4-2",]
#plot.96.spectrum(x, main="4-2", bottom=F)

#x = d['5-1',]
#plot.96.spectrum(x, main="5-1", bottom=F)
#x = d["5-2",]
#plot.96.spectrum(x, main="5-2", bottom=F)

#x = d['6-1',]
#plot.96.spectrum(x, main="6-1", bottom=F)
#x = d["6-2",]
#plot.96.spectrum(x, main="6-2", bottom=F)
#

#dev.off()
