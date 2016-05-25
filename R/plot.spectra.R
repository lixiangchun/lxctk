


plot.spectra <- function(d, frequency=TRUE, col=brewer.pal(dim(d)[2], "Dark2"), cex.tck=0.7, show.smp.name=TRUE, pdffig=NA, ...)
{
	if (!frequency) {
		d <- d / rowSums(d)
		par(mar=c(5, 4, 4, 7) + 0.1, lwd=0.1)
	} else {
		ylab = "Number of mutations"
		par(mar=c(5, 4, 4, 2) + 0.1, lwd=0.1)
	}
	mat = t(as.matrix(d))
	if (dim(mat)[1] == 8) {
		mat = mat[,order(colSums(mat), mat[8,],mat[7,], mat[6,], mat[5,], mat[4,], mat[3,], mat[2,], mat[1,], decreasing=TRUE)]
		#mat = mat[order(mat[8,], decreasing=TRUE),]
		#mat <- mat[do.call(order, as.data.frame(mat)),]
	} else {
		mat = mat[,order(colSums(mat), mat[6,], mat[5,], mat[4,], mat[3,], mat[2,], mat[1,], decreasing=TRUE)]
	}

	legend.text = gsub('\\.', ':', gsub('\\.\\.', '->', colnames(d)))
	legend.text = gsub("p:", "p*", legend.text)

	if (!is.na(pdffig)) { pdf(pdffig, width=8, height=6) }
	# Set xpd=TRUE to enable things to be drawn outside plot region
	# http://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
	p = barplot(mat, legend.text=NULL, space=0, axisnames=FALSE, col=col, las=2, ...)
	if (frequency) {
		legend("topright", legend=rev(legend.text), fill=rev(col), border='white', box.col='white', xpd=TRUE)
	} else {
		legend(max(p) + 0., 0.9, legend=rev(legend.text), fill=rev(col), border='white', box.col='white', xpd=TRUE, x.intersp=0.1)
	}
	if (show.smp.name) {
		#labels = sub("-VS-.*","",rownames(d))
		labels = rownames(d)
		axis(side=1, at=p, labels=labels, las=2, cex.axis=cex.tck, padj=0)
	}
	if (!is.na(pdffig)) { dev.off() }
	invisible(mat)
}

#spectrumplot("BC.EC.SNV.Genomic.muSpec")
#muSpec = read.table("BC.EC.SNV.Genomic.muSpec")
#muSpec8 = read.table("BC.EC.SNV.Genomic.muSpec8")
#save(muSpec, muSpec8, file="muSpec.RData")
#plot.spectra(df, frequency=T)
