
## Note: A simple R script to produce high quality figure of mutational landscape
##+for paper, developed and maintained by Li Xiangchun in BGI, all rights reserved.
## It's highly recommended NOT to include '-' in sample names.

plot.land <- function(upper.panel, middle.panel, left.panel, right.panel, legend.panel, subtype.panel=NULL, subtype.2.panel=NULL,
							   output.fig.pdffile=NULL,
							   left.panel.xlim=NULL,
							   family.font="sans", 
							   syn.nonsyn.legend.color=NULL,
							   syn.nonsyn.legend.ncol=NULL,
							   syn.nonsyn.legend.loc="top",
							   image.color=NULL,

							   subtype.labels=NULL,
							   subtype.label.cex=1,
							   subtype.panel.color=brewer.pal(7, "Dark2"),
							   subtype.panel.identifier="",
							   subtype.panel.identifier.cex=0.6,
							   subtype.legend.text.width=NULL,

							   subtype.2.labels=NULL,
							   subtype.2.label.cex=1,
							   subtype.2.panel.color=brewer.pal(7, "Set2"),
							   subtype.2.panel.identifier="",
							   subtype.2.panel.identifier.cex=0.6,
							   subtype.2.legend.text.width=NULL,

							   right.panel.color=brewer.pal(7, "Purples")[5],
							   sort.by.qvalue=FALSE,
							   group.by.subtype=FALSE,
							   left.panel.cex.axis=1.5,
							   right.panel.cex.axis=1.3,
							   upper.panel.mar=c(0, 1.3, 0.5, 1.16),
							   legend.panel.mar=c(0, 1, 2, 1),
							   legend.panel.ncol=NULL,
							   legend.panel.y=-0.05,
							   left.panel.mar=c(5.5, 1.5, 0.2, 2.9),
							   middle.panel.mar=c(5.5, 1.3, 0.1, 1.16),
							   right.panel.mar=c(5.5, 3, 0.2, 1.5),
							   left.panel.xlab="Individuals\nwith mutations",
							   left.panel.xlab.cex=0.8,
							   right.panel.xlab="-log10\n(q-value)",
							   right.panel.xlab.cex=0.8,
							   right.panel.qscore.cutoff=1,
							   right.panel.qscore.lwd=1,
							   right.panel.qscore.lty='dashed',
							   right.panel.qscore.col='red')
{
	upper.panel.samples <- as.character(upper.panel[,1])
	middle.panel.samples <- as.character(colnames(middle.panel))
	left.panel.genes <- as.character(left.panel[,1])
	right.panel.genes <- as.character(right.panel[,1])
	
	if (length(upper.panel.samples) != length(middle.panel.samples)) {
		stop("The number of samples in upper panel is NOT equal to the number of samples in middle panel!")
	}
	if (!all(sort(upper.panel.samples) == sort(middle.panel.samples))) {
		stop("Sample names in upper panel are not all of the same as those in middle panel!")
	}

	# Some constant variable
	CEXeq1.3=1.3
	CEXeq1.1=1.1
	CEXeq1.0=1
	CEXeq0.8=0.8
	SPACEeq0.06=0.06

	if (!is.null(output.fig.pdffile)) { pdf(output.fig.pdffile, width=6.5, height=sqrt(nrow(left.panel)), family=family.font) }
	else {par(family=family.font)}
	if (is.null(subtype.panel)) {
		layout(mat=matrix(c(1, 2, 3, 4, 5, 6), nrow=2, byrow=TRUE), widths=c(40, 120, 40), heights=c(40, 80))
	} else {
		if (is.null(subtype.2.panel)) {
			layout(mat=matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow=3, byrow=TRUE), widths=c(40, 120, 40), heights=c(40, 5, 80))
		} else {
			layout(mat=matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow=4, byrow=TRUE), widths=c(40, 120, 40), heights=c(40, 3, 3, 80))
		}
	}
	### Set colors for legend, e.g. Purples, BuGn, BuPu
	if (is.null(syn.nonsyn.legend.color)) {
		syn.nonsyn.legend.color <- brewer.pal(3, "Purples")[2:3]
	} else {
		syn.nonsyn.legend.color <- syn.nonsyn.legend.color
	}

	if (sort.by.qvalue) {
		right.panel <- right.panel[order(right.panel[,2], decreasing=F),]
		genes <- right.panel[,1]
	} else {                 # Sorting by mutation frequency
		left.panel <- left.panel[order(left.panel[,2], decreasing=T),]
		genes <- left.panel[,1]
	}
	genes <- as.character(genes)

	middle.panel <- middle.panel[genes,]
	rownames(left.panel) <- as.character(left.panel[,1])	
	rownames(right.panel) <- as.character(right.panel[,1])
	right.panel <- right.panel[genes,]
	left.panel <- left.panel[genes,]

	middle.panel <- t(middle.panel)
	# A dynamic way to order all of a matrix column, refer to following web link:
	# http://stackoverflow.com/questions/12077413/order-a-matrix-by-multiple-column-in-r
	middle.panel <- middle.panel[do.call(order, as.data.frame(-middle.panel)),] 
	middle.panel <- middle.panel[,rev(genes)]

	sample.order <- rownames(middle.panel)
	subtype.categs <- sort(unique(subtype.panel[,2]))
	subtype.categs.num <- length(subtype.categs)
	
	if (!is.null(subtype.panel)) {
		rownames(subtype.panel) <- subtype.panel[,1]
		subtype.data <- as.integer(subtype.panel[,2])
		names(subtype.data) <- as.character(subtype.panel[,1])
		subtype.data <- subtype.data[sample.order]
	}
	if (!is.null(subtype.2.panel)) {
		rownames(subtype.2.panel) <- subtype.2.panel[,1]
		subtype.2.data <- as.integer(subtype.2.panel[,2])
		names(subtype.2.data) <- as.character(subtype.2.panel[,1])
		subtype.2.data <- subtype.2.data[sample.order]
	}

	if (group.by.subtype) {
		out <- list()
		for (i in 1:subtype.categs.num) {
			subtype.tmp.mat <- middle.panel[names(subtype.data[subtype.data==subtype.categs[i]]),]
			subtype.tmp.mat <- subtype.tmp.mat[do.call(order, as.data.frame(-subtype.tmp.mat)),] 
			out[[i]] <- subtype.tmp.mat
		}
		middle.panel <- do.call('rbind', out)
		sample.order <- rownames(middle.panel)
		subtype.data <- subtype.data[sample.order]
		if (!is.null(subtype.2.panel)) {
			subtype.2.data <- subtype.2.data[sample.order]
		}
	}

	### Plot an empty figure and legend for Syn. and Non-syn. mutations
	plot.empty.graph <- function() plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
	par(mar=c(0, 0, 0, 0))
	plot.empty.graph()
	#legend(x=legend.syn.nonsyn.x,y=legend.syn.nonsyn.y, legend=c("Syn.", "Non syn."), fill=syn.nonsyn.legend.color, box.col=NA, cex=1.3, border=syn.nonsyn.legend.color, x.intersp=0.2)
	if (is.null(syn.nonsyn.legend.ncol)) {
		legend('bottom', legend=c("Syn.", "Non syn."), fill=syn.nonsyn.legend.color, box.col=NA, cex=CEXeq1.3, border=syn.nonsyn.legend.color, x.intersp=0.2)
	}

	# Configure and plot for the upper panel
	#par(mar=c(0, 0.16, 3, 0))
	par(mar=upper.panel.mar)
	rownames(upper.panel) <- as.character(upper.panel[,1])
	upper.panel <- upper.panel[sample.order,]
	rownames(upper.panel) <- NULL  # Remove row names to prevent automatic barplot labeling

	#par(mgp=c(0,-0.4,-1.3))
	upper.panel.dat = t(upper.panel[,c(3,2)]/(32 * 1))
	barplot(upper.panel.dat, border=NA, col=rev(syn.nonsyn.legend.color), las=2, axes=F, space=SPACEeq0.06, xaxs='i')
	# Use tcl to control tick length, lwd to set axis width, col to set axis color
	axis(side=2, cex.axis=CEXeq1.0, las=2, tcl=-0.3, lwd=0.6, col='grey20')
	mtext(text="Mutations/Mb", side=2, padj=-3.5, cex=0.75)
	if (!is.null(syn.nonsyn.legend.ncol)) {
		if (syn.nonsyn.legend.ncol == 1) {
			legend(syn.nonsyn.legend.loc, legend=c("Syn.", "Non syn."), fill=syn.nonsyn.legend.color, box.col=NA, cex=CEXeq1.3, border=syn.nonsyn.legend.color, x.intersp=0.2, horiz=F)
		}
		else {
			legend(syn.nonsyn.legend.loc, legend=c("Syn.", "Non syn."), fill=syn.nonsyn.legend.color, box.col=NA, cex=CEXeq1.3, border=syn.nonsyn.legend.color, x.intersp=0.2, horiz=T)
		}
	}

	if (is.null(image.color)) {
		#image.color = brewer.pal(7,"Set1")
		#image.color[6] = brewer.pal(7, "PuRd")[3]
		#image.color[1] = brewer.pal(7, "BuPu")[1]
		data('LandscapeColor', package="hustlxc")
		image.color=image.color.lixc
	}
	#par(mar=c(0, 1, 2, 1))
	par(mar=legend.panel.mar)
	plot.empty.graph()
	#legend(x=legend.panel.x, y=legend.panel.y, legend=legend.panel[,2][1:6], fill=rev(image.color), box.col=NA, cex=1.1, border=rev(image.color))
	if (is.null(legend.panel.ncol)) {
		legend("bottomleft", legend=legend.panel[,2][1:6], fill=rev(image.color), box.col=NA, cex=CEXeq1.0, border=rev(image.color))
	}

	if (!is.null(subtype.panel)) {
		subtype.panel.color = subtype.panel.color[1:subtype.categs.num]
		par(mar=c(0, 1, 0, 3))
		plot.empty.graph()
		par(mar=c(0, 1.3, 0., 1.16))
	
		image(as.matrix(subtype.data), col=subtype.panel.color, axes=F)
		grid(dim(middle.panel)[1], 1, col="white", lty=1, lwd=SPACEeq0.06)
		mtext(side=2, text=subtype.panel.identifier, las=2, cex=subtype.panel.identifier.cex)
	
		par(mar=c(0, 0, 0.0, 1))
		plot(0, type='n', xlim=c(0,10), ylim=c(0, 1), axes=F)
		#if (subtype.categs.num == 3) {
		#	text.width=c(0,3.5, 2.5)
		#} else if (subtype.categs.num == 2) {
		#	text.width=c(0,3)
		#}
		# Use text.width to adjust distances between adjacent strings in the legend
		# http://stackoverflow.com/questions/21262472/adjust-spacing-between-text-in-horizontal-legend
		# Use inset to adjust distance between legend and margins, see ?legend
		legend("left", text.width=subtype.legend.text.width, legend=subtype.labels, fill=subtype.panel.color, border=NA, box.col=NA, x.intersp=0.05, ncol=subtype.categs.num, xpd=T, inset=-0.04, cex=subtype.label.cex)
	}
	if (!is.null(subtype.2.panel)) {
		subtype.2.categs <- sort(unique(subtype.2.panel[,2]))
		subtype.2.categs.num <- length(subtype.2.categs)

		subtype.2.panel.color = subtype.2.panel.color[1:subtype.2.categs.num]
		par(mar=c(0, 1, 0, 3))
		plot.empty.graph()
		par(mar=c(0, 1.3, 0., 1.16))
		
		image(as.matrix(subtype.2.data), col=subtype.2.panel.color, axes=F)
		grid(dim(middle.panel)[1], 1, col="white", lty=1, lwd=SPACEeq0.06)
		mtext(side=2, text=subtype.2.panel.identifier, las=2, cex=subtype.2.panel.identifier.cex)
	
		par(mar=c(0, 0, 0.0, 1))
		plot(0, type='n', xlim=c(0,10), ylim=c(0, 1), axes=F)
		#if (subtype.2.categs.num == 3) {
		#	text.width=c(0,3,3)
		#} else if (subtype.2.categs.num == 2) {
		#	text.width=c(0,3)
		#}
		# Use text.width to adjust distances between adjacent strings in the legend
		# http://stackoverflow.com/questions/21262472/adjust-spacing-between-text-in-horizontal-legend
		# Use inset to adjust distance between legend and margins, see ?legend
		legend("left", text.width=subtype.2.legend.text.width, legend=subtype.2.labels, fill=subtype.2.panel.color, border=NA, box.col=NA, x.intersp=0.05, ncol=subtype.2.categs.num, xpd=T, inset=-0.04, cex=subtype.2.label.cex)
	}



	# The left panel
	rownames(left.panel) <- as.character(left.panel[,1])
	left.panel[,1] <- NULL
	left.panel <- left.panel[rev(genes),]
	#par(mar=c(5.5, 1.5, 0, 2.9))
	par(mar=left.panel.mar)
	#par(mgp=c(1, 0.5, 0))
	if (is.null(left.panel.xlim)) {
		left.panel.xlim = c(-1 * max(left.panel[,1]) - 1, 0)
	}
	bp <- barplot(t(-1 * left.panel[,1]), horiz=T, border=NA, col=rev(syn.nonsyn.legend.color), space=SPACEeq0.06, axes=F, xlim=left.panel.xlim, yaxs='i')
	at <- seq(left.panel.xlim[1], 0, length.out=3)
	#axis(side=1, at=c(0, -10, -20, -30), labels=c("0", "10", "20", "30"), cex.axis=1, xpd=TRUE)
	axis(side=1, at=at, labels=as.character(at * -1), cex.axis=1, xpd=TRUE, lwd=0.6, tcl=-0.3, col="grey20")
	#axis(side=1, cex.axis=1, xpd=TRUE)
	axis(side=4, at=bp, labels=left.panel[,2], las=2, tick=F, hadj=0.01, cex.axis=left.panel.cex.axis)
	mtext(text=left.panel.xlab, side=1, padj=1.6, font=1, cex=left.panel.xlab.cex)
	#text(x=-60, y=8, labels="73", col="white")

	# The middle panel
	#par(mar=c(6, 1.3, 0.5, 1.16))
	par(mar=middle.panel.mar)
	image(middle.panel, col=image.color, axes=F)
	grid(dim(middle.panel)[1], dim(middle.panel)[2], col="white", lty=1, lwd=SPACEeq0.06)
	if (!is.null(legend.panel.ncol) && legend.panel.ncol==3) {
		legend(x=0.1, y=legend.panel.y, legend=legend.panel[,2][1:6], fill=rev(image.color), box.col=NA, cex=CEXeq1.1, border=rev(image.color), xpd=T, ncol=3, x.intersp=0.2)
	} 
	else if (!is.null(legend.panel.ncol) && legend.panel.ncol==1) {
		legend(x=0, y=legend.panel.y, legend=legend.panel[,2][1:6], fill=rev(image.color), box.col=NA, cex=0.9, border=rev(image.color), xpd=T, horiz=T, x.intersp=0.2)
	} else if (!is.null(legend.panel.ncol)){
		legend(x=0.2, y=legend.panel.y, legend=legend.panel[,2][1:6], fill=rev(image.color), box.col=NA, cex=CEXeq1.1, border=rev(image.color), xpd=T, ncol=2, x.intersp=0.2)
	}

	# The right panel
	q.value <- right.panel[,2]
	smallest.q.value.clipped=FALSE
	if (any(q.value == 0)) {
		q.value[q.value==0] <- min(q.value) / 2
		smallest.q.value.clipped=TRUE
	} else {
		q.value.n = length(q.value)
		I=which(q.value == min(q.value))
		second.min.q.value=sort(q.value)[2]/100
		if (second.min.q.value/min(q.value) >100) {
			q.value[I]=sort(q.value)[2]/100
			smallest.q.value.clipped=TRUE
		}
	}
	right.panel[,2] <- q.value
	right.panel$qscore <- -log10(right.panel[,2])
	if (smallest.q.value.clipped) {
		warning('q-value was clipped to make visualization better!')
	}

	par(mgp=c(1, 0.5, 0))
	par(mar=right.panel.mar)
	right.panel.color <- right.panel.color
	bp <- barplot(rev(right.panel$qscore), horiz=TRUE, col=right.panel.color, border=NA, space=SPACEeq0.06, axes=F, yaxs='i')
	axis(side=2, at=bp, labels=rev(right.panel[,1]), las=2, font=3, tick=F, hadj=0.9, cex.axis=right.panel.cex.axis, lwd=CEXeq1.0)
	#axis(side=1, at=c(0,4,8), labels=c("0",'4','8'), xpd=TRUE)
	at <- ceiling(seq(0, max(right.panel$qscore), length.out=3))
	axis(side=1, at=at, labels=as.character(at), xpd=TRUE, lwd=0.6, tcl=-0.3, col="grey20")
	mtext(text=right.panel.xlab, side=1, padj=1.6, font=1, cex=right.panel.xlab.cex)
	if (!is.null(right.panel.qscore.cutoff) && right.panel.qscore.lwd!=0) {
		abline(v=right.panel.qscore.cutoff, lty=right.panel.qscore.lty, col=right.panel.qscore.col, lwd=right.panel.qscore.lwd)
	}
	if (!is.null(output.fig.pdffile)) {
		dev.off()
	}
}

