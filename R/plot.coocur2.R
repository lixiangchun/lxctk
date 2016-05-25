

plot.coocur2 <- function(d, genes=NULL, sort.by.freq=TRUE, k=2, gene.num=26, right.panel=NULL,right.panel.col=NULL, dist.to.left=5, discard.syn=F, num.of.smp=NULL, image.color=NULL, show.legend=T, legend.panel=NULL,subtype=NULL,subtype.color=NULL, subtype.labels=NULL,subtype.legend=NULL,...)
{
if (is.null(num.of.smp)) {
	num.of.smp=nrow(d)
}
if (discard.syn) {
	gene.freq=colSums(d>2)/num.of.smp
} else {
	gene.freq=colSums(d!=0)/num.of.smp
}

if (is.null(image.color)) {
	image.color = brewer.pal(7, "Set1")
	image.color[6] = brewer.pal(7, "PuRd")[3]
	image.color[1] = brewer.pal(7, "BuPu")[1]
}

if (is.null(genes)) {
	genes=names(sort(gene.freq, decreasing=TRUE)[1:gene.num])
	gene.num=gene.num
} else {
	genes=intersect(genes,colnames(d))
	gene.num=length(genes)
}

if (!is.null(genes) && sort.by.freq) {
	##genes=names(sort(colSums(d[,genes] > 2), decreasing=TRUE))
	genes=names(sort(colSums(d[,genes] > k), decreasing=TRUE))
}

mat=as.matrix(d[,genes])

# Version 1:
#mat=mat[do.call(order, as.data.frame(-mat)),]

# Version 2:
tmp <- mat
tmp[tmp>0] <- 1
mat <- mat[do.call(order, as.data.frame(-tmp)),]

mat=mat[,rev(genes)]
freq=paste(format(gene.freq[genes] * 100, digits=3), '%', sep="")
labels=paste(genes, freq, sep=" ")
labels=paste(labels, "", sep="  ")

if (show.legend) {
	if (is.null(subtype)) {
		layout(matrix(c(1,2,1,2), nrow=2, byrow=TRUE), widths=c(80,30))
		par(mar=c(2,dist.to.left,3,0))
	} else {
		layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE), widths=c(80,30), heights=c(12, 90))
		par(mar=c(2,dist.to.left,3,0))
	}
} else {
	par(mar=c(1,dist.to.left,3,1))
}
neither_matrix_data.frame=FALSE
if (!is.null(subtype)) { 
	if (is.matrix(subtype) || is.data.frame(subtype)) {
		subtype = as.matrix(subtype[rownames(mat),]) 
	} else {
		subtype = as.matrix(subtype[rownames(mat)]) 
		neither_matrix_data.frame=TRUE
	}
}
if (is.null(subtype.color)) {subtype.color=as.integer(unique(subtype))}
if (!is.null(subtype)) {
	if (neither_matrix_data.frame) {
		par(mar=c(0.1,dist.to.left,4,0))
	} else {
		par(mar=c(0.1,dist.to.left,ncol(subtype)-1,0))
	}
	if (is.null(subtype.labels)) {
		subtype.labels = unique(subtype)
	}
	image(subtype, col=subtype.color, axes=F)
	#legend(x=0.1, y=6, legend = subtype.labels, fill = rev(subtype.color), box.col = NA, cex=0.9,border = rev(subtype.color), xpd=TRUE, horiz=T)
	#mtext('  Subtype  ', side=4, las=2, cex=0.8, xpd=TRUE, at=0, font=2)
	mtext(subtype.labels, side=2, las=2, cex=0.6, xpd=TRUE, at=seq(0,1,length.out=3), font=2)

	par(mar=c(0,0,2,0))
	plot(0, type='n', axes=FALSE, xlab="", ylab="")
	legend('topleft', legend = subtype.legend[,2], fill=subtype.color, box.col = NA, cex=0.6,border = subtype.color, x.intersp=0.2, text.width=0.15, ncol=3, pt.cex=0.6)
	#if (is.null(subtype.labels)) {
	#	subtype.labels = unique(subtype)
	#}
	#legend('left', legend = subtype.labels, fill = rev(subtype.color), box.col = NA, cex=0.9,border = rev(subtype.color), horiz=T)
}

image.color=image.color
if (is.null(subtype)) {
	par(mar=c(2,dist.to.left,3,0))
} else {par(mar=c(2,dist.to.left,0,0))}

image(mat, col=image.color, axes=FALSE,...)
mtext(rev(labels), side=2, las=2, cex=0.5, xpd=TRUE, at=seq(0,1,length.out=gene.num), font=3)
if (is.null(legend.panel)) {
	legend.panel=data.frame(V1=c(7,6,5,4,3,2),V2=c('Nonsense','Splice site','Frame shift','Inframe shift','Missense','Syn.'))
}

if (show.legend) {
	if (is.null(right.panel)) {
		par(mar=c(2,0,3,0))
		plot(0, type='n', axes=FALSE, xlab="", ylab="")
		legend('left', legend = legend.panel[,2], fill = rev(image.color), box.col = NA, cex=0.6,border = rev(image.color), x.intersp=0.2, ncol=2)
	} else {
		legend(0, 0, legend = legend.panel[,2], fill = rev(image.color), box.col = NA, cex=0.6,border = rev(image.color), x.intersp=0.2, ncol=3, xpd=TRUE)
		par(mar=c(2,0,0,2))
		#if (all(rownames(right.panel) == rev(genes))) {
		#	stop('Right panel row names are not identical to gene symbols')
		#}	
		colnames(right.panel)=NULL
		r=barplot(right.panel, horiz=TRUE, yaxs='i', space=0.05, border=NA, axes=FALSE, col=right.panel.col)
		
		val = seq(0, max(colSums(right.panel)), length.out=3)
		val = format(val*100, digits=3)
		labels = paste(val, "%", sep="")
		axis(side=1, lwd=0.6, tcl=-0.3, at=as.numeric(val)/100, labels=labels)
	}
}

}

