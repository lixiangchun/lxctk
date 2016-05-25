

plot.coocur <- function(d, genes=NULL, gene.num=26, dist.to.left=5, discard.syn=F, num.of.smp=NULL, image.color=NULL, show.legend=T, legend.panel=NULL,subtype=NULL,subtype.color=NULL, subtype.labels=NULL,...)
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
	genes=names(sort(gene.freq, decreasing=T)[1:gene.num])
	gene.num=gene.num
} else {
	genes=genes
	gene.num=length(genes)
}
mat=as.matrix(d[,genes])
mat=mat[do.call(order, as.data.frame(-mat)),]
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
	#par(mar=c(1,dist.to.left,3,1), family="Arial")
}

if (!is.null(subtype)) { subtype = as.matrix(subtype[rownames(mat)]) }
if (is.null(subtype.color)) {subtype.color=as.integer(unique(subtype))}
if (!is.null(subtype)) {
	par(mar=c(0.1,dist.to.left,4,0))
	if (is.null(subtype.labels)) {
		subtype.labels = unique(subtype)
	}
	image(subtype, col=subtype.color, axes=F)
	legend(x=0.1, y=6, legend = subtype.labels, fill = rev(subtype.color), box.col = NA, cex=0.9,border = rev(subtype.color), xpd=TRUE, horiz=T)
	mtext('  Subtype  ', side=4, las=2, cex=0.8, xpd=TRUE, at=0, font=2)

	par(mar=c(0,0,2,0))
	plot(0, type='n', axes=FALSE, xlab="", ylab="")
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
	par(mar=c(2,0,3,0))
	plot(0, type='n', axes=FALSE, xlab="", ylab="")
	#legend('left', legend = legend.panel[,2][1:5], fill = rev(image.color), box.col = NA, cex=0.9,border = rev(image.color))
	legend('left', legend = legend.panel[,2], fill = rev(image.color), box.col = NA, cex=0.9,border = rev(image.color))
}

}

