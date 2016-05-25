
oncoprinter.sort.data.frame <- function(data, col.names=NULL, decreasing=TRUE, rev.colnames=TRUE) 
{
	d <- as.data.frame(data)
	if (ncol(d) < 2) {
		stop('A data frame of >=2 columns is required!')
	}
	if (is.null(col.names)) {
		## Method 1:
		col.names <- colnames(data)

		## Do not enable the following clauses since I want to arrange genes in user specified 
		##+order not merely in descending order. If enabled, genes are always display in
		##+descending order.

		## Method 2:
		#gene.mut.freq <- colSums(d>2)
		#gene.mut.freq <- sort(gene.mut.freq, decreasing=TRUE)
		#col.names <- names(gene.mut.freq)
		#d <- d[,col.names] 
	}
	k <- 0
	tmp <- d; tmp[tmp<=k] <- 0; tmp[tmp>k] <- 1
	if (decreasing) {
		#x <- d[do.call("order", -d[col.names]), ]
		x <- d[do.call("order", -tmp[col.names]), ]
	} else {
		#x <- d[do.call("order", d[col.names]), ]
		x <- d[do.call("order", tmp[col.names]), ]
	}
	x
}

## sort a data frame by blocks specified by idx; idx specifies the block of each entry
sort.data.frame.by.index <- function(data, idx, col.names=NULL, decreasing=TRUE, rev.colnames=TRUE)
{
	r <- by(data, idx, oncoprinter.sort.data.frame, col.names, decreasing, rev.colnames)	
	out <- list()
	for (e in names(r)) {
		out[[e]] <- rownames(r[[e]])
	}
	x <- do.call('rbind', r)
	rownames(x) <- do.call('c', out)
	x
}

land.image <- function(mat,col,legend.panel)
{
	if (!is.matrix(mat)) {
		mat <- as.matrix(mat)
	}
	genes <- rev(colnames(mat))
	mat <- mat[,genes]
	gene.mut.freq <- sprintf('%.1f%%', colSums(mat>2) / nrow(mat) * 100 )
	gene.label <- paste(genes, gene.mut.freq)
	gene.label <- paste(gene.label, '  ')
	image(mat,col=col, axes=FALSE, oldstyle=FALSE)
	mtext(gene.label, side=2, at=seq(0,1,length.out=ncol(mat)), las=1, font=3,cex=0.5)
	grid.lwd <- getOption('grid.lwd',0.001)
	if (grid.lwd!=0) { grid(dim(mat)[1], dim(mat)[2], col="white", lty=1, lwd=grid.lwd)}
	legend.ncol <- getOption('land.image.legend.ncol', 6)
	legend(x=0,y=-0.0, legend=legend.panel[,2][1:nrow(legend.panel)], fill=rev(col), box.col=NA, cex=0.7, border=rev(col), ncol=legend.ncol, xpd=TRUE,text.width=0.13,x.intersp=0.2)
}

oncoprinter <- function(mat,col,legend.panel)
{
	land.image(mat,col,legend.panel)
}

subtype.image <- function(mat,col,legend.panel)
{
	if (!is.matrix(mat)) {
		mat <- as.matrix(mat)
	}

	ele.num <- length( unique(stack(as.data.frame(mat))$values) )
	if (ele.num != nrow(legend.panel)) {
		warning('The unique number of values in subtype.mat is not equal to legend number.')
	}

	image(mat,col=col, axes=FALSE, oldstyle=FALSE)
	label <- paste(colnames(mat), '  ')
	mtext(label, side=2, at=seq(0,1,length.out=ncol(mat)), las=1, font=3,cex=0.5)
	grid.lwd <- getOption('grid.lwd',0.001)
	if (grid.lwd!=0) { grid(dim(mat)[1], dim(mat)[2], col="white", lty=1, lwd=grid.lwd)}
	
	## use options(subtype.image.legend.xycoord=c(x, y)) to set `subtype.image.legend.xycoord`
	xycoord <- getOption('subtype.image.legend.xycoord', c(0, 1.6))
	legend.ncol <- getOption('subtype.image.legend.ncol', 6)

	#legend(x=xycoord[1], y=xycoord[2], legend=legend.panel[,2], fill=rev(col), box.col=NA, cex=0.8, border=rev(col), ncol=legend.ncol, xpd=TRUE,text.width=0.13,x.intersp=0.2)
	legend(x=xycoord[1], y=xycoord[2], legend=legend.panel[,2], fill=col, box.col=NA, cex=0.7, border=col, ncol=legend.ncol, xpd=TRUE,text.width=0.13,x.intersp=0.2)
}

freq.by.subtype <- function(gene, subtype, mat)
{
	freq <- by(mat[,gene]>2, subtype, sum) / table(subtype)
	#cat(gene, by(mat[,'TP53']>2, subtype, sum), table(subtype), '\n')
	freq
}
barplot.image <- function(mat, subtype, genes=NULL, ...)
{
	if (is.null(genes)) {
		gene.mut.freq <- sort(colSums(mat>2), decreasing=FALSE)
		genes <- names(gene.mut.freq)
    }
	r <- lapply(genes, freq.by.subtype, subtype, mat)
	y <- do.call('rbind', r)
	#rownames(y) <- genes
	#print(y)
	#y <- y[names(gene.mut.freq),]
	rownames(y) <- NULL
	barplot(t(y), horiz=TRUE, las=2, ...)
}

landplot <- function(d, genes, clin.vars, idx, land.color, land.legend.panel, subtype.color, subtype.legend.panel, barplot.color, barplot.legend.panel)
{
	y <- sort.data.frame.by.index(d, idx, c(genes,clin.vars))
	subtype.mat <- y[,clin.vars]
	if (!is.null(subtype.mat)) {
		landplot.layout.heights <- getOption('landplot.layout.heights', c(40, 80))
		layout(mat=matrix(c(1, 2, 3, 4), nrow=2, byrow=TRUE), width=c(80, 15), heights=c(landplot.layout.heights[1], landplot.layout.heights[2]))
	} else {
		layout(mat=matrix(c(1, 2), nrow=1, byrow=TRUE), width=c(80, 15))
	}
	if (!is.null(clin.vars)) {	
		subtype.mat <- y[, clin.vars]
		par(mar=c(0.02, 6.1, 4.1, 0.1))
		subtype.image(subtype.mat, subtype.color, subtype.legend.panel)
		par(mar=c(0.02, 0.02, 4.1, 1))
		plot(0, type='n', axes=FALSE, xlab="", ylab="")
		par(mar=c(5.1, 6.1, 0.02, 0.1))
	} else {
		par(mar=c(5.1, 6.1, 4.1, 0.1))
	}
	land.mat <- y[,genes]
	land.image(land.mat, land.color, land.legend.panel)
	if (!is.null(clin.vars)) {
		par(mar=c(5.1, 0.02, 0.02, 1))
	}
	else {
		par(mar=c(5.1, 0.02, 4.1, 1))
	}
	barplot.image(d, idx, rev(genes), lwd=0.5, cex.axis=0.6, cex.names=0.6,xpd=TRUE,yaxs='i', las=1, border=NA, col=barplot.color, legend.text=barplot.legend.panel[,2], args.legend=list(ncol=1, box.col=NA,border=NA,text.width=0.13,x.intersp=0.2,xpd=TRUE, x='right', cex=0.7))
	mtext('Mutation frequency', side=1, cex=.7, padj=4)
}

#library(hustlxc)
#library(RColorBrewer)
#data("SMG", package="hustlxc")
#data("LandscapeColor", package="hustlxc")
##data('icgc.hcc.RData', package="hustlxc")
#load('icgc.hcc.RData')

#genes <- names(sort(colSums(icgc.hcc[,Smgs_for_nmf]>2), decreasing=TRUE))
#d <- icgc.hcc[,genes]
#image.color.lixc[1] <- 'white'
#image.color.lixc <- c('white','grey88','#644B39','forestgreen','#FF8B00','#9867CC','#DB1C00', 'black')
#y <- sort.data.frame.by.index(d, rep(0,nrow(d)))
##y <- sort.data.frame.by.index(d, data$nmf_clustid)
#legend.panel <- rbind(data.frame(V1=8,V2='TERT.promoter'), legend.panel)
##pdf('oncoprinter.pdf', width=6.5, height=6)
#oncoprinter(y, image.color.lixc, legend.panel)
##dev.off()
