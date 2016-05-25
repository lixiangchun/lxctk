scale01 <- function(x, low = min(x), high = max(x)) {
	x <- (x - low)/(high - low)
	x
}

scale02 <- function(x, scale=c('row', 'column', 'none'), na.rm=TRUE)
{
	if (scale == "row") {
		rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
		rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
	}
	x
}

plotroc <- function(train.grp, predict.prop, rocfig, rocfig.width=5, rocfig.height=5)
{
	pdf(rocfig, width=rocfig.width, height=rocfig.height)
	plot.roc(factor(train.grp), predict.prop, lwd=1, ci=TRUE, identity=TRUE, identity.lwd=0.3, identity.lty=2,print.auc=TRUE, print.auc.x=0.8, grid=1, print.auc.col='blue', las=1, col='red', main="ROC")
	dev.off()
}

CIplot <- function(m, variables, xi, yi, conf.level=0.95, fig=NULL, width=4, height=5)
{
        N <- length(variables)
        MEAN <- rep(0, N)
        LI <- rep(0, N)
        UI <- rep(0, N)
        i = 1
        for (var in variables) {
                x <- m[xi,var]
                y <- m[yi,var]
                r <- t.test(x-y, mu=0, conf.level=conf.level)
                MEAN[i] = mean(x-y)
                LI[i] = r$conf.int[1]
                UI[i] = r$conf.int[2]
        }
        library(plotrix)
        if (!is.null(fig)) { pdf(fig, width=width, height=height) }
        par(mar=c(4,8,2,1))
        plotCI(x=MEAN,y=1:N,li=LI,ui=UI,pch=20,gap=0.02,main="CI plot", las=1, cex=1.2, axes=FALSE, err='x', xlab='Confidence Interval',ylab='')
        abline(v=0, lty=2, lwd=0.4, col='black')
        axis(side=2,at=1:N, variables, las=2, cex.axis=0.4, lwd=0.4)
        axis(side=1, lwd=0.4)
        if (!is.null(fig)) dev.off()
}


runRF <- function(x, y, ..., ccol=brewer.pal(8, 'Set2'), cpoint=2, cellipse=1.5, rocfig.prefix=NULL, rocfig.width=5, rocfig.height=5, mdsfig=NULL, mdsfig.width=5, mdsfig.height=5, indexfig=NULL, indexfig.width=5, indexfig.height=4, impfig=NULL, impfig.width=6, impfig.height=6, CIfig.prefix=NULL, CIfig.width=5, CIfig.height=6, n.imp.var=10, classes=c())
{
	rf = randomForest(x, y, proximity=TRUE, importance=TRUE, ...)
	groups = rownames(rf$confusion)
	if (!is.null(rocfig.prefix)) {
		for (grp in groups) {
			rocfig = sprintf("%s.%s.pdf", rocfig.prefix, grp)
			plotroc(y, rf$votes[,grp], rocfig, rocfig.width, rocfig.height)
		}
	}
	if (!is.null(indexfig)) {
		pdf(indexfig, width=indexfig.width, height=indexfig.height)
		par(lwd=0.5)
		dat=data.frame(index=rf$votes[,2], grp=y)
		#beanplot2(index~grp, dat, las=1, horizontal=TRUE, xlab="Index")
		boxplot.jitter(index~grp, dat, las=1, horizontal=TRUE, xlab="Index", lwd=0.5, cex=0.8, dot.col=ccol)
		dev.off()
	}
	if (!is.null(mdsfig)) {	
		pdf(mdsfig, width=mdsfig.width, height=mdsfig.height)
		par(lwd=0.3)

		#d=dudi.pca(rf$proximity, scale=F, scan=F)
		#s.class(d$li, y, cpoint=cpoint, cellipse=cellipse, grid=FALSE, addaxes=F, col=1:10)
		
		# source code from MDSplot in randomForest package
		rf.mds <- stats:::cmdscale(1 - rf$proximity, eig = TRUE, k = 2)
		s.class(rf.mds$points, y, cpoint=cpoint, cellipse=cellipse, grid=FALSE, addaxes=F, col=ccol)
		

		# unrelated source code
		#d$li['g',] = d$li['g',] - c(0.04,0.014)
		#s.label(d$li+0.012, label=label, add.plot=T)
		#text(d$li+0.01, label=label, xpd=T, col=rep(1:3, each=3), cex=0.6)

		X=c(0.1, 0.3, 0.5, 0.7, 0.9)
		labels=levels(y)
		cols=ccol
		for (i in 1:length(labels)) {
			grid.points(X[i],0.1, default.units='npc', pch=20, gp=gpar(col=cols[i]))
			grid.text(labels[i],X[i]+0.04,0.1, gp=gpar(col=cols[i]))
		}
	
		Y=c(0.95, 0.92,0.89, 0.86,0.83,0.80)
		C=rownames(rf$confusion)
		C.n = length(C)+1
		class.error=rf$confusion[,C.n]
		class.error=format(class.error, digits=2)
		C = c('Class', C)
		class.error=c('Class.error', class.error)

		grid.text(C, x=rep(0.04, C.n), Y[1:C.n], gp=gpar(fontsize=7))
		grid.text(class.error, x=rep(0.13, C.n), Y[1:C.n], gp=gpar(fontsize=7))

		dev.off()
	}
	if (!is.null(impfig)) {
		pdf(impfig, width=impfig.width, height=impfig.height)
		varImpPlot(rf, n.var=min(20, nrow(rf$importance)), cex=0.7, main='Variable importance')
		dev.off()
	}
	if (!is.null(CIfig.prefix) && length(unique(y)) == 2) {
		imp=as.data.frame(importance(rf))
		MeanDecreaseAccuracy=imp$MeanDecreaseAccuracy
		names(MeanDecreaseAccuracy)=rownames(imp)
		MeanDecreaseGini=imp$MeanDecreaseGini
		names(MeanDecreaseGini)=rownames(imp)
		
		imp.var.MeanDecreaseAccuracy=names(sort(MeanDecreaseAccuracy, TRUE)[1:n.imp.var])
		imp.var.MeanDecreaseGini=names(sort(MeanDecreaseGini, TRUE)[1:n.imp.var])
		if (is.null(classes)) { classes =  unique(y)}
		xi=which(y==classes[1])
		yi=which(y==classes[2])
		CIplot(x, imp.var.MeanDecreaseAccuracy, xi, yi, conf.level=0.95, fig=paste(CIfig.prefix, 'MDA.pdf', sep=''), width=impfig.width, height=impfig.height)
		CIplot(x, imp.var.MeanDecreaseGini, xi, yi, conf.level=0.95, fig=paste(CIfig.prefix, 'MDG.pdf', sep=""), width=impfig.width, height=impfig.height)
	}
	invisible(rf)
}



