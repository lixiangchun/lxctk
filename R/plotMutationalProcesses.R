
sort.data.frame <- function(df, col.names=NULL)
{
	if (!is.null(col.names)) {  # sort by specific column(s)
		sorted.df=df[do.call("order", df[col.names]), ]
	} else {  # if col.names not supplied, sort each col one-by-one
		sorted.df=df[do.call("order", df[colnames(df)]), ]
	}
sorted.df
}

plotMutationalProcesses <- function(d, color=brewer.pal(6, "Paired"), main='Mutational Signature', ylim.max=0.3, figpdf='mutation_signature.pdf')
{

d$types=as.factor(d$types)
d$subtypes=as.factor(d$subtypes)
if (!is.data.frame(d)) {
	d <- as.data.frame(d)
}
NR=ncol(d)
if (length(color) != 6) {
	warning('Length of color vector is not equal to 6, use default colour scheme!\n')
	color = brewer.pal(6, "Paired")
}

# Sort by 'types', followed by subtypes
#d=d[do.call("order", d[c('types','subtypes')]), ]
d=sort.data.frame(d, c('types','subtypes'))

if (!is.null(figpdf)) pdf(figpdf, width=6, height=NR-0.55)

col0=color
col=rep(col0, each=16)

x=rep(18.7, 6)
par(mfrow=c(NR,1), family="sans")

par(mar=c(0,4,6,2))
r=barplot(x, xaxs='i', axes=FALSE, border=NA, space=0.02, col=col0, main=main, cex.main=1.5)
mtext(levels(d[,2]), side=3, at=as.numeric(r), cex=0.8, xpd=TRUE)
#mtext(levels(d[,2]), side=3, at=as.numeric(r), cex=0.8, xpd=TRUE, col=color)

axis.at <- seq(0, ylim.max, length.out=3)
axis.at.labels <- sprintf('%.1f',axis.at * 100)
signatureNames <- colnames(d)
YLIM=max(axis.at)

for (i in 3:NR) {
	par(mar=c(1/3,4,0,2))
	barplot(d[,i], col=col, border=NA, xaxs='i', las=2, tcl=-0.2, cex.axis=0.8, ylim=c(0,YLIM), axes=FALSE)
	signatureName <- signatureNames[i]
	text(24, 0.11, signatureName, cex=1.5)
	#axis(side=2, at=c(0,0.04,0.08,0.14), labels=c("0%","4%","8%","14%"), las=2, lwd=0.5)
	axis(side=2, at=axis.at, labels=axis.at.labels, las=2, lwd=0.5)
	par(mar=c(1/3,4,0,2))
}

par(mar=c(6,4,0,2))
r=barplot(x, xaxs='i', axes=FALSE, border=NA, space=0.02, col=col0)
#axis(side=1, at=as.numeric(r), labels=levels(d$V2), xaxt='')
mtext(levels(d[,2]), side=1, at=as.numeric(r), cex=0.8)
#mtext(levels(d[,2]), side=1, at=as.numeric(r), cex=0.8, col=color)

if (!is.null(figpdf)) {
	dev.off()
}

invisible(d)
}
