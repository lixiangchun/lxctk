

boxplot.jitter <- function(..., horizontal=FALSE, dot.col=1:10, dot.cex=1)
{
	# mpg[2] can be used to control distance between axis label and axis ticks
	#par(mgp=c(3,0.6,0))
	#boxplot(..., horizontal=horizontal, pars=list(outpch=NA, boxwex=0.6))
	#box()
	#axis(side=1, lwd=0.5, las=1)
	#axis(side=2, lwd=0.5, las=1)
	#stripchart(Petal.Length~Species, data=iris, vertical=TRUE, method='jitter', pch=20, col=1:10, add=TRUE, jitter=0.1)
	#stripchart(..., vertical=!horizontal, method='jitter', pch=20, add=TRUE, col=dot.col, cex=dot.cex)

	boxplot(..., horizontal = horizontal, pars = list(outpch=NA,boxwex = 0, staplewex = 0, outwex = 0)) 
    stripchart(..., vertical = !horizontal, method = "jitter", pch = 20, add = TRUE, col = dot.col, cex = dot.cex)
	boxplot(..., horizontal = horizontal, pars = list(outpch = NA,boxwex = 0.6), border=dot.col, add=TRUE)

}
#boxplot.jitter(Petal.Length~Species, data=iris, horizontal=T, las=1, lwd=0.6, cex=0.8)
