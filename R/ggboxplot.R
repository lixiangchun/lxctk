
## source('../data/ggboxplot_data.RData')

## use with reshape2::melt to reshape data
##ggplot(aes(y = Value, x = Gene, fill = Group), data = data) + geom_boxplot(outlier.shape = NA) + ylab("ylab") + coord_flip() + scale_color_brewer("Set2")


ggboxplot <- function(data, ylab="y-label", show.outlier=FALSE, color_brewer_name="Set2", horizontal=FALSE, ...)
{
	outlier.shape <- NA ## suppress outliers
	#if (show.outlier) {
	#	outlier.shape <- 19 ## default value in geom_boxplot(...)
	#}
	g <- ggplot(aes(y = Value, x = G2, fill = G1), data = data)
	g <- g + geom_boxplot(...) + ylab(ylab)
	if (horizontal) {
		g + coord_flip() + scale_color_brewer(color_brewer_name)
	} else {
		g + scale_color_brewer(color_brewer_name)
	}
}

