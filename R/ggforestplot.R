
# d is a data frame with 4 columns
# d$x gives variable names
# d$y gives center point
# d$ylo gives lower limits
# d$yhi gives upper limits
ggforestplot <- function(d, xlab="Odds Ratio", ylab="Study"){
	if (!is.data.frame(d)) {
		d = as.data.frame(d)
	}
	colnames(d) = c('x','y','ylo','yhi')
	rownames(d) = NULL

	d$x = factor(d$x, levels=d$x) ## Make the order looks the same as in `d`

    p <- ggplot2::ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi)) + 
		geom_pointrange() + 
		coord_flip() +
		geom_hline(aes(x=0, yintercept=1), lty=2, lwd=0.3) +
		ylab(xlab) +
		xlab(ylab) #switch because of the coord_flip() above
    return(p)
}

# Create some dummy data.
#d1 <- data.frame(x = toupper(letters[1:10]), y = rnorm(10, .05, 0.1))
#d1 <- transform(d1, ylo = y-1/10, yhi=y+1/10)
#ggforestplot(d1)

