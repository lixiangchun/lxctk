
rank.based.INT <- function(x, c=3/8)
{
	r <- rank(x)
	N <- length(x)
	qnorm((r-c)/(N-2*c+1))
}
