
#linreg.yx <- function(y, x,...)
#{
#	data <- data.frame(x=x, y=y)
#	lm.obj = lm(y ~ x, data,...)
	# Detect outliers, which is discarded subsequently.
#	lm.obj.infl <- influence.measures(lm.obj)
#	is.outlier <- apply(lm.obj.infl$is.inf, 1L, any, na.rm = TRUE)  # source code from print.infl(...)
#	if (any(is.outlier)) {
#		x <- x[is.outlier==FALSE]
#		y <- y[is.outlier==FALSE]
#		lm.obj = lm(y ~ x, data,...)
#	}
#	invisible(list(lm.obj=lm.obj, x=x, y=y, is.outlier=is.outlier))
#}
lm2 <- function(formula, data, ...)
{
	lm.obj = lm(formula, data, ...)
	# Detect outliers, which is discarded subsequently.
	lm.obj.infl <- influence.measures(lm.obj)
	is.outlier <- apply(lm.obj.infl$is.inf, 1L, any, na.rm = TRUE)  # source code from print.infl(...)
	if (any(is.outlier)) {
		##data <- data[is.outlier==FALSE,]
		lm.obj <- lm(formula, data[is.outlier==FALSE,], ...)
	}
	invisible(list(lm.obj=lm.obj, data=data, is.outlier=is.outlier))
}
