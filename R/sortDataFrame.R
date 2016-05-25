sortDataFrame <- function(data, col.names=NULL, decreasing=FALSE)
{
	orderAscending <- function(...) order(..., decreasing=FALSE)
	orderDescending <- function(...) order(..., decreasing=TRUE)
	if (!is.data.frame(data)) data <- as.data.frame(data)
	if (is.null(col.names)) col.names <- colnames(data)

	if (decreasing) {
		x <- data[ do.call('orderDescending', data[col.names]), ]
	} else {
		x <- data[ do.call('orderAscending', data[col.names]), ]
	}
	x

}
