\name{sort.data.frame.by.index}
\alias{sort.data.frame.by.index}
\title{Sort data frame by columns allowing for different block of rows.}
\description{
Sort data frame by columns allowing for different block of rows.
}
\usage{
sort.data.frame.by.index(data, idx, col.names=NULL, decreasing=TRUE, rev.colnames=TRUE)
}
\arguments{
  \item{data}{A data frame to be sorted, e.g. rows are samples and columns are genes.}
  \item{idx}{Block indicator, e.g. a vector of cluster IDs of samples in \code{data}.}
  \item{col.names}{Columns specified by names to be sorted in \code{data}. If absent, sort \code{data} sequentially from the 1st column to the last column.}
  \item{decreasing}{If TRUE sort \code{data} in decreasing order.}
  \item{rev.colnames}{If TRUE reverse column names before sorting \code{data}.}
}
\details{
See the example of \code{oncoprinter} for more info.
}
\seealso{
\code{oncoprinter}
}
\examples{

}

\keyword{utilities}
