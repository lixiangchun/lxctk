\name{sortDataFrame}
\alias{sortDataFrame}
\title{Sort a data frame by multiple columns or user-specified columns.}
\description{
Sort a data frame by multiple columns or user-specified columns.
}
\usage{
sortDataFrame(data, col.names=NULL, decreasing=FALSE)
}
\arguments{
  \item{data}{A data frame to be sorted.}
  \item{col.names}{User-specified column names to be sorted.}
  \item{decreasing}{Sort in descending or ascending order.}
}
\details{

}
\seealso{
\code{sort.data.frame}
}
\examples{
data('plot.mutation.signature.ex')

sortDataFrame(df)

# Or
sortDataFrame(df, col.names=colnames(df)[1:2,])
}

\keyword{utilities}
