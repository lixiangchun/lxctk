\name{sort.data.frame}
\alias{sort.data.frame}
\title{Sort a data frame by multiple columns or user-specified columns.}
\description{
Sort a data frame by multiple columns or user-specified columns.
}
\usage{
sort.data.frame(df, col.names=NULL)
}
\arguments{
  \item{df}{A data frame to be sorted.}
  \item{col.names}{User-specified column names to be sorted.}
}
\details{

}
\seealso{
}
\examples{
data('plot.mutation.signature.ex')

sort.data.frame(df)

# Or
sort.data.frame(df, col.names=colnames(df)[1:2,])
}

\keyword{utilities}
