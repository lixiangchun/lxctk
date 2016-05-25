\name{rank.based.INT}
\alias{rank.based.INT}
\title{Rank-based Inverse Normal Transformation}
\description{
Transfrom input vector with Rank-based Inverse Normal Transformation.
}
\usage{
rank.based.INT(x, c=3/8)
}
\arguments{
  \item{x}{vector of inputs to be transformed.}
  \item{c}{hyper-parameter can be used to control INT: Blom (c=3/8), Tukey (c=1/3), Rankit (c=1/2) and van der Waerden (c=0).}
}
\details{

}
\seealso{
\code{normalize.quantiles} in package \code{preprocessCore}.
}
\examples{
data(SLC12A7_tumor_expression)
y <- rank.based.INT(SLC12A7_tumor_expression)

## for comparison
plot(density(SLC12A7_tumor_expression), xlab='SLC12A7 mRNA level', main='Density plot of raw data')

plot(density(y), xlab='SLC12A7 mRNA level', main='Density plot of INT data')

}

\keyword{utilities}
