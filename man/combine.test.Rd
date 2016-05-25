\name{combine.test}
\alias{combine.test}
\title{Combine multiple p-values}
\description{
Combining multiple p-values. The source code was copied from combine.test(\dots) in package \code{survcomp}.
}
\usage{
combine.test(p, weight, method = c("fisher", "z.transform", "logit"), hetero = FALSE, na.rm = FALSE)
}
\arguments{
  \item{p}{A vector of p-values to be combined.}
  \item{weight}{A vector of weights for corresponding p-values.}
  \item{method}{Method used to combine p-values.}
  \item{hetero}{Always FALSE at this moment.}
  \item{na.rm}{whether to remove na numbers.}
}
\details{
The combine.test is directly copied source code of combine.test(...) implemented in R package "survcomp"!
See http://www.bioconductor.org/packages/release/bioc/html/survcomp.html.
}
\seealso{
}
\examples{
null
}

\keyword{utilities}
