\name{dist.JSD}
\alias{dist.JSD}
\title{Jensen-Shannon divergence (JSD) distance.}
\description{
Jensen-Shannon divergence (JSD) distance used in gut microbiota clustering, source code
from http://enterotype.embl.de/enterotypes.html
}
\usage{
dist.JSD(inMatrix, pseudocount=0.000001, ...)
}
\arguments{
  \item{inMatrix}{input matrix of sample x taxon}
  \item{pseudocount}{avoid zero in computing}
}
\details{

}
\seealso{
\code{dist}
}
\examples{
}

\keyword{utilities}
