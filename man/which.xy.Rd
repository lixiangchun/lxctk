\name{which.xy}
\alias{which.xy}
\title{Get indexes of a vector or an element in other unique vector.}
\description{
	Get indexes of a vector or an element in other unique vector.
}
\usage{
which.xy(x, y)
}
\arguments{
  \item{x}{A vector or an element.}
  \item{y}{A vector with no duplicate elements.}
}
\value{
}
\seealso{
\code{which},\code{match},\code{\%in\%}
}
\references{
http://stackoverflow.com/questions/5577727/is-there-an-r-function-for-finding-the-index-of-an-element-in-a-vector
}
\examples{
which.xy(c(2,21,3,10,20), 1:10)
c(2,21,3,10,20) \%in\% 1:10
}
\keyword{}
