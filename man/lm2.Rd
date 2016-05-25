\name{lm2}
\alias{lm2}
\title{Perform a simple linear regression model}
\description{
Perform a simple linear regression model by automatically discarding outliers.
}
\usage{
lm2(formula, data, ...)
}
\arguments{
  \item{formula}{an object of class of \code{formula}, see \code{?lm} for detail}
  \item{data}{an optional data frame, see \code{?lm} for detail}
  \item{...}{additional arguments to be passed to \code{lm}, see \code{?lm} for detail}
}
\details{

}
\seealso{
\code{lm},\code{print.infl},\code{influence.measures}
}
\examples{
y <- rnorm(100)
x <- rnorm(100)
d <- data.frame(x=x,y=y)
r <- lm2(y~x, d)
}

\keyword{utilities}
