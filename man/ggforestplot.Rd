\name{ggforestplot}
\alias{ggforestplot}
\title{forestplot via ggplot2}
\description{
forestplot via ggplot2
}
\usage{
ggforestplot(d, xlab="Odds Ratio", ylab="Study")
}
\arguments{
  \item{d}{A data frame with 4 columns. 1st column variable names; 2nd column center point; 3rd column lower limits; 4th column upper limits.}
  \item{xlab}{x-label.}
  \item{ylab}{y-label.}
}
\details{

}
\seealso{
\code{rmeta::forestplot}
}
\examples{
# Create some dummy data.
d1 <- data.frame(x = toupper(letters[1:10]), y = rnorm(10, .05, 0.1))
d1 <- transform(d1, ylo = y-1/10, yhi=y+1/10)
ggforestplot(d1)
}

\keyword{utilities}
