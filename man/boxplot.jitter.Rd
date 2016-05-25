\name{boxplot.jitter}
\alias{boxplot.jitter}
\title{boxplot with scatter points}
\description{
boxplot with scatter points
}
\usage{
boxplot.jitter(..., horizontal=FALSE, dot.col=1:10, dot.cex=1)
}
\arguments{
  \item{...}{Paramters passed to \code{boxplot} and/or \code{stripchart}}
  \item{horizontal}{If TRUE, boxplot with horizontal orientation}
  \item{dot.col}{colors for scatter points}
}
\details{

}
\seealso{
\code{boxplot}, \code{stripchart}
}
\examples{
# mpg[2] is used to control distance between labels and ticks.
par(mar=c(5,5,4,2), mgp=c(3,0.6,0))
boxplot.jitter(Petal.Length~Species, data=iris, horizontal=T, las=1, lwd=0.6, cex=0.8)
}

\keyword{utilities}
