\name{vioplot}
\alias{vioplot}
\title{violin plot}
\description{
  Produce violin plot(s) of the given (grouped) values.
}
\usage{
vioplot( x, ..., range=1.5, h, ylim, names, horizontal=FALSE, 
  col="magenta", border="black", lty=1, lwd=1, rectCol="black", 
  colMed="white", pchMed=19, at, add=FALSE, wex=1, 
  drawRect=TRUE)
}
\arguments{
  \item{x}{data vector}
  \item{...}{additional data vectors} 
  \item{range}{a factor to calculate the upper/lower adjacent values}
  \item{h}{the height for the density estimator, if omit as explained in sm.density, h will be set to an optimum}
  \item{ylim}{y limits}
  \item{names}{one label, or a vector of labels for the datas must match the number of datas given}
  \item{col, border, lty, lwd}{Graphical parameters for the violin passed to lines and polygon}
  \item{rectCol, colMed, pchMed}{Graphical parameters to control the look of the box}
  \item{drawRect}{logical. the box is drawn if \code{TRUE}.}
  \item{at}{position of each violin. Default to \code{1:n}}
  \item{add}{logical. if FALSE (default) a new plot is created}
  \item{wex}{relative expansion of the violin. }
  \item{horizontal}{logical. horizontal or vertical violins}
}
\details{
 A violin plot is a combination of a box plot and a kernel density plot. 
 Specifically, it starts with a box plot. It then adds a rotated kernel density plot to each side of the box plot.
}
\examples{
  # box- vs violin-plot 
  par(mfrow=c(2,1))
  mu<-2
  si<-0.6
  bimodal<-c(rnorm(1000,-mu,si),rnorm(1000,mu,si)) 
  uniform<-runif(2000,-4,4)
  normal<-rnorm(2000,0,3)
  vioplot(bimodal,uniform,normal)
  boxplot(bimodal,uniform,normal)
  
  # add to an existing plot
  x <- rnorm(100)
  y <- rnorm(100)
  plot(x, y, xlim=c(-5,5), ylim=c(-5,5))
  vioplot(x, col="tomato", horizontal=TRUE, at=-4, add=TRUE,lty=2, rectCol="gray")
  vioplot(y, col="cyan", horizontal=FALSE, at=-4, add=TRUE,lty=2)
}
\seealso{
  \code{\link{boxplot}}
  \code{\link[sm]{sm}}
}
\references{
  Hintze, J. L. and R. D. Nelson (1998).  \emph{Violin plots: a box plot-density trace
synergism.}  The American Statistician, 52(2):181-4.
}
\author{Daniel Adler \email{dadler@uni-goettingen.de}
  
  Romain Francois \email{francoisromain@free.fr} ; \url{http://francoisromain.free.fr} : 
  horizontal violins and additionnal graphical parameters }
\keyword{hplot}

