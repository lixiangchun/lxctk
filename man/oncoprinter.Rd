\name{oncoprinter}
\alias{oncoprinter}
\title{OncoPrinter - Visualizing mutational landscape of SMGs.}
\description{
Visualizing mutational landscape of SMGs.
}
\usage{
oncoprinter(mat,cols,legend.panel)
}
\arguments{
  \item{mat}{A data frame or mutation matrix to be visualized.}
  \item{cols}{A vector of colors used.}
  \item{legend.panel}{A data frame of 2 columns. (1) Mutation tag matched in \code{mat}; (2) Labels of mutation tags.}
}
\details{
See the following example for more info about inputs.
}
\seealso{
\code{...}
}
\examples{
library(lxctk)
data("SMG", package="lxctk")
data("LandscapeColor", package="lxctk")
data("icgc.hcc", package="lxctk") ## ICGC liver cancer data and SMGs selected for NMF clustering.

genes <- names(sort(colSums(icgc.hcc[,Smgs_for_nmf]>2), decreasing=TRUE)) ## Sort by frequency
d <- icgc.hcc[,genes]

#image.color.lixc <- c('white','grey88','#644B39','forestgreen','#FF8B00','#9867CC','#DB1C00', 'black')
## Set the 1st element of image.color.lixc to get rid of background color
image.color.lixc <- c(NA,'grey88','#644B39','forestgreen','#FF8B00','#9867CC','#DB1C00', 'black')

y <- sort.data.frame.by.index(d, rep(0,nrow(d)))
#y <- sort.data.frame.by.index(d, data$nmf_clustid)  ## Sort by black of rows specified in nmf_clustid.

legend.panel <- rbind(data.frame(V1=8,V2='TERT.promoter'), legend.panel)

## If there are too many samples, set the following to get colors shown up.
## Especially when figure is saved in pdf format.
#options(grid.lwd=0)  ## default is 0.001

#options(land.image.legend.ncol=3) ## reset column number of legend, default is 6.

#pdf('oncoprinter.pdf', width=6.5, height=6)
oncoprinter(y, image.color.lixc, legend.panel)
#dev.off()

}

\keyword{utilities}
