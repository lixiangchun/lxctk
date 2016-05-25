\name{plot.MH}
\alias{plot.MH}
\title{Visualize mountain-hill landscape with package \code{rgl}}
\description{
Visualizing gene mutation as mountain-hill landscape, containing \code{mountain} genes and \code{hill} genes.
}
\usage{
plot.MH(df, refgene, label.genes=NA, k=400, xmin=4, xmax=300, ymin=4, ymax=268, magnitude=40, shift=2, z.padj=0.5, cex=1)
}
\arguments{
  \item{df}{A data frame with 3 columns, the 1st col is gene name, 2nd mutation frequency, and 3rd q-value.}
  \item{refgene}{A data frame with 3 columns containing information about gene distribution across genome.}
  \item{label.genes}{gene(s) to be labeled in mountain-hill landscape.}
  \item{k}{Number of genes to visualize.}
  \item{xmin}{Start position in x-axis.}
  \item{xmax}{End position in x-axis.}
  \item{ymin}{Start position in y-axis.}
  \item{ymax}{End position in y-axis.}
  \item{magnitude}{Times to be magnified.}
  \item{shift}{Shift units along axis.}
  \item{z.padj}{Adjustment along z-axis.}
  \item{cex=1}{Text size.}
}
\details{
Use package \code{rgl} to visualize gene mutation as mountain-hill landscape.
}
\seealso{
}
\examples{
library(rgl)
data("MH",package="lxctk")
data("refgene",package="lxctk")
plot.MH(MH, refgene, label.genes=c('APC','TP53','PIK3CA'))
}

\keyword{utilities}
