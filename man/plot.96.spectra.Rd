\name{plot.96.spectra}
\alias{plot.96.spectra}
\title{Visualize mutation spectra in 96 contexts}
\description{
plot mutation spectra at 96 alphabetically sorted contexts.
}
\usage{
plot.96.spectra(x, col=brewer.pal(6, "Paired"), main="spectra.96", max.ylim=NA, draw.legend=TRUE, bottom=FALSE)
}
\arguments{
  \item{x}{The (relative) frequency of mutations at 96 alphabetically sorted contexts.}
  \item{col}{A vector of colors used.}
  \item{main}{The main title.}
  \item{max.ylim}{max ylim.}
  \item{draw.legend}{If TRUE, add legend to current panel.}
  \item{bottom}{Should legend be placed at the bottom?.}
}
\details{

}
\seealso{
}
\examples{
library(RColorBrewer)
data("spectra.96")
x <- spectra.96['67-VS-67TB',]
plot.96.spectra(x)

Multipe plots in a panel, e.g.

pdf("spectra.pdf", width=9, height=7)
par(mfrow=c(6,2), family="sans")
plot.96.spectra(spectra.96['67-VS-67TB',], main="29-VS-29TB")
plot.96.spectra(spectra.96['29-VS-29TB',], main="29-VS-29TB")
plot.96.spectra(spectra.96['3-VS-3TB',], main="3-VS-3TB")
plot.96.spectra(spectra.96['32-VS-29TB',], main="32-VS-29TB")
plot.96.spectra(spectra.96['87-VS-87-TB',], main="87-VS-87-TB")
plot.96.spectra(spectra.96['86-VS-87-TB',], main="86-VS-87-TB")
plot.96.spectra(spectra.96['69-VS-3TB',], main="69-VS-3TB")
plot.96.spectra(spectra.96['88-VS-67TB',], main="88-VS-67TB")
plot.96.spectra(spectra.96['sample_7-3-VS-sample_7-1',], main="sample_7-3-VS-sample_7-1")
plot.96.spectra(spectra.96['sample_7-2-VS-sample_7-1',], main="sample_7-2-VS-sample_7-1")
plot.96.spectra(spectra.96['sample_8-3-VS-sample_8-1',], main="sample_8-3-VS-sample_8-1")
plot.96.spectra(spectra.96['sample_8-2-VS-sample_8-1',], main="sample_8-2-VS-sample_8-1")
dev.off()
}

\keyword{utilities}
