\name{plot.spectra}
\alias{plot.spectra}
\title{Visualize mutation spectra in 6 or 8 contexts}
\description{
plot mutation spectrum.
}
\usage{
plot.spectra(d, frequency=TRUE, col=brewer.pal(dim(d)[2], "Dark2"), cex.tck=0.7, show.smp.name=TRUE, pdffig=NA, ...)
}
\arguments{
  \item{d}{A matrix with sample as row name and mutation category as columan name.}
  \item{frequency}{Plot the number of mutation instead of fraction.}
  \item{col}{Color scheme for mutation spectrum.}
  \item{cex.tck}{\code{cex} for tick label.}
  \item{show.smp.name}{if FALSE, does not show sample name.}
  \item{pdffig}{file name of figure to be saved.}
  \item{...}{Additional parameters passed to barplot, e.g. main, ylab.}
}
\details{
Set \code{border=NA} to omit bar border entirely.
}
\seealso{
\code{plot.96.spectra}, \code{plot.depth}
}
\examples{
library(RColorBrewer)
library(lxctk)
data("muSpec")
plot.spectra(muSpec8, frequency=TRUE, ylab="Number of mutations", main="Mutation Spectra")
plot.spectra(muSpec8, frequency=FALSE, ylab="Mutation type fraction(\%)", main="Mutation Spectra")
}

\keyword{utilities}
