\name{plot.depth}
\alias{plot.depth}
\title{Visualize sequencing depth and coverage at different levels}
\description{
Visualizing sequencing depth and coverage at different levels.
}
\usage{

plot.depth <- function(df, 
        left.jitter.amount=0.15, left.bp.notch=T, left.bp.ylim=c(0,300), left.bp.outline=F, left.bp.main="CCDS", left.bp.cex.main=0.8, left.bp.axis.label='CCDS',
        right.jitter.amount=0.3, right.bp.notch=T, right.bp.ylim=c(0.7,1), right.bp.outline=F, right.bp.main="Cov", right.bp.cex.main=0.8, right.bp.axis.labels=NA,
		left.bp.point.col=rgb(0.9, 0.5, 0.5, alpha=0.8),
		right.bp.point.col=rgb(0.9, 0.4, 0.5, alpha=0.8),
        pdffile=NA)
}
\arguments{
  \item{df}{A data frame containing sequencing depth coverage at different levels.}
  \item{left.jitter.amount}{Amount pass to jitter(amount=), controls noise level.}
  \item{left.bp.notch}{Pass to boxplot(notch=).}
  \item{left.bp.ylim}{Pass to boxplot(ylim=).}
  \item{left.bp.outline}{Pass to boxplot(outline=).}
  \item{left.bp.main}{The title.}
  \item{left.bp.cex.main}{Size of cex.main.}
  \item{left.bp.axis.label}{Axis label.}
  \item{right.jitter.amount}{Amount pass to jitter(amount=), controls noise level.}
  \item{right.bp.notch}{Pass to boxplot(notch=).}
  \item{right.bp.ylim}{Pass to boxplot(ylim=).}
  \item{right.bp.outline}{Pass to boxplot(outline=).}
  \item{right.bp.main}{The title.}
  \item{right.bp.cex.main}{Size of cex.main.}
  \item{right.bp.axis.label}{Axis label.}
  \item{left.bp.point.col}{color for the points.}
  \item{right.bp.point.col}{color for the points.}
  \item{pdffile}{If set, save figure to file.}
}
\details{
\code{df} has 8 columns, i.e.\cr
sample, =0X, >=1X, >=4X, >=10X, >=20X, and average_depth\cr
for example:\cr
TB6111  0.017483    0.982517    0.948189    0.900967    0.839802    112.070764\cr
27706A1 0.017566    0.982434    0.951066    0.908361    0.854066    125.642517\cr
TB5539  0.012896    0.987104    0.960002    0.920660    0.871974    129.513491

}
\seealso{
	\code{plot.landscape}, \code{plot.MH}
}
\examples{
data("ccds.depth")
plot.depth(ccds.depth)
}

\keyword{utilities}
