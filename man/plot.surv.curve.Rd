\name{plot.surv.curve}
\alias{plot.surv.curve}
\title{Survival curve visualization with additional information.}
\description{
Survival curve visualization with additional information.
}
\usage{
plot.surv.curve(gene, survdiff.formula, coxph.formula, data, ci.lab.xycoord=c(70,0.2), legend.xycoord=c(40,0.4), pngfig=NULL, main=gene, line.lwd=par("lwd"))
}
\arguments{
  \item{gene}{Variable of interest, e.g. a mutated gene.}
  \item{survdiff.formula}{An appropriate formula passed to \code{survdiff}, always put the variable-of-interest(VOI) as the last term.}
  \item{coxph.formula}{An appropriate formula for CoxPH regression analysis, alias put the VOI as the last term.}
  \item{data}{A data frame to be operated, column names must match those in the formula.}
  \item{ci.lab.xycoord}{The x and y coordinates for confidence interval.}
  \item{legend.xycoord}{The x and y coordinates for legend.}
  \item{pngfig}{A figure suffixed with .png, in which figure will be saved.}
  \item{main}{The main title for the figure.}
  \item{line.lwd}{Line width of survival curve.}
}
\details{

}
\seealso{
	\code{survdiff},\code{survfit},\code{coxph}
}
\examples{
library(survival)
data("crcdat", package='lxctk')
survdiff.formula = as.formula("Surv(Survival.month, status) ~ hypermutated>0")
plot.surv.curve('hypermutated', survdiff.formula, survdiff.formula, data=crcdat, main="")
}

\keyword{utilities}
