\name{plot.survfit.lixc}
\alias{plot.survfit.lixc}
\title{Survival curve visualization with additional information.}
\description{
Survival curve visualization with additional information.
}
\usage{
plot.survfit.lixc(survdiff.formula, data, legend.labels=c('Wild-type','Mutant'), ci.lab.xycoord=c(70,0.2), legend.xycoord=c(40,0.4), col=c("cyan4", "brown"),no.ci.lab=FALSE, pngfig=NULL,...)
}
\arguments{
  \item{survdiff.formula}{An appropriate formula passed to \code{survdiff}, always put the variable-of-interest(VOI) as the last term.}
  \item{data}{A data frame to be operated, column names must match those in the formula.}
  \item{ci.lab.xycoord}{The x and y coordinates for confidence interval.}
  \item{legend.xycoord}{The x and y coordinates for legend.}
  \item{col}{Color for KM survival curves.}
  \item{no.ci.lab}{Do not add confident interval in the plot.}
  \item{pngfig}{A figure suffixed with .png, in which figure will be saved.}
  \item{...}{Parameters passed to \code{plot.survfit}.}
}
\details{

}
\seealso{
	\code{survdiff},\code{survfit},\code{coxph}
}
\examples{
library(survival)
data("crcdat", package='lxctk')
survdiff.formula <- as.formula("Surv(Survival.month, status) ~ hypermutated>0")
#survdiff.formula <- Surv(Survival.month, status) ~ hypermutated>0
plot.survfit.lixc(survdiff.formula, data=crcdat, main="Hypermutated")
}

\keyword{utilities}
