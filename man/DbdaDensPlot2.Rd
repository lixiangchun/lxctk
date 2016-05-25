\name{DbdaDensPlot2}
\alias{DbdaDensPlot2}
\title{Density plot of posterior distribution of a given model parameter.}
\description{
Density plot of posterior distribution of a given model parameter.
}

\usage{
DbdaDensPlot2( codaObject , parName=varnames(codaObject)[1] , plColors=NULL, 
xlab='Param. Value', main=paste('Posterior distribution of\n', xlab) )
}
\arguments{
  \item{codaObject}{A coda::mcmc.list object.}
  \item{parName}{The parameter name.}
  \item{plColors}{Colors to be used for plotting.}
  \item{xlab}{xlabel.}
  \item{main}{Title of the figure.}
}
\details{

}
\seealso{
\code{best.robust_t_test},\code{best.student_t}
}
\examples{
library(lxctk)
library(coda)
library(rstan)

x <- rt(1000, df=1)
fit <- best.student_t(x, stanDso)
codaObject <- stanfit2mcmc.list(fit)
#diag.mcmc(codaObject)
DbdaDensPlot2(codaObject)


data('x')
##fit <- best.student_t(x, stanDso) ## The above error will be generated
fit <- best.student_t(x, stanDso, unifLo=0, unifHi=.Machine$double.xmax) ## No errors
codaObject <- stanfit2mcmc.list(fit)
#diag.mcmc(codaObject)
DbdaDensPlot2(codaObject)

}

\keyword{utilities}
