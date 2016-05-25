\name{diag.mcmc}
\alias{diag.mcmc}
\title{MCMC diagnostics plot.}
\description{
MCMC diagnostics plot.
}

\usage{
diag.mcmc( codaObject, parName=varnames(codaObject)[1], 
          DBDAplColors = c("skyblue","black","royalblue","steelblue"),
          figName=NULL)
##diagStanFit(stanFit,parName,saveName=NULL,saveType="jpg") ## Does not work yet
}
\arguments{
  \item{codaObject}{A coda::mcmc.list object or rstan::stanfit object.}
  \item{stanFit}{A rstan::stanfit object.}
  \item{parName}{The parameter name.}
  \item{DBDAplColors}{Colors to be used for plotting.}
  \item{figName}{A output figure name, must be ended with \code{.png}.}
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
diag.mcmc(codaObject)


data('x')
##fit <- best.student_t(x, stanDso) ## The above error will be generated
fit <- best.student_t(x, stanDso, unifLo=0, unifHi=.Machine$double.xmax) ## No errors
codaObject <- stanfit2mcmc.list(fit)
diag.mcmc(codaObject)

}

\keyword{utilities}
