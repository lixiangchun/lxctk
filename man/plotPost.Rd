\name{plotPost}
\alias{plotPost}
\title{Plot posterior distribution of a given parameter}
\description{
Plot posterior distribution of a given parameter.
}
\usage{
plotPost( paramSampleVec , cenTend=c("mode","median","mean")[1] , 
	compVal=NULL, ROPE=NULL, credMass=0.95, HDItextPlace=0.7, 
	xlab=NULL , xlim=NULL , yaxt=NULL , ylab=NULL , 
	main=NULL , cex=NULL , cex.lab=NULL ,
	col=NULL , border=NULL , showCurve=FALSE , breaks=NULL , 
	... ) 
}
\arguments{
  \item{paramSampleVec}{A vector of sample associated with a parameter estimated by MCMC.}
}
\details{

}
\seealso{
\code{best.robust_t_test},\code{stan_model},\code{HDIofMCMC}
}
\examples{

## Not Run
stanFit <- best.robust_t_test(...) ## A stanfit object
codaFit <- stanfit2mcmc.list(stanFit) ## convert stanfit object to coda object

plotPost(codaFit[,'mu[1]'],"mean",xlab=expression(mu[tumor]), main='TERT')

plotPost(codaFit[,'mu_diff'][[1]],"mean",xlab=expression(mu[tumor]), main='TERT')

plotPost(codaFit[,'mu_diff'],"mean",xlab=expression(mu[tumor]), main='TERT')
}

\keyword{utilities}
