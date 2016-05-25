\name{stanfit2mcmc.list}
\alias{stanfit2mcmc.list}
\title{Convert rstan::stanfit object to coda:mcmc.list.}
\description{
Convert rstan::stanfit object to coda:mcmc.list.
}
\usage{
stanfit2mcmc.list(stanFit)
}
\arguments{
  \item{stanFit}{A rstan::stanfit object.}
}
\value{
	A coda::mcmc.list.	
}
\details{

}
\seealso{
\code{best.robust_t_test},\code{best.student_t}
}
\examples{
stanFit <- best.robust_t_test(vals, grp, stanDso)
codaObject <- stanfit2mcmc.list(stanFit)
}

\keyword{utilities}
