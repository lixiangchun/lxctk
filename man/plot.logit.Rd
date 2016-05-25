\name{plot.logit}
\alias{plot.logit}
\title{Visualize a logit model object.}
\description{
Visualize a logit model object.
}
\usage{
plot.logit(logitObj=NULL, remove.intercept=TRUE, coef.names=NULL, main=NULL, main.xycoord=c(0.3,0.7),...)
}
\arguments{
  \item{logitObj}{A logit object return by \code{glm(...)}.}
  \item{remove.intercept}{Whether to remove logit intercept in visualization.}
  \item{coef.names}{Variable names to be displayed on \code{forestplot}.}
  \item{main}{The title of the figure.}
  \item{main.xycoord}{Coordinate to place the title.}
  \item{...}{More parameters passed to \code{forestplot}, e.g. \code{clip} etc.}
}
\details{

}
\seealso{
	\code{glm},\code{plot.survfit.lixc},\code{plot.multinom}
}
\examples{
}

\keyword{utilities}
