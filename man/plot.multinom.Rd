\name{plot.multinom}
\alias{plot.multinom}
\title{Visualize a multiple logit model object return by \code{multinom} from package \code{nnet}.}
\description{
Visualize a multiple logit model object return by \code{multinom} from package \code{nnet}.
}
\usage{
plot.multinom(multinomObj, plot=TRUE, remove.intercept=TRUE, coef.names=NULL,main=NULL, main.xycoord=c(0.4, 1),...)
}
\arguments{
  \item{multinomObj}{A multiple logit object return by \code{multinom(...)}.}
  \item{plot}{Whether to plot the figure.}
  \item{remove.intercept}{Whether to remove logit intercept in visualization.}
  \item{coef.names}{Variable names to be displayed on \code{forestplot}.}
  \item{main}{The title of the figure.}
  \item{main.xycoord}{Coordinate to place the title.}
  \item{...}{More parameters passed to \code{forestplot}, e.g. \code{clip} etc.}
}
\details{

}
\seealso{
	\code{glm},\code{plot.survfit.lixc},\code{plot.logit}
}
\examples{
}

\keyword{utilities}
