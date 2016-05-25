\name{optimcluster}
\alias{optimcluster}
\title{Optimal clustering via PAM}
\description{
Perform optimal clustering via PAM (package \code{cluster}) and maximum Silhouette index.
}
\usage{
optimcluster(x, D=dist(x), plot=FALSE)
}
\arguments{
  \item{x}{Input data, either a matrix, data frame or vector.}
  \item{D}{Distance among entities in x.}
  \item{plot}{Whether to plot clustering metrics.}
}
\value{
  \item{optim.pam}{Optimal pam object.}
  \item{optim.nclust}{Optimal cluster number.}
  \item{silhouette.index}{Silhouette index}
  \item{Calinski.Harabasz.index}{Calinski-Harabasz(CI) index}
}
\details{

}

\seealso{
\code{pam}
}
\examples{
library(cluster)
library(plotrix)
library(MASS)

x=iris[,1:4]
D=dist(x)
#pdf('clustering_metrics.pdf', height=5, width=6)
r=optimcluster(x, plot=TRUE)
#dev.off()

# Classical multidimensional scaling (CMDS) of a data matrix, or you can choose
#+to use non-metric multidimensional scaling (NMDS) provided in MASS package
mds=cmdscale(D, 2, eig=TRUE)                  # CMDS
#mds=isoMDS(D, k=2, maxit=5000, tol=1e-6)     # NMDS

# visualize MDS output
#pdf('clustering_mds_figure.pdf', height=5, width=6)
s.class2(mds$points, as.factor(r$optim.pam$clustering), grid=FALSE, clabel=0.5, cpoint=0.9, ellipse.lwd=1)
#dev.off()

}

\keyword{utilities}
