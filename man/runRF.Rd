\name{runRF}
\alias{runRF}
\title{Run random forest classifier training and visualization}
\description{
Run random forest classifier training and visualization.
}
\usage{
runRF(x, y,                                 # mandatory data passed to randomForest, e.g. randomForest(x=x, y=y)
	...,                                    # parameters for randomForest except proximity because it's set internally.
	ccol=brewer.pal(8, 'Set2')           # colours for MDS figure
	cpoint=2, cellipse=1.5,                 # pass to s.class in package ade4
	rocfig.prefix=NULL,                     # prefix for ROC plot, if NULL does not plot
	rocfig.width=5, rocfig.height=5,
	mdsfig=NULL,                            # figure name for MDS plot, if NULL does not plot
	mdsfig.width=5, mdsfig.height=5,
	indexfig=NULL,                          # figure name for index plot, if NULL does not plot
	indexfig.width=5, indexfig.height=5,
	impfig=NULL, impfig.width=6, impfig.height=6)  # Importance of variables
}
\arguments{
  \item{x}{Data passed to randomForest(x=x)}
  \item{y}{A vector of class factor passed to randomForest(y=y)}
  \item{...}{Other parameters passed to randomForest(...)}
}
\details{

}
\seealso{
\code{randomForest}
}
\examples{
library(randomForest)
library(ade4)
library(grid)
library(pROC)
data('runRF')
runRF(x=x, y=y, ntree=500, do.trace=TRUE, rocfig.prefix="roc", mdsfig='mds.pdf', indexfig='index.pdf')
}

\keyword{utilities}
