\name{ggboxplot}
\alias{ggboxplot}
\title{Fancy boxplot via ggplot2}
\description{
Fancy boxplot via ggplot2.
}
\usage{
ggboxplot(data, ylab="y-label", color_brewer_name="Set2", horizontal=FALSE, ...)
}
\arguments{
  \item{data}{See example for data format.}
  \item{...}{Additional params passed \code{geom_boxplot}.}
}
\details{

}
\seealso{
}
\examples{

## source('../data/ggboxplot_data.RData')
data('ggboxplot_data')
#library(reshape2) # if required load it
## use with reshape2::melt to reshape data
ggplot(aes(y = Value, x = Gene, fill = Group), data = data) + geom_boxplot(outlier.shape = NA) + ylab("ylab") + coord_flip() + scale_color_brewer("Set2")

colnames(data) = c('G1','G2','Value')
ggboxplot(data)

ggboxplot(data,horizontal=TRUE,size=0.1,outlier.shape=NA)

}

\keyword{utilities}
