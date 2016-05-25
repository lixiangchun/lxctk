\name{plot.mutational.exposures}
\alias{plot.mutational.exposures}
\title{Visualizing mutational exposures.}
\description{
Visualizing mutational exposures.
}
\usage{
plot.mutational.exposures(d = NULL,
           mutational.exposure.file = NULL,
           clustid = NULL,
           cols =  c('#6C4F4C','#8EA37A','#843818','#AEDD99','#F89867','#FF0000','#A79D90'),
           sorting = FALSE,
           use.colname.as.signature.name = FALSE,
           legend.ncol = 4,
           legend.ypos = 0,
           sorting.by.sample.name = FALSE,
           add.sample.name = FALSE)
}
\arguments{
  \item{d}{A data frame of mutational exposures. Check example for more info about the file format.}
  \item{mutational.exposure.file}{Mutational exposure file outputted from mutational signature deciphering framework \code{https://github.com/lixiangchun/decipherMutationalSignatures}.}
  \item{clustid}{A vector of clustid of each sample.}
  \item{cols}{A vector of colors used.}
  \item{sorting}{If TRUE, sort exposures sequentially.}
  \item{legend.ypos}{The y-position of the legend.}
}
\details{
See the following example for more info about inputs.
}
\seealso{
\code{plot.mutational.processes}
}
\examples{
library(lxctk)
data("icgc.hcc", package="lxctk") ## ICGC liver cancer data and SMGs selected for NMF clustering.

d <- icgc.hcc
mutationalSignatureNames <- paste('mutationalSignature',1:7,sep="")
clustid <- d$nmf_clustid
# clustid <- rep(0, nrow(d))

plot.mutational.exposures(d[,mutationalSignatureNames])

# Or
plot.mutational.exposures(d[,mutationalSignatureNames],clustid)

## Or
plot.mutational.exposures(NULL,system.file('data/Rank_eq_9.exposures.txt',package='lxctk'))

}

\keyword{utilities}
