\name{cmds.focal.test}
\alias{cmds.focal.test}
\title{Correlation matrix diagonal segmentation(CMDS) to identify recurrent CNV chromosome by chromosome.}
\description{
Correlation matrix diagonal segmentation(CMDS) to identify recurrent CNV.
}
\usage{
cmds.focal.test(data, out.file, wsize=30, wstep=1, chromosomes=paste("chr",1:22, sep=""))
}
\arguments{
  \item{data}{A data frame with \code{chromosome}, \code{position} and CN in each sample.}
  \item{out.file}{Output file.}
  \item{wsize}{Block size.}
  \item{wstep}{Step size.}
  \item{chromosomes}{A vector of chromosomes analyzed.}
}
\details{

Qunyuan Zhang et al. CMDS: a population-based method for identifying recurrent DNA copy number aberrations in cancer from high-resolution data Bioinformatics 2010(464-469).


The R source code was downloaded from:
https://dsgweb.wustl.edu/qunyuan/software/cmds

}
\seealso{
null
}
\examples{
data("cmds")
cmds.focal.test(chr11q, out.file="out.txt", chromosomes="11")
}

\keyword{utilities}
