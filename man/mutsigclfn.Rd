\name{mutsigclfn}
\alias{mutsigclfn}
\title{Permutation-based identification of Significantly Mutated Genes, i.e. MutSigCL or MutSigFN}
\description{
Permutation-based identification of Significantly Mutated Genes, i.e. MutSigCL or MutSigFN.
}
\usage{
mutsigclfn(bkgrSQLiteDB, obs_data, outfile='out.txt', genes=c(), type=c('CL','FN'), hotspot.alg=c('hclust','ratio'), min.cl=0.2, nperm=1000, mc.cores=4, bkgr_data=dbConnect(dbDriver('SQLite'), bkgrSQLiteDB))
}
\arguments{
  \item{bkgrSQLiteDB}{An SQLite DB storing background mutation-related info.}
  \item{obs_data}{A data frame of the observed input data.}
  \item{outfile}{Output file for mutsigclfn.}
  \item{genes}{Genes to be permutated.}
  \item{type}{Algorithm, i.e. CL - MutSigCL, FN - MutSigFN.}
  \item{hotspot.alg}{Algorithm to define hotspot statistic.}
  \item{min.cl}{Genes with fraction of hotspot >min.cl are selected for MutSigCL analysis. To disable it set min.cl to 0.}
  \item{nperm}{Number of permutation.}
  \item{mc.cores}{Number of cores used in mclapply.}
  \item{bkgr_data}{An RSQLite object, never change unless your're quite sure!}
}
\value{
}
\details{
Make sure that \code{bkgrSQLiteDB} and \code{obs_data} are consistent.

When hotspot.alg is set to 'hclust', mutsigclfn will employ the following algorithm to define hotspot statistic (fraction):

"A hotspot is defined as a 3-base-pair region of the gene containing many mutations: at least 2, and at least 2% of 
the total mutations (nature12912)."

This involves using \code{hclust} in package \code{stats} or \code{fastcluster} to perform hierarchical clustering and
call \code{cutree} to calculate hotspot statistic. It's quite time-consuming in large permutations.

When hotspot.alg is set to 'ratio' (much faster than 'hclust' above), I define the hotspot statistic as:

The ratio of the number of mutations to the total number of mutation positions, at least ≥2 mutations is required.

Change Log:
Function name \code{pbiSMG} was changed to \code{mutsigclfn} in v 1.45.

}

\seealso{
\code{mutsig.gene},\code{mutsig.pathway}
}
\examples{
library(parallel)
library(RSQLite)
library(fastcluster) # accelerate MutSigCL analysis
mutsigclfn(...)
}

\keyword{utilities}
