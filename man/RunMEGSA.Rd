\name{RunMEGSA}
\alias{RunMEGSA}
\title{Running Mutually Exclusive Gene Set Analysis.}
\description{
Running Mutually Exclusive Gene Set Analysis.
}
\usage{
RunMEGSA(mutationMat, maxSSimu=NULL, resultTableFile="resultMEGS.txt", figureDir="figure", nSimu=1000, nPairStart=10, maxSize=6, level =0.05, detail = FALSE, type="pdf", maxSSimuRdataFile='maxSSimu.RData')
funMEGSA(mutationMatFile, maxSSimuFile = NULL, resultTableFile = "resultMEGS.txt", figureDir = "figure", nSimu = 1000, nPairStart = 10, maxSize = 6, level = 0.05, detail = FALSE, type = "pdf")
}
\arguments{
  \item{mutationMat}{Binary (1: mutation, 0: no mutation) mutation matrix. The first row has the gene names and the first column has the patient IDs..}
  \item{maxSSimu}{A data frame loading from simulation file.}
  \item{mutationMatFile}{Binary (1: mutation, 0: no mutation) mutation matrix file.}
  \item{maxSSimuFile}{Simulation file.}
  \item{resultTableFile}{Output file.}
  \item{figureDir}{Output figure directory.}
  \item{nSimu}{Number of simulations (recommand 1000 or more, it may take ~ 10 hours for 1000 simulations).}
  \item{nPairStart}{We first tested all pairs of genes and then pick up the top nPairStart gene pairs (ranked by P-values) to perform multiple-path search to include more genes. Increasing nPairStart will slightly increase power but linearly increasing the computational time. Recommended nPairStart: (10, 30).}
  \item{maxSize}{The maximum size of putative MEGS..}
  \item{level}{Significant level, 0.05 by default.}
  \item{detail}{If TRUE, output detail info during simulation.}
  \item{type}{Figure format, can be either \code{pdf}, \code{png} and \code{jpg}.}
  \item{maxSSimuRdataFile}{The file name that \code{maxSSimu} will be saved.}
}
\details{
See the following example for more info about inputs.
}
\seealso{
\code{...}
}
\source{
1. Xing Hua et al. MEGSA: A Powerful and Flexible Framework for Analyzing Mutual Exclusivity of Tumor Mutations AJHG 2016.
}
\examples{
# If multi-cores available, use
library(parallel)
options(mc.cores=detectCores())
# or
options(mc.cores=min(detectCores(), 8))

RunMEGSA(...)

}

\keyword{utilities}
