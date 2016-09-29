\name{BayesSomaticSignatures}
\alias{BayesSomaticSignatures}
\title{Perform Bayesian NMF based mutational signature analysis}
\description{
Perform Bayesian NMF based mutational signature analysis. Source code borrow from Broad Institute CGA.
}
\usage{
BayesSomaticSignatures(x, tumor.type="tumor.type", hyper=FALSE, fafile=NULL, out.dir="OUTPUT_lego96", prior="L1KL", niter=100000)
}
\arguments{
  \item{x}{A data frame with at least 6 named columns, i.e. 'chr', 'start', 'end', 'ref', 'alt' and 'sample'. The seqlevels of mutation data and reference must be identical, or ERROR will occur during mutation context extraction.}
  \item{tumor.type}{A character string specify tumor type analyzed.}
  \item{hyper}{Default = FALSE ; TRUE - to reduce the effect of hyper-mutant samples in the signature discovery.}
  \item{fafile}{BWA-indexed reference sequence. If this is NULL (default), package \code{BSgenome.Hsapiens.UCSC.hg19} must be installed.}
  \item{out.dir}{The output directory.}
  \item{prior}{Default = L1KL (expoential priors); L2KL (half-normal priors).}
  \item{niter}{Number of iterations in Bayesian NMF.}
}
\details{

}
\seealso{
Mutational signature analysis framework in package \code{SomaticSignatures}
}
\examples{

}

\keyword{utilities}
