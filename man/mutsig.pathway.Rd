\name{mutsig.pathway}
\alias{mutsig.pathway}
\title{Canonical algorithms to detect significantly mutated genes}
\description{
Multiple traditional methods to identify significantly mutated pathways based on outputs (*.categs.txt, *.coverage.txt, *.mutations.txt) produced by MutSigCV.
}
\usage{
mutsig.pathway(categs.file, coverage.file, mutations.file, output.file, pathway, sep="", method=c("PCT"), exclude.noncoding=FALSE, only.mutated.gene=FALSE, trace=1, cancer.genes=c())
}
\arguments{
  \item{categs.file}{File with categories mutation rate from MutSigCV, often ended with \code{.categs.txt}.}
  \item{coverage.file}{Coverage file from MutSigCV, often ended with \code{.coverage.txt}.}
  \item{mutations.file}{Mutation file from MutSigCV, often ended with \code{.mutations.txt}.}
  \item{sep}{\code{sep=} for mutations.file, set to \code{\\t} when multiple blank columns present.}
  \item{output.file}{Output file.}
  \item{pathway}{A data frame with 3 cols containg pathway information, the same as the one in \code{mutsigcv.sig.pathway}; A data frame of pathway is used if not provided.}
  \item{method}{Method used to identify SMP, currently only PCT is supported.}
  \item{exclude.noncoding}{If TRUR, exclude \code{noncoding} mutations.}
  \item{only.mutated.gene}{If TRUR, exclude \code{noncoding} mutations in estimating expected mutation number.}
  \item{trace}{If TRUE, print tracing message.}
  \item{cancer.genes}{If provided, only considers genes that is subset of \code{cancer.genes}.}
}
\details{
The poisson convolution test (\code{PCT}) proposed by Kan et al. in \code{Diverse somatic mutation patterns and pathway alterations in human cancers}.
The Fisher's combined P-value test (\code{FCPT}) and likelihood ratio test (\code{LRT}) were 2 of the 3 tests proposed in MuSiC.
The \code{perm.score} is a permutation procedure used to compute p-value for statistic \code{Sg} proposed in page 16 of supplementary of MutSigCV.
The \code{perm.num} permutes mutation number, instead of \code{Sg}, to compute p-value.
\code{PCT}, \code{FCPT} and \code{perm.score} are recommended to use.


}
\seealso{
}
\examples{
categs.file = "TCGA.Breast_Cancer.MutSigCV.categs.txt"
coverage.file = "TCGA.Breast_Cancer.MutSigCV.coverage.txt"
mutations.file = "TCGA.Breast_Cancer.MutSigCV.mutations.txt"
data("CanonicalCancerPathway", package="lxctk")
mutsig.pathway(categs.file, coverage.file, mutations.file, pathway, sep="\t", output.file="out.PCT")
}

\keyword{utilities}
