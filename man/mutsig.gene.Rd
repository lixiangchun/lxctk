\name{mutsig.gene}
\alias{mutsig.gene}
\title{Canonical algorithms to detect significantly mutated genes}
\description{
Multiple traditional methods to identify significantly mutated genes based on outputs (*.categs.txt, *.coverage.txt, *.mutations.txt) produced by MutSigCV.
}
\usage{
mutsig.gene(categs.file, coverage.file, mutations.file, mut.df=NULL, output.file, sep="", method=c("PCT","FCPT","LRT","perm.score","perm.num"), exclude.noncoding=FALSE, N=10000, trace=0, only.CGC=FALSE, GOI=NULL, mutations.file.formatted_as_csv=FALSE)
}
\arguments{
  \item{categs.file}{File with categories mutation rate from MutSigCV, often ended with \code{.categs.txt}.}
  \item{coverage.file}{Coverage file from MutSigCV, often ended with \code{.coverage.txt}.}
  \item{mutations.file}{Mutation file from MutSigCV, often ended with \code{.mutations.txt}.}
  \item{mut.df}{A data frame of \code{mutations.file}. \code{mutations.file} is ignored if \code{mut.df} set.}
  \item{sep}{\code{sep=} for mutations.file, set to \code{\\t} when multiple blank columns present.}
  \item{output.file}{Output file.}
  \item{method}{Method used to identify SMG.}
  \item{exclude.noncoding}{If TRUR, exclude \code{noncoding} mutations.}
  \item{N}{Permutation times used in \code{perm.score} or \code{perm.num}.}
  \item{trace}{If TRUE, print tracing message.}
  \item{only.CGC}{If TRUE, only considers cancer gene census.}
  \item{GOI}{A vector of gene-of-interest(GOI), if provided, \code{only.CGC} is set to FALSE.}
  \item{mutations.file.formatted_as_csv}{Formatted \code{mutations.file} in csv.}
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
mutsig.gene(categs.file, coverage.file, mutations.file, sep="\t", output.file="out.PCT")
}

\keyword{utilities}
