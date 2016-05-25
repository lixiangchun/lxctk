\name{mut_cor_permu_test}
\alias{mut_cor_permu_test}
\title{MuSiC Mutation Relation Analysis}
\description{
Mutation relation analysis (source code from MuSiC version 0.04).
}
\usage{
mut_cor_permu_test(x=NULL,n.permu=100,seed=NULL,mut.file=NULL,out.file=NULL,out.rdata=NULL)
}
\arguments{
  \item{x}{Mutation relation matrix with gene and sample as column and row names, respectively. Non-zero cell indicates mutation, and zero wild-type.}
  \item{n.permu}{Number of permutations.}
  \item{seed}{Random seed.}
  \item{mut.file}{File to store mutation relation, if x not specified, read from this file.}
  \item{out.file}{Output text file.}
  \item{out.rdata}{Output RDATA file.}
}
\details{

}
\seealso{
}
\examples{
data("CRCMutMat", package='lxctk')
res=mut_cor_permu_test(x=mutmat[1:20, 1:20], n.permu=100)
mut_cor_permu_test(x=mutmat, n.permu=100, out.file="CRCMutMat.txt")
}

\keyword{utilities}
