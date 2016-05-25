\name{mutsignr}
\alias{mutsignr}
\title{Mutsig analysis of somatic mutations in noncoding region}
\description{
Mutsig analysis of somatic mutations in noncoding region by taking
into account mutatational categories and DNA replication timing.
}
\usage{
mutsiggr(bkgr.mut.file, obs.mut.file, bkgr.mut=NULL,obs.mut=NULL,time=10000, mc.cores=1, outfile='output.sig_elements.txt')
}
\arguments{
  \item{bkgr.mut.file}{Background mutation file stratified by categories and replication timing.}
  \item{obs.mut.file}{Observational mutation file stratified by categories and replication timing.}
  \item{bkgr.mut}{A data frame of \code{bkgr.mut.file}. Set it to \code{NULL} if \code{bkgr.mut.file} is provided.}
  \item{obs.mut}{A data frame of \code{obs.mut.file}. Set it to \code{NULL} if \code{obs.mut.file} is provided.}
  \item{time}{Permutation times.}
  \item{mc.cores}{Number of cores used in mclapply. Set it to 4 or 8 to perform computation in parallel.}
  \item{outfile}{Output file.}
}
\value{
}
\details{

The format of both input files (See \code{Examples}):
The 1st row is the header. The 1st column can be any data type (e.g. string, integer, ...).
The last column is the DNA replication timing cluster ID indicating the S_i stages to 
which this genomic region belongs.

1. Duquenne L, Huvet M, Chen C, et al: Impact of replication timing on non-CpG and CpG substitution rates in mammalian genomes447–457, 2010
2. Koren A, Polak P, Nemesh J, et al: Differential Relationship of DNA Replication Timing to Different Forms of Human Mutation and Variation. Am J Hum Genet 91:1033–1040, 2012
3. Ghosh S, Polak P, Eggan K, et al: Genetic Variation in Human DNA Replication Timing 1–12, 2014

}

\seealso{
\code{mutsigclfn},\code{mutsig.gene},\code{mutsig.pathway}
}
\examples{
library(parallel)
data('cds_start_3k_upstream') ## Data from meta-study of more than 1000 liver cancers
mutsignr(bkgr.mut=bkgr.mut, obs.mut=obs.mut, time=1000, mc.cores=2, outfile='output.sig_elements.txt')

########The format of input file###################
region	*CpG->T	Tp*C->mut	Tp*A->T	C->T	C->A	misc	reptime_cluster
chr1:66091-69091	26	1050	814	954	980	14176	3
chr1:364659-367659	10	1098	546	974	984	14388	0
chr1:622034-625034	16	1014	548	1094	1110	14218	4
chr1:858322-861322	678	660	52	1914	2592	12104	2
chr1:894620-897620	634	660	56	1842	2476	12332	1
###################################################
## NOTE: the column of input file must be the same as the above shown.

}

\keyword{utilities}

