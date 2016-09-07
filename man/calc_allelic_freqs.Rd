\name{calc_allelic_freqs}
\alias{calc_allelic_freqs}
\title{Calculate allelic frequencies.}
\description{
Calculate allelic frequencies for multiple sites of multiple bam files.
}
\usage{
calc_allelic_freqs(bamFileLists, fastaFile, OncotatorFile, outfile, lines=80, flag=FALSE, minBQ=0, minMQ=13, minDepth=8)
}
\arguments{
  \item{bamFileLists}{A vector of bam files.}
  \item{fastaFile}{Path to indexed reference sequence.}
  \item{OncotatorFile}{Oncotator annotated somatic mutation file.}
  \item{outfile}{Output file.}
  \item{lines}{Number of sites processed at a time, if there are many bam files reduced it to lower memory usage.}
  \item{flag}{If TRUE, a data frame will be returned.}
  \item{minBQ}{minBaseQuality.}
  \item{minMQ}{minMapQuality.}
  \item{minDepth}{allelic frequency is set to 0 if sequencing depth at the mutation site is < minDepth.}
}
\details{
}
\seealso{
\code{calc_allelic_freqs_core}
}
\examples{

library(Rsamtools)
fastaFile="/ifshk1/BC_CANCER/01bin/01bin/lixiangchun/db/hg19Virus/bwa-0.7.12/hg19Virus.fasta"
bamFileLists=read.table("/ifshk5/BC_CANCER/PROJECT/HKC11086_HUMunaX/20150716_intelligent/limiao_paper/recapseg/bam.list",header=FALSE,stringsAsFactors=FALSE)$V1
OncotatorFile="CRC_CapSeg.maf"
calc_allelic_freqs(bamFileLists, fastaFile, OncotatorFile, outfile='CRC_CapSeg_Allelic_freqs.maf', lines=20)
}

\keyword{utilities}
