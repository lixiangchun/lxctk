\name{pileup_freqs}
\alias{pileup_freqs}
\title{Create a table of A,C,G,T nucleotide frequencies for pileup result.}
\description{
Create a table of A,C,G,T nucleotide frequencies for pileup result. Source code from:
https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
}
\usage{
pileup_freqs(plp) 
}
\arguments{
  \item{plp}{Pileup result returned by \code{Rsamtools::pileup(...)}.}
}
\details{
}
\value{
	A data frame (or list, need to check by yourself) of A,C,G,T nucleotide frequencies.
}
\seealso{
\code{calc_allelic_freqs}
}
\examples{

library(Rsamtools)

# 1. Not run:
plp <- Rsamtools::pileup(...)

# 2. Not run -  computing the allelic counts for several SNPs
library(VariantAnnotation)

bamfile <- 'Tumor.rmdup.bam'
bf <- BamFile(bamfile)

#read list of SNPs (several)
vcf <- readVcf("vcfExample_chr17.vcf", "hg19")
vcf.ranges <- rowData(vcf) # You can also create a list of SNPs with GRanges

#compute pileup
#param <- ScanBamParam(which=GRanges("chr17",IRanges(start=685640, end=685650)))
param <- ScanBamParam(which=vcf.ranges)

p_param <- PileupParam(max_depth=1000, ignore_query_Ns=FALSE)
plp <- pileup(bf, scanBamParam=param, pileupParam=p_param)
dim(plp)

#create table of nucleotide frequencies
pileup_freqs(plp)

close(bf)

}

\keyword{utilities}
