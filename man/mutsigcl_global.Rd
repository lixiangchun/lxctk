\name{mutsigcl_global}
\alias{mutsigcl_global}
\title{Permutation-based identification of hotspot-enriched genes without taking into account mutational contexts.}
\description{
Permutation-based identification of hotspot-enriched genes without taking into account mutational contexts.
}
\usage{
mutsigcl_global(Oncotator_file, NCBI_CCDS_GeneID2CDS_Positions_Rdata_file, d=NULL, NCBI_CCDS_GeneID2CDS_Positions=NULL, nperm=20000, outfile='mutsigcl_global_output.txt', exclude_indel=TRUE)
generate_mutsigcl_global_bkgr_info(CCDS_current_file)
}
\arguments{
  \item{Oncotator_file}{Oncotator annotated file.}
  \item{NCBI_CCDS_GeneID2CDS_Positions_Rdata_file}{Gene info derived from \code{CCDS.current.txt}, obtained from \code{ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/}.}
  \item{d}{A data frame refers to \code{Oncotator_file}.}
  \item{NCBI_CCDS_GeneID2CDS_Positions}{The R object stored in \code{NCBI_CCDS_GeneID2CDS_Positions_Rdata_file}, obtained from \code{generate_mutsigcl_global_bkgr_info(...)}.}
  \item{nperm}{Number of permutation.}
  \item{outfile}{Output file.}
  \item{exclude_indel}{If TRUE exclude indels.}
}
\value{
}
\details{
}

\seealso{
\code{mutsig.gene},\code{mutsig.pathway},\code{mutsigclfn}
}
\examples{
library(parallel)
library(fastcluster) # accelerate MutSigCL analysis
options(mc.cores=4)  # Use 4 cores in permutation.

# NCBI_CCDS_GeneID2CDS_Positions <- generate_mutsigcl_global_bkgr_info("CCDS.current.txt") ## Generate bkgr info for CCDS genes if required
NCBI_CCDS_GeneID2CDS_Positions_Rdata_file="/ifshk1/BC_CANCER/03user/lixiangchun/iCGA/v0.02/NCBI_CCDS_GeneID2CDS_Positions.RData"
# mutsigcl_global("exomic.maf", NCBI_CCDS_GeneID2CDS_Positions_Rdata_file)

}

\keyword{utilities}
