\name{ggPlotMutationalSignatures}
\alias{ggPlotMutationalSignatures}
\title{Plotting mutational signatures}
\description{
Plotting mutational signatures.
}
\usage{
ggPlotMutationalSignatures(mutational.processes.file, d=NULL, color_pallete="Dark2", xlabel="Mutational contexts", ylabel="Percentage of mutations", y.upper.value=NA)
}
\arguments{
  \item{mutational.processes.file}{Mutational process file obtained after running \code{decipherMutationalSignatures}. See github \code{https://github.com/lixiangchun/decipherMutationalSignatures} for more information.}
  \item{d}{A data frame refers to \code{mutational.processes.file}.}
  \item{color_pallete}{The name of color pallete supported by ColorBrewer.}
  \item{xlabel}{X-label.}
  \item{ylabel}{Y-label.}
  \item{y.upper.value}{Truncate y-value to this value if provided (range: 0-100).}
}
\details{
}
\value{
A ggplot object.
}
\seealso{
\code{SomaticSignatures::plotSignatures}
}
\examples{

infile='/Users/lixiangchun/Public/WorkSpace/Project/TJMUCH-Prof.ChenKexin-CRC/lixc/revision/decipherMutationalSignatures/output/Rank_eq_3.processes.txt'
ggPlotMutationalSignatures(infile)

# Or

d=read.table(infile,header=TRUE,stringsAsFactors=FALSE)
x=d[,c(1,2,4,5,3)]
colnames(x)[3:5] = c('Signature1','Signature2', 'Signature3')
ggPlotMutationalSignatures(d=x)

}
\keyword{utilities}
