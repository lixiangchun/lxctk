\name{plot.CoocurExclus}
\alias{plot.CoocurExclus}
\title{Visualize mutation coocurrence and exclusivity.}
\description{
Visualize mutaiton coocurrence and exclusivity.
}
\usage{
plot.CoocurExclus(d, genes=NULL, gene.num=26, discard.syn=F, num.of.smp=NULL, image.color=NULL,show.legend=T, legend.panel=NULL,...)
}
\arguments{
  \item{d}{Mutation matrix with gene and sample names as columan and row names, respectively. The integer in each cell represents different mutation types.}
  \item{genes}{A vector of gene names.}
  \item{gene.num}{The number of top recurrently mutated genes selected for visualization.}
  \item{discard.syn}{If TRUE discard syn mutation in calculating gene mutation frequency.}
  \item{num.of.smp}{The number of samples used in calculating gene mutation frequency.}
  \item{image.color}{A vector of 8 elements representing a color scheme.}
  \item{show.legend}{If TRUE draw legend panel}
  \item{legend.panel}{Legend panel, e.g. legend.panel=data.frame(V1=c(7,6,5,4,3,2),V2=c('Nonsense','Splice site','Frame shift','Inframe shift','Missense','Syn.'))}
  \item{...}{Parameters pass to image(...).}
}
\details{

}
\seealso{
	\code{mut_cor_permu_test}
}
\examples{
library(RColorBrewer)
data("CRCMutMat", package='lxctk')
plot.CoocurExclus(mutmat, gene.num=26)
}

\keyword{utilities}
