\name{plotMutationalProcesses}
\alias{plotMutationalProcesses}
\title{Plot mutation signatures.}
\description{
Plot mutation signatures.
}
\usage{
plotMutationalProcesses(d, color=brewer.pal(6, "Paired"), main='Mutational Signature', ylim.max=0.3, figpdf='mutation_signature.pdf')
}
\arguments{
  \item{d}{Mutational signatures in data.frame.}
  \item{color}{Colors for each mutational signature.}
  \item{main}{Title for mutational signature plot.}
  \item{ylim.max}{Max of ylim.}
  \item{figpdf}{Output fig if provided.}
}
\details{

}
\seealso{
\code{sort.data.frame}
}
\examples{
data('plot.mutation.signature.ex')
plotMutationalProcesses(df, figpdf=NULL)

# Or
plotMutationalProcesses(df)
}

\keyword{utilities}
