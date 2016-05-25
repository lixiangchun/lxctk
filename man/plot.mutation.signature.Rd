\name{plot.mutation.signature}
\alias{plot.mutation.signature}
\title{Plot mutation signatures.}
\description{
Plot mutation signatures. Deprecated, use \code{plot.mutational.processes} instead.
}
\usage{
plot.mutation.signature(d, color=brewer.pal(6, "Paired"), main='Mutational Signature', sig.labs=NULL, figpdf='mutation_signature.pdf')
}
\arguments{
  \item{d}{Mutational signatures in data.frame.}
  \item{color}{Colors for each mutational signature.}
  \item{main}{Title for mutational signature plot.}
  \item{sig.labs}{Labels of mutational signatures.}
  \item{figpdf}{Output fig if provided.}
}
\details{

}
\seealso{
\code{sort.data.frame}
}
\examples{
data('plot.mutation.signature.ex')
plot.mutation.signature(df, figpdf=NULL)

# Or
plot.mutation.signature(df)
}

\keyword{utilities}
