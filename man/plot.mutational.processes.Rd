\name{plot.mutational.processes}
\alias{plot.mutational.processes}
\title{Plot mutational processes obtaining from mutational signature deciphering framework (https://github.com/lixiangchun/decipherMutationalSignatures).}
\description{
Plot mutational processes obtaining from mutational signature deciphering framework (https://github.com/lixiangchun/decipherMutationalSignatures).
}
\usage{
plot.mutational.processes(d=NULL, mutational.processes.file = NULL,
  color=brewer.pal(6, "Paired"), main='Mutational Signature', 
  sig.labs=NULL, figpdf='mutation_signature.pdf')
}
\arguments{
  \item{d}{Mutational signatures in data.frame. Check example for data format. Note that the first and second column names must be \code{subtypes} and \code{types}, respectively.}
  \item{mutational.processes.file}{Mutational process file obtained after running \code{decipherMutationalSignatures}. See github \code{https://github.com/lixiangchun/decipherMutationalSignatures} for more information.}
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

## Not run
plot.mutational.processes(mutational.processes.file = system.file('data/Rank_eq_5.processes.txt',package='lxctk'), figpdf=NULL)

plot.mutational.processes(mutational.processes.file = system.file('data/Rank_eq_5.processes.txt',package='lxctk'))

}

\keyword{utilities}
