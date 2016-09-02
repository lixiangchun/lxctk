\name{AlexandrovEtAl_signatures_corr}
\alias{AlexandrovEtAl_signatures_corr}
\title{Calculate cosine similarity of mutational signature with those reported by Alexandrov et al.}
\description{
Calculate cosine similarity of mutational signature with those reported by Alexandrov et al.
}
\usage{
AlexandrovEtAl_signatures_corr(mutational_processes_file, plot=FALSE, outfile, cutoff=0.7, annot_outfile, ...)
}
\arguments{
  \item{mutational_processes_file}{Mutational processes file obtained from \code{decipherMutationalSignatures} pipeline.}
  \item{plot}{If TRUE, visualize similarity with a heatmap.}
  \item{outfile}{If provided, similarity result is exported to \code{outfile}.}
  \item{cutoff}{Cutoff of cosine similarity for downstream annotation.}
  \item{annot_outfile}{If provided, annotation against Alexandrov et al. is written to \code{annot_outfile}.}
  \item{...}{Additional params passed to \code{NMF::aheatmap}.}
}
\details{

}
\seealso{
}
\examples{
}

\keyword{utilities}
