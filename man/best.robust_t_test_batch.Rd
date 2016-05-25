\name{best.robust_t_test_batch}
\alias{best.robust_t_test_batch}
\title{Batch running robust t-test.}
\description{
}
\usage{
best.robust_t_test_batch(tumor_expression, normal_expression, best.robust_t_test_model_stanDso, outfile, ...)
}
\arguments{
  \item{tumor_expression}{Gene expression matrix where rows are genes and columns are samples.}
  \item{normal_expression}{Gene expression matrix where rows are genes and columns are samples.}
  \item{best.robust_t_test_model_stanDso}{An object of class \code{stanmodel}.}
  \item{outfile}{Output file.}
  \item{...}{Parameters passed to \code{sampling(...)} in package \code{rstan}.}
}
\details{
Running \code{best.robust_t_test} in batch. Note that the row and column names of 
\code{tumor_expression} and \code{normal_expression} must be consistent.
}
\seealso{
\code{best.robust_t_test}
}
\examples{
# Check man page of 'best.robust_t_test' for related info.
}

\keyword{utilities}
