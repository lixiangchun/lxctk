\name{glmfit}
\alias{glmfit}
\title{Plot posterior distribution of a given parameter}
\description{
Plot posterior distribution of a given parameter.
}
\usage{
glmfit(..., clip = c(0, 8), user_specified_variables, plotting_for_sig_only = FALSE,
        main = NULL, main.xycoord = c(0.3, 1.05), use.stats.glm=FALSE)
}
\arguments{
  \item{...}{Parameters passed to \code{stats::glm} or \code{arm::bayesglm}.}
  \item{clip}{Clip used in \code{rmeta::forestplot(...)}.}
  \item{user_specified_variables}{Variables to display on the forest plot.}
  \item{plotting_for_sig_only}{If TRUE, forest plot is generated if there is at least one significant coefficient.}
  \item{main}{The title of the forest plot.}
  \item{main.xycoord}{The coordinate of the title.}
  \item{use.stats.glm}{If TRUE, use \code{stats::glm}. \code{arm::bayesglm} is used by default.}
}
\details{

}
\seealso{
\code{stan_lm},\code{stan_glm}
}
\examples{

## Not Run
gene <- "TP53"
f <- sprintf("%s ~ msig1 + msig2 + msig3 + msig4 + msig5 + msig6 + msig7", gene)
fit = glmfit(
	as.formula(f),
	data = regular,
	family = binomial,
	clip = c(0, 9),
	main = sprintf("%s vs. mutational exposures", gene))

}

\keyword{utilities}
