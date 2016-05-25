\name{bayes.student_t}
\alias{bayes.student_t}
\title{Bayesian inference to estimate the parameters of student's t distribution.}
\description{
Bayesian inference to estimate the parameters of student's t distribution.

An example of student's t distirbution:

y ~ student_t( nu , mu , sigma )

\code{nu}, \code{mu} and \code{sigma} will be estimated via bayesian inference.

Deprecated and will be removed in future release, so use \code{best.student_t} instead.

}
\usage{
bayes.student_t(vals, stanDso, ...)
}
\arguments{
  \item{vals}{A vector of metric values assumed to be generated from student's t distirbution.}
  \item{stanDso}{An object of class \code{stanmodel}.}
  \item{...}{Parameters passed to \code{sampling(...)} in package \code{rstan}.}
}
\details{
Estimate student's t distribution via Bayesian statistics provided by \code{rstan} package.

Please refer to:
Stan-Ymet-Xnom1grp-Mrobust.R and Stan-Ymet-Xnom1grp-Mrobust-Example.R
Accompanies the book:
   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
   A Tutorial with R and JAGS, 2nd Edition. Academic Press / Elsevier.


}
\seealso{
\code{...}
}
\examples{
library(lxctk)

## Step 1: compile the model
library(rstan)

# Translate to C++ and compile to DSO. This will take a few mins. Select one of the model_code.
model_code <- student_t_distribution_model()
stanDso <- stan_model( model_code=model_code ) 
## save(stanDso, file='stanDso.RData')  ## For later use

## If you have already had a stanmodel object stanDso saved as a RData file, 
## e.g. stanDso.RData, just load in R with load("stanDso.RData").

## Step 2: running an example
##rstan_options(auto_write = TRUE)
##options(mc.cores = 4)

x <- rt(1000, df=1)
fit <- bayes.student_t(x, stanDso)

## Important notes:
## If bayes.student_t was run in a for-loop (e.g. >128), there will be ERROR:
## "all connections are in use Stan model  does not contain samples."
## Because "A maximum of 128 connections can be allocated..." (?close.connection)
##+for more info. In this case, we need to close connections by hand. For example:

## for (i in 1000)
## {
##	  fit <- bayes.student_t(...)
##	  showConnections(all=TRUE) or ## print(showConnections(all=TRUE))
## }


}

\keyword{utilities}
