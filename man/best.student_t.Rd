\name{best.student_t}
\alias{best.student_t}
\title{Bayesian inference to estimate the parameters of student's t distribution.}
\description{
Bayesian inference to estimate the parameters of student's t distribution.

An example of student's t distirbution:

y ~ student_t( nu , mu , sigma )

model \{

  sigma ~ uniform( unifLo , unifHi ); // prior on sigma

  mu ~ normal( meanY , normalSigma ); // prior on mu

  nuMinusOne ~ exponential( expLambda ); // prior on nu

  y ~ student_t( nu , mu , sigma ); // nu = nuMinusOne + 1;

\}

\code{nu}, \code{mu} and \code{sigma} will be estimated via bayesian inference.

}
\usage{
best.student_t(y, stanDso, meanY=mean(y), sdY=sd(y), 
               unifLo=sdY/1000, unifHi=sdY*1000, 
               normalSigma=sdY*100,
               expLambda=1/abs(meanY-1.0), ...)
}
\arguments{
  \item{y}{A vector of metric values assumed to be generated from student's t distirbution.}
  \item{stanDso}{An object of class \code{stanmodel}.}
  \item{unifLo}{A hyper-parameter for prior distribution for \code{sigma}.}
  \item{unifHi}{A hyper-parameter for prior distribution for \code{sigma}.}
  \item{normalSigma}{A hyper-parameter for prior distribution for \code{mu}.}
  \item{expLambda}{A hyper-parameter for prior distribution for \code{nu}. Must be > 0.}
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
\code{best.robust_t_test}
}
\examples{
library(methods)
library(utils)
library(lxctk)
library(Rcpp)
library(rstan)

## Step 1: compile the model

# Translate to C++ and compile to DSO. This will take a few mins. Select one of the model_code.
model_code <- best.student_t_model()
stanDso <- stan_model( model_code=model_code ) 
## save(stanDso, file='stanDso.RData')  ## For later use

## If you have already had a stanmodel object stanDso saved as a RData file, 
## e.g. stanDso.RData, just load in R with load("stanDso.RData").

## Step 2: running an example
##rstan_options(auto_write = TRUE)
##options(mc.cores = 4)

x <- rt(1000, df=1)
fit <- best.student_t(x, stanDso)

## Important note 1:
## If best.student_t was run in a for-loop (e.g. >128), there will be ERROR:
## "all connections are in use Stan model  does not contain samples."
## Because "A maximum of 128 connections can be allocated..." (?close.connection)
##+for more info. In this case, we need to close connections by hand. For example:

## for (i in 1000)
## {
##	  fit <- best.student_t(...)
##	  showConnections(all=TRUE) or ## print(showConnections(all=TRUE))
## }

## Important note 2:
## If the following error occurs, my recipe is to set 'unifLo' to small value and 
##+'unifHi' to very large number.
#SAMPLING FOR MODEL 'ede2285102cdfeb5c5e92e95855f30cc' NOW (CHAIN 4).
# [1] "  Log probability evaluates to log(0), i.e. negative infinity."                                       
# [2] "  Stan can't start sampling from this initial value."                                                 
# [3] "Rejecting initial value:"                                                                             
# [4] "  Log probability evaluates to log(0), i.e. negative infinity."                                       
# [5] "  Stan can't start sampling from this initial value."                                                 
# [6] "Rejecting initial value:"                                                                             
# [7] "  Log probability evaluates to log(0), i.e. negative infinity."                                       
# [8] "  Stan can't start sampling from this initial value."                                                 
# [9] "Initialization between (-2, 2) failed after 100 attempts. "                                           
#[10] " Try specifying initial values, reducing ranges of constrained values, or reparameterizing the model."
#error occurred during calling the sampler; sampling not done

data('x')
print(x)
fit <- best.student_t(x, stanDso) ## The above error will be generated
fit <- best.student_t(x, stanDso, unifLo=0, unifHi=.Machine$double.xmax) ## No errors

}

\keyword{utilities}
