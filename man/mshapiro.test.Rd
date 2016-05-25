%	$Id: mshapiro.test.Rd,v 1.10 2009-02-09 22:25:24 sjarek Exp $	
\name{mshapiro.test}
\alias{mshapiro.test}
\encoding{latin1}
\title{Shapiro-Wilk Multivariate Normality Test. Source code from package \code{mvnormtest}.}
\description{
  Performs the Shapiro-Wilk test for multivariate normality.
}
\usage{
mshapiro.test(U)
}
\arguments{
  \item{U}{a numeric matrix of data values, the number of which must be 
    for each sample between 3 and 5000. }
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the Shapiro-Wilk statistic.}
  \item{p.value}{the p-value for the test.}
  \item{method}{the character string \code{"Shapiro-Wilk normality test"}.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\author{ Slawomir Jarek (\email{slawomir.jarek@gallus.edu.pl})}
\references{
  Czeslaw Domanski (1998)
  Wlasnosci testu wielowymiarowej normalnosci Shapiro-Wilka i jego 
  zastosowanie.
  \emph{Cracow University of Economics Rector's Lectures}, \bold{No. 37}.
  
  Patrick Royston (1982)
  An Extension of Shapiro and Wilk's \eqn{W} Test for Normality to Large
  Samples.
  \emph{Applied Statistics}, \bold{31}, 115--124.

  Patrick Royston (1982)
  Algorithm AS 181: The \eqn{W} Test for Normality.
  \emph{Applied Statistics}, \bold{31}, 176--180.

  Patrick Royston (1995)
  A Remark on Algorithm AS 181: The \eqn{W} Test for Normality.
  \emph{Applied Statistics}, \bold{44}, 547--551.
}
\seealso{
  \code{\link{shapiro.test}} for univariate samples,
  \code{\link[stats:qqnorm]{qqnorm}} for producing a normal quantile-quantile plot.  
}
\examples{
##library(mvnormtest)
data(EuStockMarkets)

C <- t(EuStockMarkets[15:29,1:4])
mshapiro.test(C)

C <- t(EuStockMarkets[14:29,1:4])
mshapiro.test(C)

R <- t(diff(t(log(C))))
mshapiro.test(R)

dR <- t(diff(t(R)))
mshapiro.test(dR)

}
\keyword{htest}
