\name{plot.coxph}
\alias{plot.coxph}
\title{Draw a forestplot to visualize multivariate Cox regression analysis results.}
\description{
Draw a forestplot to visualize multivariate Cox regression analysis results.
}
\usage{
plot.coxph(coxTable, coxObj=NULL, clip=c(0,8), variable.display=NULL, ...)
}
\arguments{
  \item{coxTable}{A multivariate Cox regression table in user specified format.}
  \item{coxObj}{A multivariate Cox regression object from \code{coxph}.}
  \item{clip}{CI out of which is clipped.}
  \item{variable.display}{Variable names expected to display.}
  \item{...}{More variables passed to \code{forestplot}, e.g. \code{xlog}.}
}
\details{
See \code{examples} to check how to prepare a coxTable.
}
\seealso{
	\code{forestplot},\code{survdiff},\code{survfit},\code{coxph}
}
\examples{
data("plot.coxph", package='lxctk')
plot.coxph(coxTable=coxTable)
plot.coxph(coxObj=coxObj)

########Set user specified variable names############
data("crcdat", package='lxctk')
crcdat$MSI <- as.factor(crcdat$MSI>0)
crcdat$Gender <- as.character(crcdat$Gender)
crcdat$Gender <- relevel(as.factor(crcdat$Gender), ref='F')
crcdat$TNM <- relevel(as.factor(crcdat$TNM), ref=1)
crcdat$Tumor.localisation <- relevel(as.factor(crcdat$Tumor.localisation), ref=1)
crcdat$signature <- as.factor(crcdat$signature)
crcdat$Kras_G12_13X <- relevel(as.factor(crcdat$Kras_G12_13X), ref=0)

f <- Surv(Survival.month, status) ~ Age+Gender+MSI+TNM+Tumor.localisation+hypermutated+Kras_G12_13X+signature
r <- coxph(f, crcdat)
plot.coxph(coxObj=r)

########Set user specified variable names############
variable.display <- c('Age', 'Gender (M vs. F)', 'MSI (+ vs. -)', 'TNM (2 vs. 1)', 'TNM (3 vs. 1)', 'TNM (4 vs. 1)', 'Tumor site (2 vs. 1)', 'Hypermutated (True vs. False)', 'Kras_G12_13X (True vs. False)', 'Signature (True vs. False)')
plot.coxph(coxObj=r, variable.display=variable.display)
}

\keyword{utilities}
