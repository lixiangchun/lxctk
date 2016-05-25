\name{unicoxph}
\alias{unicoxph}
\title{Perform univariate Cox Regression Analysis for multipe variables sequentially.}
\description{
Perform univariate Cox Regression Analysis for multipe variables sequentially.
}
\usage{
unicoxph(cox.f,data,variable.display=NULL, plot=TRUE, ...)
}
\arguments{
  \item{cox.f}{A formula passed to \code{coxph(...)}, multiple variables on the right side supported.}
  \item{data}{A data frame includes all variables in the \code{cox.f}.}
  \item{variable.display}{Variable names expected to display.}
  \item{plot}{Whether to plot univariate Cox regression analysis results.}
}
\details{
See \code{examples} to check how to prepare a coxTable.
}
\seealso{
	\code{plot.coxph},\code{forestplot},\code{survdiff},\code{survfit},\code{coxph}
}
\examples{

data("crcdat", package='lxctk')

## Note: relevel all related variables in the data frame; relevel variables in the
##+formula does not work because use `all.vars(...)` to extract variables from a
##+formula. However, the relevel manipulation cannot extract with `all.vars`. So
##+the use must relevel variables in the data frame, or the default rule will be
##+used.
crcdat$MSI <- as.factor(crcdat$MSI>0)
crcdat$Gender <- as.character(crcdat$Gender)
crcdat$Gender <- relevel(as.factor(crcdat$Gender), ref='F')
crcdat$TNM <- relevel(as.factor(crcdat$TNM), ref=1)
crcdat$Tumor.localisation <- relevel(as.factor(crcdat$Tumor.localisation), ref=1)
crcdat$signature <- as.factor(crcdat$signature)
crcdat$Kras_G12_13X <- relevel(as.factor(crcdat$Kras_G12_13X), ref=0)

f <- Surv(Survival.month, status) ~ Age+Gender+MSI+TNM+Tumor.localisation+hypermutated+Kras_G12_13X+signature
unicoxph(f, crcdat)

variable.display <- c('Age', 'Gender (M vs. F)', 'MSI (+ vs. -)', 'TNM (2 vs. 1)', 'TNM (3 vs. 1)', 'TNM (4 vs. 1)', 'Tumor site (2 vs. 1)', 'Hypermutated (True vs. False)', 'Kras_G12_13X (True vs. False)', 'Signature (True vs. False)')
unicoxph(f, crcdat, variable.display)
title(main='An example of unicoxph', sub='Unicoxph')

}

\keyword{utilities}
