\name{beanplot}
\alias{beanplot}

\title{Beanplot: A Boxplot/Stripchart/Vioplot Alternative}

\description{
  Plots beans to compare the distributions of different groups;
  it draws one bean per group of data. A bean consists of a
  one-dimensional scatter plot, its distribution as a density
  shape and an average line for the distribution. Next to that,
  an overall average for the whole plot is drawn per default.
}

\usage{
beanplot(..., bw = "SJ-dpi", kernel = "gaussian", cut = 3, cutmin = -Inf, 
         cutmax = Inf, grownage = 10, what = c(TRUE, TRUE, TRUE, TRUE), 
         add = FALSE, col, axes = TRUE, log = "auto", handlelog = NA, 
         ll = 0.16, wd = NA, maxwidth = 0.8, maxstripline = 0.96, 
         method = "stack", names, overallline = "mean", beanlines = overallline, 
         horizontal = FALSE, side = "no", jitter = NULL, beanlinewd = 2, 
         frame.plot = axes, border = NULL, innerborder = NA, at = NULL, 
         boxwex = 1, ylim = NULL, xlim = NULL, show.names = NA)
}

\arguments{
  \item{\dots}{data which to perform the beanplot on. This data can consist of dataframes, vectors and/or formulas. For each formula, a dataset can be specified with data=[dataset], and a subset can be specified with subset=[subset]. If subset/data arguments are passed, but there are not enough subset/data arguments, they are reused. Additionally, na.action, drop.unused.levels and xlev can be passed to \code{\link{model.frame}} in the same way. Also, parameters for \code{\link{axis}} and \code{\link{title}} can be passed.}
  \item{bw}{the bandwidth (method) being used, used by \code{\link{density}}. In case of a method, the average computed bandwidth is used.}
  \item{kernel}{see \code{\link{density}}.}
  \item{cut}{the beans are cut beyond \code{cut*bw}}
  \item{cutmin}{the low-ends of the beans are cut below \code{mincut*bw}. Defaults to \code{cut}.}
  \item{cutmax}{the high-ends of the beans are cut beyond \code{maxcut*bw}. Defaults to \code{cut}.}
  \item{grownage}{the width of a bean grows linearly with the count of points, until \code{grownage} is reached.}
  \item{what}{a vector of four booleans describing what to plot. In the following order, these booleans stand for the total average line, the beans, the bean average, and the beanlines. For example, \code{what=c(0,0,0,1)} produces a \code{\link{stripchart}}}
  \item{add}{if true, do not start a new plot}
  \item{col}{the colors to be used. A vector of up to four colors can be used. In the following order, these colors stand for the area of the beans (without the border, use \code{border} for that color), the lines inside the bean, the lines outside the bean, and the average line per bean. Transparent colors are supported. \code{col} can be a list of color vectors, the vectors are then used per bean, and reused if necessary.}
  \item{axes}{if false, no axes are drawn.}
  \item{log}{use \code{log="y"} or \code{log=""} to force a log-axis. In case of \code{log="auto"}, a log-transformation is used if appropriate}
  \item{handlelog}{if \code{handlelog} then all the calculations are done using a log-scale. By default this is determined using the log parameter.}
  \item{ll}{the length of the beanline per point found.}
  \item{wd}{the linear transformation that determines the width of the beans. By default determined using \code{maxwidth}, and returned.}
  \item{maxwidth}{the maximum width of a bean.}
  \item{maxstripline}{the maximum length of a beanline.}
  \item{method}{the method used when two points on a bean are the same. \code{"stack"}, \code{"overplot"} and \code{"jitter"} are supported.}
  \item{names}{a vector of names for the groups.}
  \item{overallline}{the method used for determining the overall line. Defaults to \code{"mean"}, \code{"median"} is also supported.}
  \item{beanlines}{the method used for determining the average bean line(s). Defaults to \code{"mean"}, \code{"median"} and \code{"quantiles"} are also supported.}
  \item{horizontal}{if true, the beanplot is horizontal}
  \item{side}{the side on which the beans are plot. Default is \code{"no"}, for symmetric beans. \code{"first"}, \code{"second"} and \code{"both"} are also supported.}
  \item{jitter}{passed to \code{\link{jitter}} as amount in case of \code{method="jitter"}.}
  \item{beanlinewd}{the width used for the average bean line}
  \item{frame.plot}{if true, plots a frame.}
  \item{border}{the color for the border around a bean. \code{NULL} for \code{par("fg")}, \code{NA} for no border. If border is a vector, it specifies the color per bean, and colors are reused if necessary.}
  \item{innerborder}{a color (vector) for the border inside the bean(s). Especially useful if \code{side="both"}. Use \code{NA} for no inner border. Colors are reused if necessary. The inner border is drawn as the last step in the drawing process.}
  \item{at}{the positions at which a bean should be drawn.}
  \item{boxwex}{a scale factor applied to all beans. Compatible with \code{\link{boxplot}}.}
  \item{ylim}{the range to plot.}
  \item{xlim}{the range to plot the beans at.}
  \item{show.names}{if true, plots the names as axis labels}
}

\details{
  Most parameters are compatible with \code{\link{boxplot}} and
  \code{\link{stripchart}}. For compatibility, arguments with the
  name "formula" or "x" are used as data. However, data or formulas
  do not need to be named "x" or "formula". The function handles
  (combinations of) dataframes, vectors and/or formulas.
}

\value{
  \item{bw}{The bandwith (bw) used.}
  \item{wd}{The bean width (wd) used.}
}

\references{
  Kampstra, P. (2008) Beanplot: A Boxplot Alternative for Visual Comparison of Distributions.
  \emph{Journal of Statistical Software, Code Snippets}, \bold{28}(1), 1-9.
  URL \url{http://www.jstatsoft.org/v28/c01/}
}

\author{
  Peter Kampstra <pkampst=beanplot@cs.vu.nl>
}

\note{
  In case of more than 5000 values per bean, the autodetection of log fails. In such cases, use
  \code{log=""} or \code{log="y"}.
}

\seealso{
  \code{\link{boxplot}}, \code{\link{stripchart}}, \code{\link{density}},
  \code{\link{rug}}, \code{\link[vioplot]{vioplot}} in package \code{vioplot}
}

\examples{
beanplot(rnorm(22),rnorm(22),rnorm(22),main="Test!",rnorm(3))

#mostly examples taken from boxplot:
par(mfrow = c(1,2))
boxplot(count ~ spray, data = InsectSprays, col = "lightgray")
beanplot(count ~ spray, data = InsectSprays, col = "lightgray", border = "grey", cutmin = 0)

boxplot(count ~ spray, data = InsectSprays, col = "lightgray")
beanplot(count ~ spray, data = InsectSprays, col = "lightgray", border = "grey", overallline = "median")

boxplot(decrease ~ treatment, data = OrchardSprays,
        log = "y", col = "bisque", ylim = c(1,200))
beanplot(decrease ~ treatment, data = OrchardSprays,
        col = "bisque", ylim = c(1,200))

par(mfrow = c(2,1))
mat <- cbind(Uni05 = (1:100)/21, Norm = rnorm(100),
             T5 = rt(100, df = 5), Gam2 = rgamma(100, shape = 2))
par(las=1)# all axis labels horizontal
boxplot(data.frame(mat), main = "boxplot(*, horizontal = TRUE)",
        horizontal = TRUE, ylim = c(-5,8))
beanplot(data.frame(mat), main = "beanplot(*, horizontal = TRUE)",
        horizontal = TRUE, ylim = c(-5,8))

par(mfrow = c(1,2))
boxplot(len ~ dose, data = ToothGrowth,
        boxwex = 0.25, at = 1:3 - 0.2,
        subset = supp == "VC", col = "yellow",
        main = "Guinea Pigs' Tooth Growth",
        xlab = "Vitamin C dose mg",
        ylab = "tooth length", ylim = c(-1, 40), yaxs = "i")
boxplot(len ~ dose, data = ToothGrowth, add = TRUE,
        boxwex = 0.25, at = 1:3 + 0.2,
        subset = supp == "OJ", col = "orange")
legend("bottomright", bty="n", c("Ascorbic acid", "Orange juice"),
       fill = c("yellow", "orange"))
allplot <- beanplot(len ~ dose+supp, data = ToothGrowth, 
    what=c(TRUE,FALSE,FALSE,FALSE),show.names=FALSE,ylim=c(-1,40), yaxs = "i")
beanplot(len ~ dose, data = ToothGrowth, add=TRUE,
        boxwex = 0.6, at = 1:3*2 - 0.9,
        subset = supp == "VC", col = "yellow",border="yellow2",
        main = "Guinea Pigs' Tooth Growth",
        xlab = "Vitamin C dose mg",
        ylab = "tooth length", ylim = c(3, 40), yaxs = "i",
        bw = allplot$bw, wd = allplot$wd, what = c(FALSE,TRUE,TRUE,TRUE))
beanplot(len ~ dose, data = ToothGrowth, add = TRUE,
        boxwex = 0.6, at = 1:3*2-0.1,
        subset = supp == "OJ", col = "orange",border="darkorange",
        bw = allplot$bw, wd = allplot$wd, what = c(FALSE,TRUE,TRUE,TRUE))
legend("bottomright", bty="n", c("Ascorbic acid", "Orange juice"),
       fill = c("yellow", "orange"))

par(mfrow = c(1,2))
boxplot(len ~ dose, data = ToothGrowth,
        boxwex = 0.25, at = 1:3 - 0.2,
        subset = supp == "VC", col = "yellow",
        main = "Guinea Pigs' Tooth Growth",
        xlab = "Vitamin C dose mg",
        ylab = "tooth length", ylim = c(-1, 40), yaxs = "i")
boxplot(len ~ dose, data = ToothGrowth, add = TRUE,
        boxwex = 0.25, at = 1:3 + 0.2,
        subset = supp == "OJ", col = "orange")
legend("bottomright", bty="n",c("Ascorbic acid", "Orange juice"),
       fill = c("yellow", "orange"))
beanplot(len ~ reorder(supp, len, mean) * dose, ToothGrowth, 
        side = "b", col = list("yellow", "orange"), border = c("yellow2", 
            "darkorange"), main = "Guinea Pigs' Tooth Growth", 
        xlab = "Vitamin C dose mg", ylab = "tooth length", ylim = c(-1, 
            40), yaxs = "i")
legend("bottomright", bty="n",c("Ascorbic acid", "Orange juice"),
       fill = c("yellow", "orange"))

#Example with multiple vectors and/or formulas
par(mfrow = c(2,1))
beanplot(list(all = ToothGrowth$len), len ~ supp, ToothGrowth, len ~ dose)
title("Tooth growth length (beanplot)")
#Trick using internal functions to do this with other functions:
mboxplot <- function(...){
    graphics::boxplot(beanplot:::getgroupsfromarguments(), ...)
}
mstripchart <- function(..., method = "overplot", jitter = 0.1, offset = 1/3,
           vertical = TRUE, group.names, add = FALSE,
           at = NULL, xlim = NULL, ylim = NULL,
           ylab = NULL, xlab=NULL, dlab = "", glab = "",
           log = "", pch = 0, col = par("fg"), cex = par("cex"), 
           axes = TRUE, frame.plot = axes) {
    graphics::stripchart(beanplot:::getgroupsfromarguments(),
           method, jitter, offset, vertical, group.names, add,
           at, xlim, ylim, ylab, xlab, dlab, glab, log, pch, col, cex, 
           axes, frame.plot)
}
mstripchart(list(all = ToothGrowth$len), len ~ supp, ToothGrowth, len ~ dose,
            xlim = c(0.5,6.5))
title("Tooth growth length (stripchart)")
}
\keyword{hplot}
