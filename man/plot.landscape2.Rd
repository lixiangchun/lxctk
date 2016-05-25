\name{plot.landscape2}
\alias{plot.landscape2}
\title{Visualize mutation landscape}
\description{
A function to plot mutational landscape of significantly mutated genes.
}
\usage{
plot.landscape2(upper.panel, middle.panel, left.panel, right.panel, legend.panel, subtype.panel=NULL, subtype.2.panel=NULL,
    output.fig.pdffile=NULL,
    left.panel.xlim=NULL,
    family.font="sans", 
    syn.nonsyn.legend.color=NULL,
    syn.nonsyn.legend.ncol=NULL,
    syn.nonsyn.legend.loc="top",
    image.color=NULL,

    subtype.labels=NULL,
    subtype.label.cex=1,
    subtype.panel.color=brewer.pal(7, "Dark2"),
    subtype.panel.identifier="",                        # subtype identifier, e.g. Lauren's subtype, Clonality status
    subtype.panel.identifier.cex=0.6,                   # cex for subtype identifier
    subtype.legend.text.width=NULL,                     # text.width passed to legend(...) to adjust distances between 2 adjacent legend strings

	subtype.2.labels=NULL,
    subtype.2.label.cex=1,
    subtype.2.panel.color=brewer.pal(7, "Set2"),
    subtype.2.panel.identifier="",                      # subtype identifier, e.g. Lauren's subtype, Clonality status
    subtype.2.panel.identifier.cex=0.6,                 # cex for subtype identifier
    subtype.2.legend.text.width=NULL,                   # text.width passed to legend(...) to adjust distances between 2 adjacent legend strings

    right.panel.color=brewer.pal(7, "Purples")[5],
    sort.by.qvalue=FALSE,
    group.by.subtype=FALSE,
    left.panel.cex.axis=1.5,
    right.panel.cex.axis=1.3,
    upper.panel.mar=c(0.2, 1.3, 0.5, 1.16),
    legend.panel.mar=c(0, 1, 2, 1),
    legend.panel.ncol=NULL,
    legend.panel.y=-0.05,
    left.panel.mar=c(5.5, 1.5, 0.2, 2.9),
    middle.panel.mar=c(5.5, 1.3, 0.2, 1.16),
    right.panel.mar=c(5.5, 3, 0.2, 1.5),
    left.panel.xlab="Individuals\nwith mutations",
    left.panel.xlab.cex=1,
    right.panel.xlab="-log10\n(q-value)",
    right.panel.xlab.cex=1,
    right.panel.qscore.cutoff=1,
    right.panel.qscore.lwd=0.5,
    right.panel.qscore.lty='dashed',
    right.panel.qscore.col='red')
}
\arguments{
  \item{upper.panel}{A data frame for upper panel barplot, see \code{details} for detail description of file format.}
  \item{middle.panel}{A matrix for middle panel, attached with row and col names, see \code{details} for detail description of file format.}
	\item{left.panel}{A data frame for the left panel, see \code{details} for detail description of file format.}
	\item{right.panel}{A data frame for the left panel, see \code{details} for detail description of file format.}
	\item{legend.panel}{A data frame contains information about mutation type mapping.}
	\item{subtype.panel}{The main subtype panel, a data frame of sample sutbypes, 1st col is sample name, whereas 2nd col is subtype.}
	\item{subtype.2.panel}{The 2nd subtype panel, a data frame of sample sutbypes, 1st col is sample name, whereas 2nd col is subtype.}
	\item{output.fig.pdffile}{figure file name, ended with \code{.pdf}.}
	\item{subtype.labels}{text labels for sample subtypes, e.g. \code{Intestinal} and \code{Diffuse}.}
	\item{subtype.2.labels}{text labels for sample subtypes, e.g. \code{Intestinal} and \code{Diffuse}.}
	\item{group.by.subtype}{if TRUE, group samples by subtypes.}
	\item{syn.nonsyn.legend.color}{A vector of 2 elements indicating colors used for syn. and non-syn. mutations.}
	\item{syn.nonsyn.legend.ncol}{\code{ncol} passes to \code{legend(...)}.}
	\item{syn.nonsyn.legend.loc}{Legend location, e.g. \code{top}, \code{right} etc.}
	\item{image.color}{A vector of elements representing different mutation types.}
	\item{subtype.panel.color}{A vector representing colors used for subtypes.}
	\item{right.panel.color}{Color for barplot on the right panel.}
	\item{sort.by.qvalue}{Sort right panel data frame according to q-value.}
	\item{left.panel.xlim}{xlim for left panel.}
	\item{family.font}{font family used set in par(family=family.font).}
	\item{left.panel.cex.axis}{cex.axis for the left panel at side=4, control the size of percentage.}
	\item{right.panel.cex.axis}{cex.axis for right panel at side=1, control the size of gene names.}
	\item{upper.panel.mar}{mar in par(mar=c(...)) for upper panel.}
	\item{legend.panel.mar}{mar in par(mar=c(...)) for top right legend panel.}
	\item{legend.panel.ncol}{\code{ncol} passes to \code{legend(...)}.}
	\item{legend.panel.y}{Y coordinate for legend.panel.}
	\item{left.panel.mar}{mar in par(mar=c(...)) for left panel.}
	\item{middle.panel.mar}{mar in par(mar=c(...)) for middle panel.}
	\item{right.panel.mar}{mar in par(mar=c(...)) for right panel.}
	\item{legend.syn.nonsyn.x,legend.syn.nonsyn.y}{set for top left legend via legend(x=\code{legend.syn.nonsyn.x},y=\code{legend.syn.nonsyn.y}).}
	\item{left.panel.xlab}{left panel \code{xlab}.}
	\item{left.panel.xlab.cex}{\code{cex} for left panel \code{xlab}.}
	\item{right.panel.xlab}{right panel \code{xlab}.}
	\item{right.panel.xlab.cex}{\code{cex} for right panel \code{xlab}.}
	\item{right.panel.qscore.cutoff}{Q-score cutoff.}
	\item{right.panel.qscore.lwd}{lwd for qscore cutoff, set 0 to get rid of it.}
	\item{right.panel.qscore.lty}{lty for qscore cutoff.}
	\item{right.panel.qscore.col}{col for qscore cutoff.}
}
\details{
\code{upper.panel} has 3 columns, the 1st is sample names, the 2nd and 3rd are the number of silent and non-silent mutations, e.g.
*Note:It's highly recommended NOT to include '-' in sample names

sample  silence nonsilence\cr
GC03_T  167     551\cr
GC10_T  179     471\cr
GC04_T  160     388\cr
GC05_T  156     391

\code{middle.panel} is a matrix, attached with row and columns, e.g.

GC118_T GC191_T GC178_T GC101_T\cr
TP53 0 0 0 0\cr
ARID1A 0 0 2 0\cr
CDH1 6 0 0 3\cr
PTK2B 0 0 0 0\cr
RPL22 0 0 5 0\cr
BAX 0 0 0 0\cr
CYC1 0 0 0 0\cr
KRAS 0 0 0 0\cr
ITGAV 0 0 0 0\cr
HIST1H2AD 0 0 3 0

\code{left.panel} is a data frame with 3 columns, e.g.

gene    nonsilence_sample_num   nonsilent_sample_percent\cr
TP53    45      37.2\%\cr
ARID1A  28      23.1\%\cr
CDH1    21      17.4\%\cr
PTK2B   17      14.0\%\cr
RPL22   12      9.9\%\cr
BAX     11      9.1\%

\code{right.panel} is a data frame with 2 columns, e.g.

gene    q\cr
TP53    1.73999e-07\cr
ARID1A  0.002062569\cr
CDH1    2.875046e-06\cr
PTK2B   0.001464807\cr
RPL22   1.73999e-07\cr
BAX     1.73999e-07

\code{legend.panel} is a data frame with 2 columns, e.g.

7	Nonsense\cr
6	Splice site\cr
5	Frame shift\cr
4	Inframe indel\cr
3	Missense\cr
2	Syn.

You can also use the pre-specified image color as:

data("LandscapeColor")\cr
ls()

}
\seealso{
  \code{plot.depth},\code{plot.MH}
}
\examples{
library(RColorBrewer)
data("SMG", package="lxctk")
data("LandscapeColor")

plot.landscape2(upper.panel, middle.panel, left.panel, right.panel, legend.panel, subtype.panel, subtype.labels=c('Intestinal','Diffuse'), image.color=image.color.lixc)

# Use different color scheme
plot.landscape2(upper.panel, middle.panel, left.panel, right.panel, legend.panel, subtype.panel, subtype.labels=c('Intestinal','Diffuse'), image.color=image.color.Set2)


plot.landscape2(upper.panel, middle.panel, left.panel, right.panel, legend.panel, subtype.panel, subtype.labels=c('Intestinal','Diffuse'),output.fig.pdffile="GastricCancer.SMG.Landscape.pdf")
}

\keyword{utilities}
