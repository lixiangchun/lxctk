
# filename
#++At least has two columns, named "gene" and "value"

# Number of genes to be plotted, set to NA to plot all genes
#++k <- NA

plot.MH <- function(df, refgene, label.genes=NA, k=400, xmin=4, xmax=300, ymin=4, ymax=268, magnitude=40, shift=2, z.padj=0.5, cex=1)
{
names(df) <- c('gene', 'value', 'qvalue')
df <- df[order(df$value, decreasing=TRUE),]
df$value <- df$value * magnitude
df$gene <- as.character(df$gene)

if (is.na(label.genes) || length(label.genes) == 0) {
	#label.genes <- df$gene[1:6]
	label.genes <- df$gene[df$value>0.09 * magnitude]
}

if (!is.na(k) & nrow(df) > k) {
	df <- df[1:k,]
}
df.n <- nrow(df)

if (k > xmax * ymax) {
	warning("k > xmax * ymax")
}

# The length of ech chromosome, e.g. chromosome.length[1] is the length of chr1, ...
chromosome.length <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
chromosome.index = 1:length(chromosome.length)

col <- matrix(nrow=xmax+shift, ncol=ymax+shift, "forestgreen")
z <- matrix(nrow=xmax+shift, ncol=ymax+shift, 0)

VALUE <- c()
X <- c()
Y <- c()
Z <- c()
texts <- c()

for (i in 1:df.n) {
	gene <- df$gene[i]
	value <- df$value[i]
	e <- refgene[refgene$gene==gene, c("chrom", "pos")]
	chrom0 <- as.integer(e[1])
	pos0 <- as.integer(e[2])
	chrom0.len <- chromosome.length[chrom0]

	x <- (chrom0 - min(chromosome.index)) / (max(chromosome.index) - min(chromosome.index)) * (xmax - xmin) + xmin
	y <- (pos0 - 0) / (chrom0.len - 0) * (ymax - ymin) + ymin
	xi <- as.integer(round(x))
	yi <- as.integer(round(y))
	
	z[xi,yi] <- value
	if (any(label.genes==gene)) {            # if gene is the one that is to be plotted
		#col[c(xi-1, xi, xi+1), c(yi-1,yi,yi+1)] <- colors()[459]
		VALUE <- c(VALUE, value)
		X <- c(X, xi)
		Y <- c(Y, yi)
		Z <- c(Z, value+z.padj)
		texts <- c(texts, gene)
		#cat(gene, value, "\n", sep="\t")
	} else {
		col[xi, yi] <- "orange"
		#col[c(xi-1, xi, xi+1), c(yi-1,yi,yi+1)] <- "orange"
	}
	#cat(chrom0, pos0, gene, value, chrom0.len,  xi, yi, "\n", sep="\t")
}

# It is necessary to assign value to genes to be colored and labeled at the last step!
for (i in 1:length(VALUE)) {
	z[X[i], Y[i]] <- VALUE[i]
	col[c(X[i]-1, X[i], X[i]+1), c(Y[i]-1,Y[i],Y[i]+1)] <- colors()[459]
}

library(rgl)
#surface3d(1:(xmax+shift), 1:(ymax+shift), z, color=col, back='cull')
surface3d(1:dim(z)[1], 1:dim(z)[2], z, color=col, back='cull')
texts3d(X, Y, Z, texts=texts, cex=cex,col='black', font=3, useFreeType=TRUE)

}

