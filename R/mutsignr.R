#!/ifshk4/BC_PUB/biosoft/PIPE_RD/Package/R-3.1.1/bin/Rscript

## Access significance level for a noncoding region
## y: number of observed mutations in different mutational categories
## N: number of background mutations in different mutational categories
## bmr: background mutation rate for different mutational categories
## time: number of permutations

core_sig <- function(y, N, bmr, time=500)
{
	#Sg <- dbinom(y, N, bmr, log = TRUE) * (-1)
	#l <- lapply(1:length(N), function(i) dbinom(rbinom(time, N[i], bmr[i]), N[i], bmr[i], log = TRUE) * (-1))
	#permutation.Sg <- rowSums(as.data.frame(l))
	#pval <- sum(permutation.Sg >= Sg) / time

	y.sum <- sum(y)
	l <- lapply(1:length(N), function(i) rbinom(time, N[i], bmr[i]))
	n <- rowSums(as.data.frame(l))
	permutation.p.value <- sum(n>=y.sum) / time
	lambda <- mean(n)
	list(mutation.num=y.sum, lambda=lambda, pois.pvalue=ppois(y.sum,lambda,lower.tail=FALSE), permutation.p.value=permutation.p.value)
}

evaluate.sig <- function(N, y, identifier, mc.cores=4, time=100)
{
	## Background mutation rates
	bmr <- colSums(y) / colSums(N)
	cat("categs:", colnames(y), file=stderr()); cat("\n", file=stderr())
	cat("   bmr:", bmr, file=stderr());cat("\n", file=stderr())
	if (any(bmr == 0)) {
		stop("One of the bmr is ZERO.")
	}
	r <- mclapply(1:nrow(y), function(i,y,N,bmr) core_sig(as.numeric(y[i,]),as.numeric(N[i,]),bmr,time),y,N,bmr,mc.cores=mc.cores)
	d <- as.data.frame(do.call('rbind', r))
	d <- cbind(identifier, d)
	d
}

## Mutsig analysis of somatic mutations for genomic regions
mutsignr <- function(bkgr.mut.file, obs.mut.file, bkgr.mut=NULL,obs.mut=NULL,time=100, mc.cores=4, outfile='output.sig_elements.txt')
{
	## Column names of both input files must be:
	## region  *CpG->T Tp*C->mut       Tp*A->T C->T    C->A    misc    reptime_cluster
	if (!is.null(bkgr.mut) && !is.null(obs.mut)) {
		y <- obs.mut
		N <- bkgr.mut
	} else {
		y <- read.table(obs.mut.file, header=TRUE, stringsAsFactors=FALSE, check.name=FALSE)		
		N <- read.table(bkgr.mut.file, header=TRUE, stringsAsFactors=FALSE, check.name=FALSE)		
	}

	if (any(N$region != y$region)) {
		stop("The `region` cols of 2 input files are not consistent.")
	}
	if (any(N$reptime_cluster != y$reptime_cluster)) {
		stop("The `reptime_cluster` cols of 2 input files are not consistent.")
	}
	region <- N$region
	reptime_cluster <- N$reptime_cluster

	y$region <- NULL
	y$reptime_cluster <- NULL
	N$region <- NULL
	N$reptime_cluster <- NULL

	results <- list()
	clusters <- unique(reptime_cluster)

	for (i in 1:length(clusters)) {
		k <- clusters[i]
		cat("\nreptime_cluster: ", k, ', Ambiguous region:', ifelse(k==0, "Yes", "No"), file=stderr()); cat("\n", file=stderr())
		results[[i]] <- evaluate.sig(N[reptime_cluster==k,], y[reptime_cluster==k,], region[reptime_cluster==k], mc.cores, time)
	}
	x <- do.call('rbind', results)
	x$pois.qvalue <- p.adjust(x$pois.pvalue,method='fdr')
	x <- x[order(x$pois.qvalue,decreasing=FALSE),]
	x <- data.frame(lapply(x, as.character), stringsAsFactors=FALSE)
	write.table(x, file=outfile, quote=FALSE, sep="\t", row.names=FALSE)
}

#library(parallel)
#bkgr.mut <- read.table('bkgr.txt',header=TRUE,stringsAsFactors=FALSE,check.name=FALSE)
#obs.mut <- read.table('obs.txt',header=TRUE,stringsAsFactors=FALSE,check.name=FALSE)
#save(bkgr.mut, obs.mut, file='cds_start_3k_upstream.RData')
#load('cds_start_3k_upstream.RData')
#mutsiggr('bkgr.txt', 'obs.txt', time=100, outfile='test.txt', mc.cores=4)
#mutsiggr(bkgr.mut=bkgr.mut, obs.mut=obs.mut, time=100, outfile='test.txt', mc.cores=4)
