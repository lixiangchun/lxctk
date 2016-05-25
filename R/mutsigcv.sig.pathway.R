

combine.test <- function (p, weight, method = c("fisher", "z.transform", "logit"), hetero = FALSE, na.rm = FALSE)
{
	# The combine.test is directly copied source code of combine.test(...) implemented in R package "survcomp"!
	# See http://www.bioconductor.org/packages/release/bioc/html/survcomp.html
    if (hetero) {
        stop("function to deal with heterogeneity is not implemented yet!")
    }
    method <- match.arg(method)
    na.ix <- is.na(p)
    if (any(na.ix) && !na.rm) {
        stop("missing values are present!")
    }
    if (all(na.ix)) {
        return(NA)
    }
    p <- p[!na.ix]
    k <- length(p)
    if (k == 1) {
        return(p)
    }
    if (missing(weight)) {
        weight <- rep(1, k)
    }
    switch(method, fisher = {
        cp <- pchisq(-2 * sum(log(p)), df = 2 * k, lower.tail = FALSE)
    }, z.transform = {
        z <- qnorm(p, lower.tail = FALSE)
        cp <- pnorm(sum(weight * z)/sqrt(sum(weight^2)), lower.tail = FALSE)
    }, logit = {
        tt <- (-sum(log(p/(1 - p))))/sqrt(k * pi^2 * (5 * k + 
            2)/(3 * (5 * k + 4)))
        cp <- pt(tt, df = 5 * k + 4, lower.tail = FALSE)
    })
    return(cp)
}

mutsigcv.sig.pathway <- function(mutsigcv.sig_gene.filename, pathway, out.filename=NA, exclude.genes=c(), min.n_nonsilent=0)
{
	d <- read.table(mutsigcv.sig_gene.filename, header=TRUE)	
	if (min.n_nonsilent>0) {
		d <- d[d$n_nonsilent>=min.n_nonsilent,]
	}
	if (any(d$p == 0)) {
		min.nonzero.p <- min(d$p[d$p!=0])
		d$p[d$p == 0] <- min.nonzero.p / length(d$p)
	}
	mutsigcv.gene.data <- d[, c('n_nonsilent', 'n_silent', 'n_noncoding', 'p')]
	mutsigcv.genes <- as.character(d$gene)
	rownames(mutsigcv.gene.data) <- mutsigcv.genes

	df = c()
	combine.p.values <- c()
	for (i in 1:nrow(pathway)) {
		s <- as.character(pathway[i, 3])
		pathway.genes <- unlist(strsplit(s, ","))
		genes <- intersect(pathway.genes,mutsigcv.genes)
		genes <- setdiff(genes, exclude.genes)
		p.values <- mutsigcv.gene.data[genes, 'p']
		combine.p.value <- combine.test(p.values)
		mutation.number <- colSums(mutsigcv.gene.data[genes, c('n_nonsilent', 'n_silent', 'n_noncoding')])
		
		pathway.name <- as.character(pathway[i, 1])
		mutated.genes <- paste(genes, collapse=',')
		e <- c(pathway.name, mutated.genes, mutation.number, sum(mutation.number))
		df <- rbind(df, e)
		combine.p.values <- c(combine.p.values, combine.p.value)
	}
	colnames(df) <- c("pathway_name", 'mutated.genes', 'n_nonsilent', 'n_silent', 'n_noncoding', 'total_mutation_number')
	df <- as.data.frame(df)
	df$p.value <- combine.p.values
	df$q.value <- p.adjust(df$p.value, method="fdr")
	df <- df[order(df$q.value),]
	if (!is.na(out.filename)) { write.table(df, out.filename, sep="\t", quote=FALSE, row.names=FALSE) }
	invisible(df)
}

#mutsigcv.sig_gene.filename = "/ifshk4/BC_CANCER/PROJECT/F12HPCNCSZ0222_HUMrduX/lixc/muTect/ana-original-muTect/MutSigCV/EC.sig_genes.txt"
#mutsigcv.sig_gene.filename = "/ifshk4/BC_CANCER/PROJECT/F12HPCNCSZ0222_HUMrduX/lixc/muTect/MutSigCV/TCGA/TCGA.Breast_Cancer.MutSigCV.sig_genes.txt"
#pathway.filename = "/ifshk1/BC_CANCER/03user/lixiangchun/db/pathway/BroadCuratedPathway/classic_cancer_pathway.lxc4"
#pathway.filename = "/ifshk1/BC_CANCER/03user/lixiangchun/db/pathway/BroadCuratedPathway/c2.cp.v4.0.symbols.bin"
#out.filename = "EC.sig_pathways.txt"
#mutsigcv.sig.pathway(mutsigcv.sig_gene.filename, pathway.filename, out.filename)
#mutsigcv.sig.pathway(mutsigcv.sig_gene.filename, "/ifshk1/BC_CANCER/03user/lixiangchun/db/pathway/DNA_repair/DNA_repair_categs", "EC.DNA_repair.sig_genesets.txt")
