	
	
	
	
#categs.file = "TCGA.Breast_Cancer.MutSigCV.categs.txt"
#coverage.file = "TCGA.Breast_Cancer.MutSigCV.coverage.txt"
#mutations.file = "TCGA.Breast_Cancer.MutSigCV.mutations.txt"
	
mutsig.pathway <- function(categs.file, coverage.file, mutations.file, output.file, pathway=NULL, sep="", method=c("PCT"), exclude.noncoding=FALSE, only.mutated.gene=FALSE, trace=1, cancer.genes=c())
{
	cat("INFO - ", date(), " - Begin to read files...", "\n", sep="", file=stderr())
	mut.df = try(read.table(mutations.file, header=TRUE, sep=sep, stringsAsFactors=FALSE))
	if (class(mut.df) == 'try-error') {
		require(data.table)
		mut.df <- fread(mutations.file)
		setDF(mut.df)
	}
	cov.df = read.table(coverage.file, header=TRUE, stringsAsFactors=FALSE)
	categ.df = read.table(categs.file, header=TRUE, stringsAsFactors=FALSE)
	cat("INFO - ", date(), " - Finished reading files.", "\n", sep="", file=stderr())

	if (exclude.noncoding) {
		mut.df = mut.df[mut.df$effect != "noncoding",]
		cov.df = cov.df[cov.df$effect != "noncoding",]
	}
	mutated.genes=NULL
	if (only.mutated.gene) {
		if (exclude.noncoding) {
			mutated.genes=mut.df$gene[mut.df$effect!="noncoding"]
		} else {
			mutated.genes=mut.df$gene	
		}
	}
	
	sample.number <- length(unique(mut.df$patient))
	
	#genes = c("TP53", "PIK3CA", "GATA3", "PTEN")
	
	# Set BMR from categs.file
	BMR = categ.df$rate * sample.number
	contexts = as.character(categ.df$name)
	names(BMR) = contexts
	context.num = length(contexts)
	
	cat('INFO -', date(), "- contexts: ", contexts, "\n", file=stderr())
	cat("INFO -", date(), "- BMR: ", BMR, "\n", file=stderr())
	
	out.df.names <- c("gene", "pathway", contexts, "obs.total.mut.num", "exp.total.mut.num", "p", "q")
	out.df.attr <- c(character(0),character(0), rep(numeric(0), 8))
	out.df <- as.matrix(as.data.frame(setNames(replicate(length(out.df.names), out.df.attr, simplify=F), out.df.names)))
	
	method <- match.arg(method)
	if (grepl(pattern='perm', method)) {
		perm = rep(0, N)
	}

if (is.null(pathway)) {
	data('LixcCuratedPathway', package='lxctk')
}
for (i in 1:nrow(pathway)) { 
	pathway.name = pathway[i,2]
	pathway.name2 = sub(".html","",basename(pathway[i,1]))
	E = 0
	O = 0
	Q = rep(0, context.num)
	names(Q) = contexts
	s <- as.character(pathway[i, 3])
	genes = unlist(strsplit(s, ","))

	if (only.mutated.gene && !is.null(mutated.genes)) {
		genes=intersect(genes, mutated.genes)
	}

	if (!is.null(cancer.genes)) {
		genes = intersect(genes, cancer.genes)
	}
	if (is.null(genes)) {next;}
	for (gene in genes) {
		cov.gene = cov.df[cov.df$gene==gene,c("categ","coverage")]
		size = try(by(cov.gene$coverage, cov.gene$categ, sum), silent=TRUE)
		if (class(size) == 'try-error') next;
		q = table(mut.df$categ[mut.df$gene==gene])
		if (any(is.na(q[contexts]))) {
			tmp.q <- q[contexts]
			names(tmp.q) = contexts
			tmp.q[is.na(tmp.q)] <- 0
			q <- tmp.q
		}
		if (any(is.na(size[contexts])) && any(is.na(size[contexts]))) {
			cat("WARNING - background coverage for context(s) is NA or 0, skipping...!\n", file=stderr())
			next;
		}
		
		for (e in names(q)) {
			Q[e] = Q[e] + q[e]
		}
		
		if (method == 'PCT') {
			E = E + sum(BMR[contexts] * size[contexts])
			O = O + sum(q)
			#p = ppois(sum(BMR[contexts] * size[contexts]), sum(q), lower.tail=F)  # gene-level significance
		} 	
	}
	if (O==0) {next;}
	p = ppois(O, E, lower.tail=F)
	if (trace) {
		cat(paste("INFO  - ", date()), pathway.name, Q[contexts], O, E, p, sep="\t", file=stderr())
		cat("\n", sep="", file=stderr())
	}
	#out.df <- rbind( out.df, c(gene, q[contexts], sum(q), p, 1) )
	out.df <- rbind( out.df, c(pathway.name2, pathway.name, Q[contexts], O, E, p, 1) )
}
	out.df <- as.data.frame(out.df, stringsAsFactors=FALSE)
	out.df$p <- as.double(out.df$p)
	out.df$q <- p.adjust(out.df$p, method="fdr")
	out.df <- out.df[order(out.df$q), ]
	if (!is.na(output.file)) {
		write.table(out.df, output.file, col.names=out.df.names, row.names=FALSE, quote=FALSE, sep="\t")
	}
	invisible(out.df)
}

