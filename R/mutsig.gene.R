	
	
	
	
#categs.file = "TCGA.Breast_Cancer.MutSigCV.categs.txt"
#coverage.file = "TCGA.Breast_Cancer.MutSigCV.coverage.txt"
#mutations.file = "TCGA.Breast_Cancer.MutSigCV.mutations.txt"
	
mutsig.gene <- function(categs.file, coverage.file, mutations.file, mut.df=NULL, output.file, sep="", method=c("PCT","FCPT","LRT","perm.score","perm.num"), exclude.noncoding=FALSE, N=10000, trace=FALSE, only.CGC=FALSE, GOI=NULL, mutations.file.formatted_as_csv=FALSE)
{
	if (is.null(mut.df)) {
		if (mutations.file.formatted_as_csv) {
			mut.df <- read.csv(mutations.file)
		} else {
			mut.df = try(read.table(mutations.file, header=TRUE, sep=sep, stringsAsFactors=FALSE))
			if (class(mut.df)=='try-error') {
				print("You may need to set sep=`\t`, or formatted mutations.file as csv and set mutations.file.formatted_as_csv=TRUE.")
				quit('no')
			}
		}
	}

	cov.df = read.table(coverage.file, header=TRUE, stringsAsFactors=FALSE)
	categ.df = read.table(categs.file, header=TRUE, stringsAsFactors=FALSE)

	if (exclude.noncoding) {
		mut.df = mut.df[mut.df$effect != "noncoding",]
		cov.df = cov.df[cov.df$effect != "noncoding",]
	}
	
	sample.number <- length(unique(mut.df$patient))
	genes = unique(mut.df$gene)
	if (!is.null(GOI)) {
		only.CGC=FALSE
		genes = intersect(genes, GOI)
	}
	if (only.CGC) {
		data('CGC20140319', package='lxctk')
		genes = intersect(genes, CGC$Symbol)
	}
	#genes = c("TP53", "PIK3CA", "GATA3", "PTEN")
	
	# Set BMR from categs.file
	BMR = categ.df$rate * sample.number
	contexts = as.character(categ.df$name)
	names(BMR) = contexts
	context.num = length(contexts)
	
	out.df.names <- c("gene", contexts, "total.mut.num", "p", "q")
	out.df.attr <- c(character(0), rep(numeric(0), 8))
	out.df <- as.matrix(as.data.frame(setNames(replicate(length(out.df.names), out.df.attr, simplify=F), out.df.names)))
	
	if (trace) {
		cat("log.time", "gene", contexts, "total.mut.num\tp", sep="\t")
		cat("\n")
	}
	
	method <- match.arg(method)
	if (grepl(pattern='perm', method)) {
		perm = rep(0, N)
	}
	
	for (gene in genes) {
		cov.gene = cov.df[cov.df$gene==gene,c("categ","coverage")]
		size = try(by(cov.gene$coverage, cov.gene$categ, sum), silent=TRUE)
		if (class(size) == 'try-error') next;
		#q = table(mut.df$categ[mut.df$gene==gene])[2:6]
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
			#tmp.size <- size[contexts]
			#names(tmp.size) = contexts
			#tmp.size[is.na(tmp.size)] <- 0
			#size <- tmp.size
		}
		
		if (method == 'PCT') {
			E = sum(BMR[contexts] * size[contexts])
			O = sum(q)
			p= ppois(O, E, lower.tail=F)
		} else if (method == 'FCPT') {
			####################MuSic (Fisher’s combined P-value test (FCPT))###########################
			X.c = -2 * sum(pbinom(q[contexts], size[contexts], BMR[contexts], lower.tail=FALSE, log.p=TRUE))
			p= pchisq(X.c, context.num * 2, lower.tail=FALSE)
			#######################################################################
		} else if (method == 'LRT') {
			####################MuSiC (Likelihodd ratio test (i.e. LRT))###########
			bmr = q[contexts] / size[contexts]
			log.lik.bmr = dbinom(q[contexts], size[contexts], bmr[contexts], log=TRUE)
			log.lik.BMR = dbinom(q[contexts], size[contexts], BMR[contexts], log=TRUE)
			X.l = 2 * (log.lik.bmr - log.lik.BMR)
			p= pchisq(X.l, context.num, lower.tail=FALSE)
			#######################################################################
		} else if (method == 'perm.score') {
			#######################MutSig1.0#######################################
			# MutSig1.0 statistic, defined in supplementary of MutSigCV in page 16.
			# Sg = Σc [–log10 binomial(nc, Nc, μc)
			# Method 1:
			#probs = dbinom(q[contexts], size[contexts], BMR[contexts], log=FALSE)
			#S.g = sum(-1 * log(probs, 10))
	
			# Method 2:
			# log10(p) = log(p)/log(10)
			probs.log = dbinom(q[contexts], size[contexts], BMR[contexts], log=TRUE)
			#probs.log10 = probs.log / log(10)
			#S.g = sum(-1 * probs.log10)
			## In this function, I use S.g = sum(-1 * probs.log), instead of the original
			##+ one (S.g = sum(-1 * probs.log10)) for ease of computation.
			S.g = sum(-1 * probs.log)
			perm[] <- 0
			for (context in contexts) {
				#perm.score = perm.score + rbinom(N, size[context], BMR[context])
				perm.prob.log = dbinom(rbinom(N, size[context], BMR[context]), size[context], BMR[context], log=TRUE)
				perm = perm + (-1 * perm.prob.log)
			}
			p= sum(perm >= S.g) / N
			#######################################################################
		} else if (method == 'perm.num') {
			n0 = sum(q)
			for (context in contexts) {
				perm = perm + rbinom(N, size[context], BMR[context])	
			}
			p= sum(perm >= n0) / N
		} else {
			cat("Unknown method name provided - ", method, file=stderr())
			cat("\n", file=stderr())
			q(save="no", status=1)
		}
	
		if (trace) {
			cat(paste("INFO  - ", date()), gene, q[contexts], sum(q), p, sep="\t")
			cat("\n", sep="")
		}
		out.df <- rbind( out.df, c(gene, q[contexts], sum(q), p, 1) )
	}
	
	out.df <- as.data.frame(out.df, stringsAsFactors=FALSE)
	out.df$p<- as.double(out.df$p)
	out.df$q<- p.adjust(out.df$p, method="fdr")
	out.df <- out.df[order(out.df$q), ]
	if (!is.na(output.file)) {
		write.table(out.df, output.file, col.names=out.df.names, row.names=FALSE, quote=FALSE, sep="\t")
	}
	invisible(out.df)
}

#mutsig(categs.file, coverage.file, mutations.file, output.file="out.PCT", method="PCT",trace=T, exclude.noncoding=F, sep="\t")
#mutsig(categs.file, coverage.file, mutations.file, output.file="out.FCPT", method="FCPT",trace=T, exclude.noncoding=F, sep="\t")
#mutsig(categs.file, coverage.file, mutations.file, output.file="out.LRT", method="LRT",trace=T,exclude.noncoding=F, sep="\t")
#mutsig(categs.file, coverage.file, mutations.file, output.file="out.perm.score", method="perm.score",trace=T, exclude.noncoding=F, sep="\t")
#mutsig(categs.file, coverage.file, mutations.file, output.file="out.perm.num", method="perm.num",trace=T, exclude.noncoding=F, sep="\t")
