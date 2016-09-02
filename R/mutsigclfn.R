#args <- commandArgs(TRUE)
#if (length(args) != 5) {
#	message('\n\tUsage: pbiSMG <bkgrSQLiteDB> <obsInfile> <outfile> <nperm> <type/(FN/CL)>\n')
#	quit('no')
#}
#print(args)

#bkgrSQLiteDB = args[1]
#obsInfile = args[2]
#outfile = args[3]
#nperm = as.integer(args[4])
#type = args[5]

#library(parallel)
#library(RSQLite)
#library(fastcluster)      # use hclust in this package to accelerate clustering

is.SQLiteConnection <- function(db) {
	attributes(db)$class[1] == 'SQLiteConnection'
}

# This routine is a little bit faster than `sample` in R
sample_with_replace <- function(x, k, replace=TRUE) x[.Internal(sample(length(x), k, replace, NULL))]

# Extract a data frame for `gene` from `d`
# @d: a data frame of 3 columns, for example:
# @gene: a gene symbol
#SSPO	T->C	149475969,149477156,149512431
#SSPO    T->A	0
#SSPO    T->G	149484522,149493556,149493731
#SSPO    C->T	149477168,149482588,149482655,149482688,149485688
#SSPO    *TpA->mut	149477527
#SSPO    TCW->mut	149476736
#SSPO    C->A	0
get.gene.spectra <- function(d, gene) {
	out = NA
	if (is.SQLiteConnection(d)) {
		s = sprintf("SELECT * from bkgr_table where gene == '%s'", gene)	
		data = dbGetQuery(d, s)
	} else {
		data=d[d[,1]==gene,]
	}
	#data = data[!grepl("NULL", data[,3]),] # this can be done for obs_data in pbiSMG 
	if (nrow(data) != 0) {
		out=list()
		for (i in 1:nrow(data)) {
			spectrum=data[i,2]	
			out[[spectrum]]=as.numeric(strsplit(data[i,3],',')[[1]])
		}
	}
	out
}

# @dl: a list returned by get.gene.spectra
# @type: CL - MutSigCL, FN - MutSigFN
# @hotspot.alg: algorithm to define hotspot statistic, 'ratio' is much faster than 'hclust'
MutSigXX.Statistic <- function(dl, type=c('CL', 'FN'), hotspot.alg=c('hclust','ratio'), report.hotspot.num=FALSE) {
	#type <- match.arg(type)
	# match.arg takes extra time to parse input parameters. In a many permutations, aggregated time
	#+ is is not trivival. So I use the 2 following clauses to replace match.arg.
	type = type[1]
	hotspot.alg = hotspot.alg[1]
	hotspot.num = 0
	if (type == 'FN') {
		# 1: access list with name, i.e. spectra, or names(dl)
		y = sapply(names(dl), function(i, dl) mean(dl[[i]]), dl)
		
		# 2: access list with index, i.e. 1:length(dl)
		#y = sapply(1:length(dl), function(i, dl) mean(dl[[i]]), dl)

		y = y[!is.na(y)]
		statistic = mean(y)      # FN-statistic for conservative score or functional impact.
	} else {
		#y = stack(dl)$values    # This is very time consuming. y=unlist(dl) is much faster!!! 		
		y = unlist(dl)           
		#y = y[!is.na(y)]
		y.n = length(y)
		statistic = 0            # CL-statistic: the hotspot statistic (fraction).
		if (y.n > 1) {
			if (hotspot.alg == 'hclust') {
				# A hotspot is defined as a 3-base-pair region of the gene containing many
				#+ mutations: at least 2, and at least 2% of the total mutations (nature12912).
				clust.table = try(table( cutree(fastcluster::hclust(dist(y), method='complete'), h=3)  ), silent=TRUE)
				if (class(clust.table) != 'try-error') {
					statistic = sum(clust.table[clust.table >=2 & clust.table / sum(clust.table) >=0.02]) / y.n
					if (report.hotspot.num) {
						hotspot.num = sum(clust.table >=2 & (clust.table / sum(clust.table) >=0.02))
					}
				}
			} else {
				statistic = y.n / length(unique(y))	 - 1.0
			}
		}
	}
	list(statistic=statistic, hotspot.num=hotspot.num)
}

# Permutating for only 1 gene.
# `bkgr` and `obs`: Lists returned from get.gene.spectra(...)
# `type`: Analysis-type, CL - MutSigCL, FN - MutSigFN
# `hotspot.alg`: Algorithm to calculate CL-statistic.
# `nperm`: Number of permutations
# `mc.cores`: Number of cores used in parallel processing
one.gene.permute <- function(bkgr, obs, type=c('CL', 'FN'), hotspot.alg=c('hclust','ratio'), nperm=1e+4, mc.cores=4) {
	spectra=names(obs)
	# `i`: used as permutation index, being the 1st parameter in lapply(...)
	one.gene.permute.core <- function(i, bkgr, obs, spectra, type, hotspot.alg) {
		bkgr_perm = list()
		for (spectrum in spectra) {
			x0 = obs[[spectrum]]
			x0.n = length(x0)
			if (x0.n == 0) next;
			#if (is.null(x0) | any(is.na(x0))) next;    # indicates no mutation at this context
			x = bkgr[[spectrum]]
			#xp = sample(x, x0.n, TRUE)
			xp = sample_with_replace(x, x0.n, TRUE)
			bkgr_perm[[spectrum]] = xp
		}
		#statistic = MutSigXX.Statistic(bkgr_perm, type, hotspot.alg)
		MutSigXX.Statistic(bkgr_perm, type, hotspot.alg)
	}
	#statistic0 = MutSigXX.Statistic(obs, type, hotspot.alg)
	res = MutSigXX.Statistic(obs, type, hotspot.alg, report.hotspot.num=TRUE)
	statistic0 = res$statistic
	hotspot.num = res$hotspot.num

	#perm.statistics <- lapply(1:nperm, one.gene.permute.core, bkgr=bkgr, obs=obs, spectra=spectra, type=type, hotspot.alg=hotspot.alg)
	r <- lapply(1:nperm, one.gene.permute.core, bkgr=bkgr, obs=obs, spectra=spectra, type=type, hotspot.alg=hotspot.alg)
	# lapply returns an array of lists, concatenate all lists into a vector using unlist.
	#perm.statistics <- unlist(perm.statistics)
	perm.statistics <- as.numeric(do.call("rbind", r)[,1])

	counter = sum(perm.statistics >= statistic0)
	p=(counter + 1) / (nperm + 1)
	list(statistic=statistic0, hotspot.num=hotspot.num, p.value=p)
}

# `bkgrSWLiteDB`: An SQLite object to bkgr feature data
# `obs_data`: A data frame object to the observed feature data
# `outfile`: Output file
# `genes`: Specified genes to be analyzed, if NULL use all genes in obs_data.
# `type`: Analysis-type, CL - MutSigCL, FN - MutSigFN
# `hotspot.alg`: Algorithm to define hotspot statistic (or CL-statistic)
# `min.cl`: Min CL-statistic set for 2nd large permutations. Genes with ≥min.cl in the 1st MutSigCL were kept in the 2nd MutSigCL.
# `nperm`: Number of permutations
# `mc.cores`: Number of cores used in parallel processing
# `bkgr_data`: Object to bkgrSQLiteDB, do not change it always.
mutsigclfn <- function(bkgrSQLiteDB, obs_data, outfile='out.txt', genes=c(), type=c('CL','FN'), hotspot.alg=c('hclust','ratio'), min.cl=0.2, nperm=1000, mc.cores=4, bkgr_data=dbConnect(dbDriver('SQLite'), bkgrSQLiteDB)) {
	tck0=proc.time()[3]
	if (type[1] != "CL" &&  type[1] != "FN") {
		cat(sprintf('\nERROR - Type must be in c("CL","FN"), input type is "%s".\n\n', type[1]), file=stderr())
		quit('no')
	}
	if (hotspot.alg[1] != "hclust" &&  hotspot.alg[1] != "ratio") {
		stop(sprintf('\nERROR - hotspot.alg must be in c("hclust","ratio"), input type is "%s".\n\n', hotspot.alg[1]), file=stderr())
	}

	mutsigclfn_core <- function(gene, bkgr_data, obs_data, type, hotspot.alg, nperm=nperm, mc.cores=mc.cores) {
		bkgr = get.gene.spectra(bkgr_data, gene)
		obs = get.gene.spectra(obs_data, gene)
		if (is.list(bkgr) && is.list(obs)) {
			r = one.gene.permute(bkgr, obs, type, hotspot.alg, nperm, mc.cores)
			INFO = sprintf("- mutsigclfn_core gene - %s, statistic = %.3f, hotspot.num = %g, p-value = %g\n", gene, r$statistic, r$hotspot.num, r$p.value)
			#cat(as.character(Sys.time()), INFO, sep="", file=stderr())
			cat("PID:", Sys.getpid(), INFO, sep=" ", file=stderr())
			res = list(gene=gene, statistic=r$statistic, hotspot.num=r$hotspot.num, p.value=r$p.value)
		} else {
			res = list(gene=gene, statistic=NA, hotspot.num=NA, p.value=NA)
		}
		res
	}
	# Remove obs_data entries with NULL features; obs_data[,3] feature column
	obs_data = obs_data[!grepl("NULL", obs_data[,3]),]
	if (is.null(genes)) {
		genes = unique(as.character(obs_data[,1]))
	}
	if (type == 'CL' && min.cl > 0) {
		cat('INFO - mutsigclfn - To save time, mutsigclfn selects candidate genes to perform MutSigCL analysis.\n', file=stderr())
		ri = mclapply(genes, function(gene, obs_data, type, hotspot.alg) MutSigXX.Statistic(get.gene.spectra(obs_data, gene), type, hotspot.alg), obs_data, type, hotspot.alg, mc.cores=mc.cores)	
		genes = genes[ri>min.cl]
	}

	cat('Number of genes to be permutated:', length(genes), '\n', file=stderr())
	cat('INFO - mutsigclfn - Begin permutation in parallel.\n', file=stderr())
	# mclapply returns an array of lists.
	cl = mclapply(genes, mutsigclfn_core, bkgr_data, obs_data, type, hotspot.alg, nperm, mc.cores=mc.cores)

	###############Output to file########################
	cl.n <- length(cl)
	genes = rep(0, cl.n)
	statistics = rep(0, cl.n)
	hotspot.nums = rep(0, cl.n)
	p.values = rep(1, cl.n)
	for (i in 1:cl.n) {
		genes[i] = cl[[i]]$gene
		statistics[i] = cl[[i]]$statistic
		hotspot.nums[i] = cl[[i]]$hotspot.num
		p.values[i] = cl[[i]]$p.value
	}
	data = data.frame(gene=genes, statistic=statistics, hotspot.num=hotspot.nums, p.value=p.values)
	data = na.omit(data)
	data$q.value = p.adjust(data$p.value, method='fdr')    # BH-FDR adjusted
	data = data[order(data$q.value, decreasing=FALSE),]    # order data by q.value in ascending order.
	write.table(data, file=outfile, row.names=FALSE, quote=FALSE, sep="\t")
	tck=proc.time()[3]
	cat(as.character(Sys.time()), sprintf('- Total time elapsed: %s secs.\n', tck-tck0), file=stderr())
	if (is.SQLiteConnection(bkgr_data)) {
		dbDisconnect(bkgr_data)
	}
	invisible(data)
}

#obs_data=read.table(obsInfile, stringsAsFactors=FALSE,header=FALSE)
#genes=unique(obs_data[,1])

# 1: This will generate scheduled cores XXX encountered errors in user code, all values of the jobs will be affected
#sqlite <- dbDriver('SQLite')
#bkgr_data <- dbConnect(sqlite, 'bkgr.db')
#pbiSMG(bkgrSQLiteDB, obs_data, 'out.txt',genes, 'FN', 10, 2, dbConnect(sqlite, bkgrSQLiteDB))

# 2: The correct approach
#pbiSMG(bkgrSQLiteDB, obs_data, outfile, genes, type, nperm, 4, dbConnect(dbDriver('SQLite'), bkgrSQLiteDB))

