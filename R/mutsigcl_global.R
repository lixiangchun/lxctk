

generate_mutsigcl_global_bkgr_info <- function(CCDS_current_file)
{
	ccds_loc_parser <- function(ccds_location) {
		a <- strsplit(substring(ccds_location, 2, nchar(ccds_location)-1),',')[[1]]
		a <- lapply(a, function(e) as.numeric(strsplit(e,'-')[[1]]))
		a = lapply(1:length(a), function(i, a) seq(a[[i]][1]+1,a[[i]][2]+1), a)
		do.call('c',a)
	}
	d <- read.table(CCDS_current_file,comment.char="",header=TRUE,sep="\t",stringsAsFactors=FALSE)
	d <- subset(d, ccds_status == 'Public')
	geneid2cds_locations <- list()	
	gene_ids <- as.character(d$gene_id)
	for (i in 1:nrow(d)) {
		gene_id <- gene_ids[i]
		cds_locations <- geneid2cds_locations[[gene_id]]
		if (is.null(cds_locations)) {
			cds_locations <- ccds_loc_parser(d$cds_locations[i])
			geneid2cds_locations[[gene_id]] <- cds_locations
		} else {
			cds_locations1 <- ccds_loc_parser(d$cds_locations[i])
			if (length(cds_locations1) > length(cds_locations)) {
				geneid2cds_locations[[gene_id]] <- cds_locations1
			}
		}
	}
	geneid2cds_locations
}

## x=generate_mutsigcl_global_bkgr_info('CCDS.current.txt.bz2')

## Faster than sample(...)
sample_with_replace <- function(x, k, replace=TRUE) x[.Internal(sample(length(x), k, replace, NULL))]

mutsigcl_global_core <- function(mutated_positions, bkgr_mutated_positions, nperm, gene_id, Hugo_Symbol)
{
	compute_hotspot_statistic <- function(y) { # y: mutated positions
		statistic <- -1 
		hotspot_num <- NA
		# A hotspot is defined as a 3-base-pair region of the gene containing many
		#+ mutations: at least 2, and at least 2% of the total mutations (nature12912).
		if (length(y) > 1) {
			clust.table = try(table( cutree(fastcluster::hclust(dist(y), method='complete'), h=3)  ), silent=TRUE)
			if (class(clust.table) != 'try-error') {
				statistic = sum(clust.table[clust.table >=2 & clust.table / sum(clust.table) >=0.02]) / length(y)
				hotspot_num <- sum(clust.table >=2 & (clust.table / sum(clust.table) >=0.02))
			}   
		}
		list(statistic=statistic, hotspot_num=hotspot_num)
	}
	r <- compute_hotspot_statistic(mutated_positions)
	obs_hotspot_statistic <- r$statistic

	n <- length(mutated_positions)
	simu_hotspot_statistics <- sapply(1:nperm, function(i) compute_hotspot_statistic(sample_with_replace(bkgr_mutated_positions, n))$statistic)
	k <- sum(simu_hotspot_statistics >= obs_hotspot_statistic)
	
	#list(hotspot_statistic=obs_hotspot_statistic, k=k, n=nperm, p.value=k/nperm)
	p = (k + 1) / (nperm + 1)
	x <- c(gene_id, Hugo_Symbol, obs_hotspot_statistic, r$hotspot_num, k, nperm, p)
	x
}

mutsigcl_global <- function(Oncotator_file, NCBI_CCDS_GeneID2CDS_Positions_Rdata_file, d=NULL, NCBI_CCDS_GeneID2CDS_Positions=NULL, nperm=20000, outfile='mutsigcl_global_output.txt', exclude_indel=TRUE)
{
	if (is.null(NCBI_CCDS_GeneID2CDS_Positions)) {
		load(NCBI_CCDS_GeneID2CDS_Positions_Rdata_file)
		message(">>>>> Finished loading NCBI_CCDS_GeneID2CDS_Positions_Rdata_file. <<<<<<")
	}
	if (is.null(d)) {
		d <- data.table::fread(Oncotator_file)		
		data.table::setDF(d)
		d <- d[,c('Hugo_Symbol','Entrez_Gene_Id','Start_position','Variant_Classification','Variant_Type')]
		#d <- d[!grepl("Silent|IGR|Intron|UTR|lincRNA|RNA", d$Variant_Classification),]
		d <- d[!grepl("Silent|IGR|Intron", d$Variant_Classification),]
	}
	if (exclude_indel) {
		d <- d[d$Variant_Type == 'SNP',]
	}

	d$Entrez_Gene_Id <- as.character(d$Entrez_Gene_Id)
	message(">>>>>> Finished reading mutation data from Oncotator_file <<<<<<.")

	gene_ids <- unique(d$Entrez_Gene_Id)	
	Start_positions <- d$Start_position

	Hugo_Symbols <- d$Hugo_Symbol
	names(Hugo_Symbols) = d$Entrez_Gene_Id
	
	obs_geneid2mutated_positions <- by(d$Start_position, d$Entrez_Gene_Id, c)
	recurrent_mutated_gene_ids <- c()
	for (gene_id in names(obs_geneid2mutated_positions)) {
		if (length(obs_geneid2mutated_positions[[gene_id]]) >= 3 && !is.null(NCBI_CCDS_GeneID2CDS_Positions[[gene_id]])) {
			recurrent_mutated_gene_ids <- c(recurrent_mutated_gene_ids, gene_id)
		}
	}
	message(sprintf(">>>>>> Recurrently mutated gene number = %d <<<<<<\n", length(recurrent_mutated_gene_ids)))

	results <- mclapply(recurrent_mutated_gene_ids, function(gene_id, obs_geneid2mutated_positions, NCBI_CCDS_GeneID2CDS_Positions, Hugo_Symbols) {
		mutated_positions <- obs_geneid2mutated_positions[[gene_id]]
		bkgr_mutated_positions <- NCBI_CCDS_GeneID2CDS_Positions[[gene_id]]
		r <- rep(NA, 6)
		if (!is.null(mutated_positions) && !is.null(bkgr_mutated_positions)) {
			r <- mutsigcl_global_core(mutated_positions, bkgr_mutated_positions, nperm, gene_id, Hugo_Symbols[gene_id])
			cat(sprintf("PID: %s", Sys.getpid()), gene_id, Hugo_Symbols[gene_id], r)
			cat("\n")
		}
		r
	}, obs_geneid2mutated_positions, NCBI_CCDS_GeneID2CDS_Positions, Hugo_Symbols)

	results <- do.call("rbind", results)
	results <- na.omit(as.data.frame(results))
	my.df <- data.frame(lapply(results, as.character), stringsAsFactors=FALSE)
	names(my.df) = c("gene_id",'Hugo_Symbol','hotspot_statistic','hotspot_num','extreme_nperm', 'permutation_nperm', 'p.value')
	my.df$q.value <- p.adjust(as.numeric(my.df$p.value, method='fdr'))
	my.df <- my.df[order(my.df$q.value, decreasing=FALSE),]
	write.table(my.df, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
}

#library(fastcluster)
#library(data.table)
#library(parallel)
#options(mc.cores=8)
#mutsigcl_global('../MutSigCV/exomic.maf','/ifshk1/BC_CANCER/03user/lixiangchun/prj/BGI_Data/Precision_Medicine/Kidney/MutSigCL_Global/NCBI_CCDS_GeneID2CDS_Positions.RData', 20000, 'MutSigCL_global_output2.txt')

