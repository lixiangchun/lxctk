
# bamFileLists: bamfile list
# fastaFile: bwa-indexed reference sequence
# OncotatorFile: Oncotator annotated file
# outfile: output file
# lines: process # lines each time. Reduce this number if too many bam files provided as input
# flag: If TRUE return a data.frame, may be very large if there huge number of mutations
# minBQ: minBaseQuality
# minMQ: minMapQuality
# minDepth: allelic frequency is set to 0 if sequencing depth at the mutation site is < minDepth

calc_allelic_freqs <- function(bamFileLists, fastaFile, OncotatorFile, outfile, lines=80, flag=FALSE, minBQ=0, minMQ=13, minDepth=8)
{
	fls <- PileupFiles(bamFileLists)
	fa <- open(FaFile(fastaFile))

	d <- try(data.table::fread(OncotatorFile), silent=TRUE)
	if (class(d) == 'try-error') {
		d <- read.table(OncotatorFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
		warning(sprintf("Exception in calling `data.table` to read `%s`.", OncotatorFile))
	} else {
		data.table::setDF(d)
	}

	d <- subset(d, Variant_Type=="SNP")
	if (any(grepl('^chr',seqnames(seqinfo(fa)))) && !any(grepl('^chr', d$Chromosome))) {
		d$Chromosome <- paste('chr', d$Chromosome, sep="")
	}
	mafFields <- c('Chromosome','Start_position','End_position','Reference_Allele','Tumor_Seq_Allele2')
	intervals <- d[,mafFields]
	colnames(intervals) <- c('chr','start','end','ref','alt')

	message("Finished initializing inputs.")
	
	x <- list()
	idx <- split(1:nrow(d), ceiling(seq_along(1:nrow(d))/lines))
	n <- length(idx)
	for (i in 1:n) {
		r <- calc_allelic_freqs_core(fls=fls, fa=fa, intervals=intervals[idx[[i]],], minBQ=minBQ, minMQ=minMQ, minDepth=minDepth)	
		xi <- cbind(d[idx[[i]],], r)
		if (flag) {
			x[[i]] <- xi
		}
		# If there is error: unimplemented type 'list' in 'EncodeElement' (often occurs in data.frame returned by do.call(...)), 
		#+use the following clause to process xi:
		# xi <- data.frame(lapply(xi, as.character), stringsAsFactors=FALSE)
		if (i == 1) {
			write.table(xi, file=outfile, quote=FALSE, sep="\t")
		} else {
			write.table(xi, file=outfile, quote=FALSE, sep="\t", col.names=FALSE, append=TRUE)
		}
		message(sprintf("Finished processing chunk `%s/%s`.", i, n))
	}
	close(fa)
	if (flag) {
		x <- do.call("rbind", x)
		x <- data.frame(lapply(x, as.character), stringsAsFactors=FALSE)
		return(x)
	}
}

# bamFileLists: an array of bam files
# fls: object returned by PileupFiles
# fastaFile: path to indexed reference
# fa: object returned by FaFile
# OncotatorFile: Oncotator annotated somatic mutations
# intervals: a data frame refers to somatic mutations. The data frame must have 5 columns and the column names must be
#+ c('chr','start','end','ref','alt'). Note that only SNV is supported.
# minBQ: minBaseQuality
# minMQ: minMapQuality
# minDepth: allelic frequency is set to 0 if sequencing depth at the mutation site is < minDepth

calc_allelic_freqs_core <- function(bamFileLists, fls, fastaFile, fa, OncotatorFile=NULL, intervals=NULL, minBQ=0, minMQ=13, minDepth=8)
{
	#calcInfo <- function(x, ref_base, alt_base) {
	#	x1 <- x[[1]]$seq[,,]
	#	n_ref_counts <- x1[ref_base,]	
	#	n_alt_counts <- x1[alt_base,]
	#	af <- n_alt_counts / (n_alt_counts + n_ref_counts)
	#	af
	#}
	if (is.null(OncotatorFile) && is.null(intervals))
		stop("At least `OncotatorFile` or `intervals` must be provided.")
	
	calcInfo <- function(x) x

	if (is.null(intervals)) {
		d <- data.table::fread(OncotatorFile)
		data.table::setDF(d)
		d <- subset(d, Variant_Type=="SNP")
		mafFields <- c('Chromosome','Start_position','End_position','Reference_Allele','Tumor_Seq_Allele2')
		intervals <- d[,mafFields]
		colnames(intervals) <- c('chr','start','end','ref','alt')
		message(sprintf(">>>>>> Finished reading '%s' <<<<<<.", OncotatorFile))
	}
	idx <- (intervals$start != intervals$end)	
	if (any(idx)) {
		warning("Only SNV supported, the other variants is discarded.")
		intervals <- intervals[!idx,]
	}

	cat(">>>>>> Interval size:", dim(intervals))
	cat("\n")

	if (missing(fls)) {
		fls <- PileupFiles(bamFileLists)	
	}
	which <- GenomicRanges::makeGRangesFromDataFrame(intervals)
	param <- ApplyPileupsParam(which=which, yieldBy="position", what="seq", maxDepth=1e6,minBaseQuality=minBQ,minMapQuality=minMQ)
	pls <- applyPileups(fls, calcInfo, param=param)
	cat(">>>>>> Finished pileup for all sites. <<<<<<\n")

	if (missing(fa)) {
		fa <- open(FaFile(fastaFile)) 
	}
	refseqs <- toupper(scanFa(fa, param=which))
	if (missing(fa)) {
		close(fa)
	}

	cat(">>>>>> Finished extracting refseqs. <<<<<<\n")
	
	allelic_freqs <- 
	lapply(1:nrow(intervals), function(i, pls, intervals) {
		alt_base <- intervals$alt[i]
		#ref_base <- intervals$ref[i]
		n_alt_counts <- pls[[i]]$seq[,,1][alt_base,]
		#n_ref_counts <- pls[[i]]$seq[,,1][ref_base,]
		#depths <- colSums(pls[[1]]$seq[,,1][c(ref_base, alt_base),])
		depths <- colSums(pls[[i]]$seq[,,1])
		afs <- n_alt_counts / depths
		afs[depths < minDepth] <- 0 
		afs
	}, pls, intervals)
	#z <- cbind(d, do.call("rbind", allelic_freqs))
	#write.table(z, file=outfile, quote=FALSE, sep="\t")
	#list(pls=pls, refseqs=refseqs, afs=allelic_freqs)
	#cbind(intervals, do.call("rbind", allelic_freqs))
	do.call("rbind", allelic_freqs)
}

#library(Rsamtools)
#bamFileLists=read.table("/ifs4/BC_CANCER/PROJECT/HUMfxfX/limiao/review/mutect/new.bam.list",header=FALSE,stringsAsFactors=FALSE)$V1
#fastaFile="/ifs4/BC_CANCER/PROJECT/HUMfxfX/limiao/review/hg19_fasta_GATK/hg19.fasta"
#fastaFile="/ifshk1/BC_CANCER/01bin/01bin/lixiangchun/db/hg19Virus/bwa-0.7.12/hg19Virus.fasta"
#OncotatorFile="CRC_WES.maf"

#bamFileLists=read.table("/ifshk5/BC_CANCER/PROJECT/HKC11086_HUMunaX/20150716_intelligent/limiao_paper/recapseg/bam.list",header=FALSE,stringsAsFactors=FALSE)$V1
#OncotatorFile="CRC_CapSeg.maf"
#calc_allelic_freqs(bamFileLists, fastaFile, OncotatorFile, outfile='CRC_CapSeg_Allelic_freqs.maf', lines=20)


