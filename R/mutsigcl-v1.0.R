
## maf: an object returned by read.maf(...), read.table(...), or data.table::fread(...)
makeVRangesFromMaf <- function(maf)
{
  vr <- VRanges(maf@data$Chromosome, IRanges(maf@data$Start_Position,maf@data$End_Position), 
                ref = maf@data$Reference_Allele, alt=maf@data$Tumor_Seq_Allele2, 
                sampleNames = maf@data$Tumor_Sample_Barcode,
                totalDepth=maf@data$t_ref_count + maf@data$t_alt_count,
                refDepth = maf@data$t_ref_count, altDepth = maf@data$t_alt_count,
                Variant_Type = as.character(maf@data$Variant_Type),
                Variant_Classification = as.character(maf@data$Variant_Classification),
                Hugo_Symbol = as.character(maf@data$Hugo_Symbol),
                Entrez_Gene_Id = maf@data$Entrez_Gene_Id)
  return(vr)
}

annotateVRanges <- function(vr) {
  if (!exists("Hsapiens"))
    library(BSgenome.Hsapiens.UCSC.hg19)
  if (!exists("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  vr <- predictCoding(vr, TxDb.Hsapiens.UCSC.hg19.knownGene, Hsapiens, varAllele = DNAStringSet(vr@alt))
  return(unique(vr)) ## I don't know why predictCoding output duplicated annotations??
}

annotate.dbnsfp <- function(vr, dbnsfp.fl, verbose=FALSE) {
  #f <- open.TabixFile(dbnsfp.fl)
  tbx <- TabixFile(dbnsfp.fl)
  h <- headerTabix(tbx)
  tbx.seqnames <- h$seqnames
  tbx.headers <- strsplit(h$header, "\t")[[1]]
  ref_i = which(tbx.headers=="ref")
  alt_i = which(tbx.headers=="alt")
  
  SIFT_score_i = which(tbx.headers=="SIFT_score")
  Polyphen2_HDIV_score_i = which(tbx.headers=="Polyphen2_HDIV_score")
  CADD_raw_i = which(tbx.headers=="CADD_raw")
  
  vr <- vr[start(vr) == end(vr)] ## remove indels
  
  gr <- GRanges(sub("^chr","",seqnames(vr)), ranges(vr))
  #gr <- unique(gr)
  message(sprintf("Number of unique GRanges = %d.", length(gr)))
  
  #a <- unlist(scanTabix(tbx, param=gr))
  #d <- lapply(a, function(s) strsplit(s,"\t")[[1]])
  #d <- do.call("rbind", d)
  #colnames(d) <- tbx.headers
  
  # call scanTabix to scan dbnsfp database for every mutation
  d <- lapply(1:length(vr), function(k) {
    e <- vr[k]
    r <- rep(NA, length(tbx.headers))
    #a <- scanTabix(tbx, param = GRanges(sub("^chr","",seqnames(e)), ranges(e)))[[1]]
    # Or
    a <- unlist(scanTabix(tbx, param = GRanges(sub("^chr","",seqnames(e)), ranges(e))))
    if (length(a) < 1)
      return(r)
    for (i in 1:length(a)) {
      b <- strsplit(a[i],"\t")[[1]]
      p <- as.numeric(b[2])
      if (b[ref_i] == ref(e) && b[alt_i] == alt(e) && start(e) == p) {
        r <- b
        break
      }
    }
    if (verbose && k >= 100 && k %% 100 == 0)
      message(sprintf(" Processed %d mutations.", k))
    return(r)
    })
  d <- do.call("rbind",d)
  #d <- as.data.frame(d)
  #colnames(d) <- tbx.headers
  #####################
  #close(f)
  #d <- d[,c(SIFT_score_i, Polyphen2_HDIV_score_i, CADD_raw_i)]
  vr$SIFT_score <- as.numeric(as.character(d[, SIFT_score_i]))
  vr$Polyphen2_HDIV_score <- as.numeric(as.character(d[, Polyphen2_HDIV_score_i]))
  vr$CADD_raw <- as.numeric(as.character(d[, CADD_raw_i]))
  return(vr)
}

## maf: object returned by read.maf(...), read.table(...) or fread(...)
## vr: VRanges object required by SomaticSignatures::mutationContext(...)

collapseMutationContext <- function(maf=NULL, vr=NULL, removeSilent=FALSE) {
  if (is.null(vr)) {
    vr <- makeVRangesFromMaf(maf)
    vr <- vr[vr$Variant_Type=="SNP"]
  }
  if (removeSilent) {
    vr <- annotateVRanges(vr)
    vr <- vr[vr$CONSEQUENCE %in% c("nonsynonymous", "nonsense", "frameshift")]
  }
  
  vr <- mutationContext(vr, BSgenome.Hsapiens.UCSC.hg19)
  context <- paste(vr$alteration, vr$context, sep=":")
  ## Mutation signatures identified from TCGA nonhypermutated samples
  msig1 <- c("CA:C.A")
  msig2 <- c("CT:T.C")
  CtoTinCpG <- c("CT:A.G","CT:C.G","CT:G.G","CT:T.G")
  APOBEC <- c("CG:T.A","CG:T.T","CT:T.A","CT:T.T") # -> TCW[A/T]
  ## Kidney and utrothelial carcinomas specific signatures, especially in Asian countries.
  Aristolochic_Acid_Signature <- c("TA:C.G")
  hyper.msig19 <- c("CA:T.T","CT:T.C")
  others <- setdiff(context, c(msig1,msig2,CtoTinCpG,APOBEC,Aristolochic_Acid_Signature,hyper.msig19))
  
  collapseContext <- rep("others", length(vr))
  collapseContext[context %in% msig1] <- "Signature1"
  collapseContext[context %in% msig2] <- "Signature2"
  collapseContext[context %in% CtoTinCpG] <- "*CpG"
  collapseContext[context %in% APOBEC] <- "APOBEC"
  collapseContext[context %in% Aristolochic_Acid_Signature] <- "aristolochic acid"
  collapseContext[context %in% hyper.msig19] <- "hyper.Signature19"
  vr$collapseContext <- collapseContext
  return(vr)
}

## vr: object returned by collapseMutationContext
getObsGeneSigature <- function(vr, HugoSymbol=NULL, EntrezGeneId=NULL) {
  if (!is.null(EntrezGeneId)) {
    gene.vr <- subset(vr, Entrez_Gene_Id == EntrezGeneId)
  }
  else if (!is.null(HugoSymbol)) {
    gene.vr <- subset(vr, Hugo_Symbol == HugoSymbol)
  } else {
    stop("Missing Hugo_Symbol or Entrez_Gene_Id.")
  }
  d <- data.frame(start=start(gene.vr), collapseContext=gene.vr$collapseContext)
  return(d)
}

## maf <- read.maf(fn)
## vr <- collapseMutationContext(maf)
## obsSignature <- getObsGeneSigature(vr,"VHL")


getLongestCCDS <- function(d, Entrez_Gene_Id=NULL, Hugo_Symbol=NULL) 
{
  parse_cds_info <- function(a) {
    cds_locations <- a[10]
    chr <- as.character(a[1])
    strand <- as.character(a[7])
    a <- strsplit(substring(cds_locations, 2, nchar(cds_locations) - 1), ",")[[1]]
    a <- lapply(a, function(e) as.numeric(strsplit(e, "-")[[1]]))
    a <- do.call("rbind", a) + 1 ## index of NCBI CCDS file is zero based.
    if (!grepl("^chr", chr)) {
      chr <- paste("chr", chr, sep = "")
    }
    GRanges(chr, IRanges(start = a[,1], end = a[,2]), strand = strand)
  }
  #d <- subset(d, ccds_status == "Public")
  if (is.null(Hugo_Symbol))
    d <- subset(d, gene_id == Entrez_Gene_Id)
  else if (is.null(Entrez_Gene_Id))
    d <- subset(d, gene == Hugo_Symbol)
  if (nrow(d) == 0)
    return(NULL)
  
  gr <- parse_cds_info(d[1,])

  if (nrow(d) >= 2) {
    for (i in 2:nrow(d)) {
      gr1 <- parse_cds_info(d[i,])
      if (sum(width(ranges(gr1))) > sum(width(ranges(gr)))) {
        gr <- gr1
      }
    }
  }
  return(gr)
}

parse_cds_info <- function(a) {
  cds_locations <- a[10]
  chr <- as.character(a[1])
  strand <- as.character(a[7])
  a <- strsplit(substring(cds_locations, 2, nchar(cds_locations) - 1), ",")[[1]]
  a <- lapply(a, function(e) as.numeric(strsplit(e, "-")[[1]]))
  a <- do.call("rbind", a) + 1 ## index of NCBI CCDS file is zero based.
  if (!grepl("^chr", chr)) {
    chr <- paste("chr", chr, sep = "")
  }
  GRanges(chr, IRanges(start = a[,1], end = a[,2]), strand = strand)
}

getLongestIsoformPerPosNucleotieGRanges <- function(d, Entrez_Gene_Id=NULL, Hugo_Symbol=NULL) 
{
  #d <- subset(d, ccds_status == "Public")
  if (!is.null(Entrez_Gene_Id))
    d <- subset(d, gene_id == Entrez_Gene_Id)
  else if (!is.null(Hugo_Symbol))
    d <- subset(d, gene == Hugo_Symbol)
  if (nrow(d) == 0) {
    message(sprintf("geneid = %s not found in background reference genes.", Entrez_Gene_Id))
    return(NULL)
  }
  gr <- parse_cds_info(d[1,])
  
  if (nrow(d) >= 2) {
    for (i in 2:nrow(d)) {
      gr1 <- parse_cds_info(d[i,])
      if (sum(width(ranges(gr1))) > sum(width(ranges(gr)))) {
        gr <- gr1
      }
    }
  }
  
  ## Per position per nucleotie GRanges
  starts <- start(gr)
  ends <- end(gr)
  positions <- lapply(1:length(gr), function(i) starts[i]: ends[i])
  positions <- do.call("c", positions)
  GRanges(seqnames(gr)[1], IRanges(positions, positions))
}

## d: a data.frame refers to NCBI CCDS file with colname genes, e.g.
## CCDS_current_file="/Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2"
## d <- read.table(CCDS_current_file, comment.char = "", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

getBkgrGeneSignature <- function(d, HugoSymbol="VHL", EntrezGeneId=7428, removeSilent=FALSE) {
  mutTbl <- list(A=c("C","G","T"), C=c("A","G","T"), G=c("A","C","T"), T=c("A","C","G"))
  gr <- getLongestIsoformPerPosNucleotieGRanges(d, Hugo_Symbol = HugoSymbol, Entrez_Gene_Id = EntrezGeneId)
  if (is.null(gr))
    return(data.frame(start=NA, collapseContext=NA))
  
  seq <- getSeq(Hsapiens, gr)
  seq <- do.call("c", seq)
  s <- as.character(seq)
  bases <- sapply(1:nchar(s), function(i) substring(s,i,i))
  ref <- rep(bases, each=3)
  alt <- do.call("c",sapply(bases, function(b) mutTbl[b]))
  positions <- rep(start(gr), each=3)
  vr <- VRanges(seqnames(gr)[1], IRanges(positions, positions), ref = ref, alt = alt, sampleNames = "sampleId")
  cmc <- collapseMutationContext(vr=mutationContext(vr, Hsapiens), removeSilent = removeSilent)
  bkgrSignature <- data.frame(start=start(cmc), collapseContext=cmc$collapseContext)
  return(bkgrSignature)
}


## How to extract longest isoforms from TxDb.Hsapiens.UCSC.hg19.knownGene??
getBkgrGeneSignature_dev <- function(d, EntrezGeneId=7428, removeSilent=FALSE) {
  A <- c("C","G","T")
  C <- c("A","G","T")
  G <- c("A","C","T")
  T <- c("C","G","T")
  mutTbl <- list(A=c("C","G","T"), C=c("A","G","T"), G=c("A","C","T"), T=c("A","C","G"))
  
  #gr <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene, filter=list(gene_id=EntrezGeneId))
  ## if you wan to return gene_id and tx_name in gr:
  #gr <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene, c("gene_id","tx_name"), filter=list(gene_id=EntrezGeneId))
  #strand(gr) = "+" ## No reverse complement for sequence required.
  #chr <- as.character(seqnames(gr))[1]
  gr <- getLongestCCDS(d, Entrez_Gene_Id = EntrezGeneId)
  seqs <- getSeq(Hsapiens, gr)
  seq <- do.call("c", seqs)
  s <- as.character(seq)
  bases <- sapply(1:nchar(s), function(i) substring(s,i,i))
  ref <- rep(bases, each=3)
  alt <- do.call("c",sapply(bases, function(b) mutTbl[b]))
  positions <- rep(start(ranges(gr)), each=3)
  vr <- VRanges(Rle(chr), IRanges(positions, positions), ref = ref, alt = alt, sampleNames = "sampleId")
  cmc <- collapseMutationContext(vr=mutationContext(vr, Hsapiens), removeSilent = removeSilent)

  bkgrSignature <- data.frame(start=start(cmc), collapseContext=cmc$collapseContext)
  return(bkgrSignature)
}


getHotspotStatistic <- function(y, y.n=length(y)) {
  # A hotspot is defined as a 3-base-pair region of the gene containing many
  #+ mutations: at least 2, and at least 2% of the total mutations (nature12912).
  clust.table = try(table( cutree(fastcluster::hclust(dist(y), method='complete'), h=3)  ), silent=TRUE)
  statistic <- NA
  hotspot.num <- NA
  if (class(clust.table) != 'try-error') {
    statistic = sum(clust.table[clust.table >=2 & clust.table / sum(clust.table) >=0.02]) / y.n 
    hotspot.num = sum(clust.table >=2 & (clust.table / sum(clust.table) >=0.02))
  }
  return(c(statistic, hotspot.num))
}
getGrpId <- function(y, y.n=length(y)) {
  grp <- cutree(fastcluster::hclust(dist(y), method='complete'), h=3)
  return(grp)
}

mutsigcl_core <- function(bkgrSignature, obsSignature, global=FALSE, nsim=1000) {
  
  mutation_sampling <- function(obsSignature, bkgrSignature, global) {
    if (global) { ## Sampling without taking into account mutation context
      x <- bkgrSignature$start
      k <- nrow(obsSignature)
      y <- sample(x, k, replace = TRUE)
    } else { ## Sampling according to mutation context
      contexts <- table(obsSignature$collapseContext)
      y <- unlist(lapply(names(contexts), function(context) {
        x <- bkgrSignature$start[bkgrSignature$collapseContext == context]
        k <- contexts[context]
        sample(x, k, replace = TRUE)  
      }))
    }
    r <- getHotspotStatistic(y)
    return(r)
  }
  y0 <- obsSignature$start
  r <- getHotspotStatistic(y0)
  statistic0 <- r[1]
  hotspot.num0 <- r[2]
 
  r <- lapply(1:nsim, function(i) {
    mutation_sampling(obsSignature, bkgrSignature, global)
  })
  dat <- do.call("rbind", r)
  k <- sum(dat[, 1] >= statistic0)
  p.value <- (k + 1) / (nsim + 1)
  
  return(c(statistic0, hotspot.num0, p.value))
}

# maf.file: the filename of mutation file annotated by oncotator
# ccds.file: NCBI CCDS file
# maf: an object returned by read.maf(...)
# d: a data.frame referred to ccds.file
# k: the minimum mutation frequency of genes to be analyzed.
# removeSilent: If TRUE, silent mutations are excluded from analysis.
# nsim: number of simulations.
# mc.cores: number of cpu cores used.
mutsigcl <- function(maf.file, ccds.file, out.file=NULL, maf=NULL, d=NULL, removeSilent=FALSE, nsim=100, mc.cores=1) {
  ## loading required packages
  library(maftools)
  library(VariantAnnotation)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(SomaticSignatures)
  library(fastcluster)
  
  if (is.null(maf)) {
    maf <- read.maf(maf.file, removeSilent=removeSilent)
    message(sprintf("Finished reading %s.", maf.file))
  }
  if (is.null(d)) {
    d <- read.table(ccds.file, comment.char = "", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }
  
  ## Since read.maf(...) has already removed silent mutations, so there is no need to do
  ##+it here. But if you are serious about consistent annotation results, I recommend to
  ##+set removeSilent to TRUE and mutations will be re-annotated by VariantAnnotation
  ##+package which is used in background in silico mutation annotation.
  vr <- collapseMutationContext(maf, removeSilent = FALSE)
  
  x <- by(start(vr),vr$Entrez_Gene_Id, getHotspotStatistic)
  hotspot.num <- do.call("rbind", x)[,2]
  hotspotMutatedGeneIds <- names(hotspot.num)[!is.na(hotspot.num) & hotspot.num > 0]

  geneToMutatedSamples <- by(vr@sampleNames, vr$Entrez_Gene_Id, function(a) length(unique(a)))
  recurrentlyMutatedGeneIds <- names(geneToMutatedSamples)[geneToMutatedSamples>=3]
  hotspotMutatedGeneIds <- intersect(hotspotMutatedGeneIds, recurrentlyMutatedGeneIds)

  ###############--MUC16---TTN--
  long_genes <- c("94025","7273")
  #recurrentlyMutatedGeneIds <- setdiff(recurrentlyMutatedGeneIds, long_genes)
  hotspotMutatedGeneIds <- setdiff(hotspotMutatedGeneIds, long_genes)
  hotspotMutatedGeneIds <- hotspotMutatedGeneIds[hotspotMutatedGeneIds!="0"]
  hotspotMutatedGeneIds <- c("7157","7428")
  
  message(sprintf("Number of genes with hotspot mutations = %d.", length(hotspotMutatedGeneIds)))
  
  dat <- data.frame(Hugo_Symbol=maf@data$Hugo_Symbol, Entrez_Gene_Id=maf@data$Entrez_Gene_Id)
  dat <- apply(dat, 2, function(a) sub("\\s+","",as.character(a)))
  
  r <- mclapply(hotspotMutatedGeneIds, function(geneId, vr, d) {
    message(geneId)
    #library(BSgenome.Hsapiens.UCSC.hg19) ## For running in parallel, loading BSgenome in each worker is required. I have not tested this. NOT WORK!!!!!!
    bkgrSignature <- getBkgrGeneSignature(d, EntrezGeneId = geneId, removeSilent = removeSilent)
    if (any(is.na(bkgrSignature)))
      return(c(NA,NA,NA,NA))
    obsSignature <- getObsGeneSigature(vr, EntrezGeneId = geneId)
    res <- mutsigcl_core(bkgrSignature, obsSignature, global = FALSE, nsim = nsim)
    global.p.value <- mutsigcl_core(bkgrSignature, obsSignature, global = TRUE, nsim = nsim)[3]
    statistic <- res[1]
    hotspot.num <- res[2]
    p.value <- res[3]
    message(sprintf("## Entrez_Gene_Id = %s, statistic = %g, hotspot.num = %g, p.value = %g, global.p.value = %g", geneId, statistic, hotspot.num, p.value, global.p.value))
    return(c(statistic, hotspot.num, p.value, global.p.value))
  }, vr, d, mc.cores=mc.cores)
  Hugo_Symbols <- sapply(hotspotMutatedGeneIds, function(geneId) {
    dat[,1][dat[,2] == geneId][1]
  })
  r <- do.call("rbind", r)
  r <- as.data.frame(r)
  colnames(r) <- c("statistic", "hotspot.num", "p.value", "global.p.value")
  r <- cbind(Hugo_Symbol=Hugo_Symbols, Entrez_Gene_Id=names(Hugo_Symbols), r)
  my.df <- data.frame(lapply(r, as.character), stringsAsFactors=FALSE)
  my.df$q.value <- p.adjust(as.numeric(my.df$p.value), method="fdr")
  my.df$global.q.value <- p.adjust(as.numeric(my.df$global.p.value), method="fdr")
  my.df$MutatedSamples <- sapply(hotspotMutatedGeneIds, function(e) geneToMutatedSamples[[e]])
  my.df <- my.df[order(my.df$q.value, decreasing = FALSE), ]
  if (!is.null(out.file)) {
    write.table(my.df, file = out.file, quote=FALSE, sep = "\t", row.names = FALSE)
  }
  invisible(my.df)
}

maf.file="/Users/lixiangchun/Public/WorkSpace/Project/Precision_Medicine/Phase1_Mutations/Kidney_mutect2_PoN2_filtered.maf"
ccds.file="/Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2"
r=mutsigcl(maf=maf, d=d)
