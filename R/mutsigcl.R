
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

## maf: object returned by read.maf(...), read.table(...) or fread(...)
## vr: VRanges object required by SomaticSignatures::mutationContext(...)

collapseMutationContext <- function(maf=NULL, vr=NULL, removeSilent=FALSE) {
  if (is.null(vr)) {
    vr <- makeVRangesFromMaf(maf)
    vr <- vr[mcols(vr)$Variant_Type=="SNP"]
  }
  if (removeSilent) {
    vr <- annotateVRanges(vr)
    vr <- vr[mcols(vr)$CONSEQUENCE %in% c("nonsynonymous", "nonsense", "frameshift")]
  }
  #library(BSgenome.Hsapiens.UCSC.hg19)
  vr <- mutationContext(vr, Hsapiens)
  
  context <- paste(mcols(vr)$alteration, mcols(vr)$context, sep=":")
  ## Mutation signatures identified from TCGA nonhypermutated samples
  msig1 <- c("CA:C.A")
  msig2 <- c("CT:T.C")
  CtoTinCpG <- c("CT:A.G","CT:C.G","CT:G.G","CT:T.G")
  APOBEC <- c("CG:T.A","CG:T.T","CT:T.A","CT:T.T") # -> TCW[A/T]
  msig5 <- setdiff(context, c(APOBEC, CtoTinCpG, msig2, msig1))
  ## Kidney and utrothelial carcinomas specific signatures, especially in Asian countries.
  Aristolochic_Acid_Signature <- c("TA:C.G")
  
  collapseContext <- rep("Signature5", length(vr))
  collapseContext[context %in% msig1] <- "Signature1"
  collapseContext[context %in% msig2] <- "Signature2"
  collapseContext[context %in% CtoTinCpG] <- "*CpG"
  collapseContext[context %in% APOBEC] <- "APOBEC"
  collapseContext[context %in% Aristolochic_Acid_Signature] <- "aristolochic acid"
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
  d <- data.frame(start=start(ranges(gene.vr)), collapseContext=mcols(gene.vr)$collapseContext)
  return(d)
}

## maf <- read.maf(fn)
## vr <- collapseMutationContext(maf)
## obsSignature <- getObsGeneSigature(vr,"VHL")

## d: a data.frame refers to NCBI CCDS file with colname genes, e.g.
## CCDS_current_file="/Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2"
## d <- read.table(CCDS_current_file, comment.char = "", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

getBkgrGeneSignature <- function(d, HugoSymbol="VHL", EntrezGeneId=7428, removeSilent=FALSE) {
  #cds_locations <- subset(d, gene=="VHL", c(1,"cds_locations"))[,1]
  if (!is.null(EntrezGeneId)) {
    a <- subset(d, gene_id==EntrezGeneId & ccds_status == "Public")
  } else if (!is.null(HugoSymbol)) {
    a <- subset(d, gene==HugoSymbol & ccds_status == "Public")
  } else {
    stop("Either Hugo_Symbol or Entrez_Gene_Id must be provided.")
  }
  if (nrow(a) == 0) {
    message(sprintf("geneid = %s not found in background reference genes.", EntrezGeneId))
    return(data.frame(start=NA, collapseContext=NA))
  }
  
  cds_locations <- a[,10]
  chrs <- a[,1]
  
  if (all(grepl("chr",chrs) == FALSE)) {
    chrs <- paste("chr",chrs,sep="")
  }
  r <- strsplit(gsub("\\[|\\]", "", cds_locations), ",")
  
  ## get the longest isoform index
  longest_isoform_index <- 1
  if (length(r) > 1) {
    message(sprintf("More than 1 transcripts for %s, identifying the longest one.", HugoSymbol))
    longest_isoform_index <- try(
      which.max(sapply(1:length(r), function(i) {
      temp <- do.call("rbind", strsplit(r[[i]], "-"))
      temp <- apply(temp, 2, as.numeric)
      sum(temp[, 2] - temp[, 1])
    })), silent = FALSE)
    if (class(longest_isoform_index) == "try-error")
      print(r)
  }
  
  A <- c("C","G","T")
  C <- c("A","G","T")
  G <- c("A","C","T")
  T <- c("C","G","T")
  mutTbl <- list(A=c("C","G","T"), C=c("A","G","T"), G=c("A","C","T"), T=c("A","C","G"))

  #for (i in 1:length(r)) {
  #a <- r[[i]]
  a <- r[[longest_isoform_index]]
  b <- do.call("rbind",strsplit(a,"-"))
  cds_starts <- as.numeric(b[,1]) + 1
  cds_ends <- as.numeric(b[,2]) + 1
  
  grs <- lapply(1:length(cds_starts), function(j) {
    p <- cds_starts[j]: cds_ends[j]
    makeGRangesFromDataFrame(data.frame(chr=chrs[longest_isoform_index], start=p, end=p))
  })
  
  gr <- do.call("c", grs)
  seq <- getSeq(Hsapiens, gr)
  seq <- do.call("c", seq)
  s <- as.character(seq)
  bases <- sapply(1:nchar(s), function(i) substring(s,i,i))
  ref <- rep(bases, each=3)
  alt <- do.call("c",sapply(bases, function(b) mutTbl[b]))
  positions <- rep(start(ranges(gr)), each=3)
  vr <- VRanges(chrs[longest_isoform_index], IRanges(positions, positions), ref = ref, alt = alt, sampleNames = "sampleId")
  cmc <- collapseMutationContext(vr=mutationContext(vr, Hsapiens), removeSilent = removeSilent)
  #}
  bkgrSignature <- data.frame(start=start(ranges(cmc)), collapseContext=mcols(cmc)$collapseContext)
  return(bkgrSignature)
}

mutsigcl_core <- function(bkgrSignature, obsSignature, nsim=1000) {
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
  
  mutation_sampling <- function(obsSignature, bkgrSignature) {
    contexts <- table(obsSignature$collapseContext)
    y <- unlist(lapply(names(contexts), function(context) {
      x <- bkgrSignature$start[bkgrSignature$collapseContext == context]
      k <- contexts[context]
      sample(x, k, replace = TRUE)  
    }))
    r <- getHotspotStatistic(y)
    return(r)
  }
  y0 <- obsSignature$start
  r <- getHotspotStatistic(y0)
  statistic0 <- r[1]
  hotspot.num0 <- r[2]
  
  r <- lapply(1:nsim, function(i) {
    mutation_sampling(obsSignature, bkgrSignature)
  })
  dat <- do.call("rbind", r)
  k <- sum(dat[,1] >= statistic0)
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
mutsigcl <- function(maf.file, ccds.file, out.file, maf=NULL, d=NULL, k=2, removeSilent=FALSE, nsim=100, mc.cores=1) {
  ## loading required packages
  library(maftools)
  library(VariantAnnotation)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(SomaticSignatures)
  library(fastcluster)
  
  if (is.null(maf)) {
    maf <- read.maf(maf.file, removeSilent=removeSilent)
    message(sprintf("Finished reading %s.", maf.file))
  }
  if (is.null(d)) {
    d <- read.table(ccds.file, comment.char = "", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    message(sprintf("Finished reading %s.", ccds.file))
  }
  if (removeSilent) {
    message("Input param removeSilent=TRUE. Mutations will be re-annotated with additional time.")
  }
  vr <- collapseMutationContext(maf, removeSilent = removeSilent)
  
  geneToMutatedSamples <- by(vr@sampleNames, mcols(vr)$Entrez_Gene_Id, function(a) length(unique(a)))
  
  Entrez_Gene_Ids <- maf@data$Entrez_Gene_Id
  r <- table(Entrez_Gene_Ids)
  #recurrentlyMutatedGeneIds <- names(r)[r>2]
  #message(sprintf("Number of genes with >= 3 mutations is %d.", length(recurrentlyMutatedGeneIds)))
  recurrentlyMutatedGeneIds <- names(geneToMutatedSamples)[geneToMutatedSamples>=k]
  message(sprintf("Genes mutated in >= %d samples is %d.", k, length(recurrentlyMutatedGeneIds)))
  
  #recurrentlyMutatedGeneIds <- c("7157","7428")
  dat <- data.frame(Hugo_Symbol=maf@data$Hugo_Symbol, Entrez_Gene_Id=maf@data$Entrez_Gene_Id)
  dat <- apply(dat, 2, function(a) sub("\\s+","",as.character(a)))
  
  r <- mclapply(recurrentlyMutatedGeneIds, function(geneId, vr, d) {
    bkgrSignature <- getBkgrGeneSignature(d, HugoSymbol = NULL, EntrezGeneId = geneId, removeSilent = removeSilent)
    if (any(is.na(bkgrSignature)))
      return(c(NA,NA,NA))
    obsSignature <- getObsGeneSigature(vr, EntrezGeneId = geneId)
    res <- mutsigcl_core(bkgrSignature, obsSignature, nsim=nsim)
    message(sprintf("## Entrez_Gene_Id = %s, statistic = %g, hotspot.num = %g, p.value = %g.", geneId, res[1],res[2],res[3]))
    return(res)
  }, vr, d, mc.cores=mc.cores)
  Hugo_Symbols <- sapply(recurrentlyMutatedGeneIds, function(geneId) {
    dat[,1][dat[,2] == geneId][1]
  })
  r <- do.call("rbind", r)
  r <- as.data.frame(r)
  colnames(r) <- c("statistic", "hotspot.num", "p.value")
  r <- cbind(Hugo_Symbol=Hugo_Symbols, Entrez_Gene_Id=names(Hugo_Symbols), r)
  my.df <- data.frame(lapply(r, as.character), stringsAsFactors=FALSE)
  my.df$q.value <- p.adjust(as.numeric(my.df[,5]), method="fdr")
  my.df$MutatedSamples <- sapply(recurrentlyMutatedGeneIds, function(e) geneToMutatedSamples[[e]])
  my.df <- my.df[order(my.df$q.value, decreasing = FALSE), ]
  write.table(my.df, file = out.file, quote=FALSE, sep = "\t", row.names = FALSE)
  invisible(my.df)
}

maf.file="/Users/lixiangchun/Public/WorkSpace/Project/Precision_Medicine/Phase1_Mutations/Kidney_mutect2_PoN2_filtered.maf"
ccds.file="/Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2"
r=mutsigcl(maf.file, ccds.file,out.file = "out.txt")
