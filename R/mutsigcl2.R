

# This routine is a little bit faster than `sample` in R
faster.sample <- function(x, k, replace=TRUE) x[.Internal(sample(length(x), k, replace, NULL))]

## Combining multiple p-values. The source code was copied from combine.test(...) in package survcomp.
combine.p.values <- function (p, weight, method = c("fisher", "z.transform", "logit"), 
          hetero = FALSE, na.rm = TRUE) 
{
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


read.ncbi.ccds.file <- function(ccds.file) {
  d <- read.table(ccds.file, comment.char = "", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  d <- subset(d, ccds_status == "Public")
  return(d)
}

## maf: an object returned by read.maf(...), read.table(...), or data.table::fread(...)
makeVRangesFromMaf <- function(maf, removeIndel=TRUE)
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
  if (removeIndel) {
    vr <- vr[vr$Variant_Type=="SNP"]
  }
  return(vr)
}

makeVRangesFromOncotator <- function(maf.file, removeIndel = TRUE, removeSilent = TRUE, hyperCutoff=NULL) {
  d <- suppressWarnings(data.table::fread(maf.file, sep = "\t"))
  totalMutation <- table(d$Tumor_Sample_Barcode)
  Tumor_Sample_Barcodes <- names(totalMutation)
  message(sprintf("Total sample n = %d.", length(Tumor_Sample_Barcodes)))
  if (!is.null(hyperCutoff)) {
    message(sprintf("Here I use optimal univariate k-means clustering to identify hypermutated samples. Samples belonging to clusters with centers > %g are considered to be hypermutated.", hyperCutoff))
    r <- Ckmeans.1d.dp(totalMutation)
    hyperMutatedClustIds <- which(r$centers > hyperCutoff)
    hyperMutatedSamples <- names(totalMutation)[r$cluster %in% hyperMutatedClustIds]
    d <- subset(d, !(Tumor_Sample_Barcode %in% hyperMutatedSamples))
    message(sprintf("Removed hypermutated samples (n = %d):", length(hyperMutatedSamples)))
    print(totalMutation[hyperMutatedSamples], file=stderr())
  }
  if (removeIndel) {
    d <- subset(d, Variant_Type == "SNP")
  }
  # Silent mutatios defined in maftools::read.maf(...)
  silentMutations <- c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron","RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA")
  if (removeSilent) {
    d <- subset(d, !(Variant_Classification %in% silentMutations))
  }
  vr <- VRanges(d$Chromosome, IRanges(d$Start_position, d$End_position),
                ref = d$Reference_Allele, alt = d$Tumor_Seq_Allele2,
                sampleNames = d$Tumor_Sample_Barcode,
                ##totalDepth = d$t_ref_count + d$t_alt_count, refDepth = d$t_ref_count, altDepth = d$t_alt_count, ## Some maf.file don't have this columns.
                Variant_Type = as.character(d$Variant_Type),
                Variant_Classification = as.character(d$Variant_Classification),
                Hugo_Symbol = as.character(d$Hugo_Symbol),
                Entrez_Gene_Id = as.character(d$Entrez_Gene_Id))
  return(vr)
}

# Annotate a list of mutations in vr (VRanges) by VariantAnnotation::predictCoding(...).
# Currently, predictCoding(...) does not support indel annotation.
# vr: mutations in VRanges format.
annotateVRanges <- function(vr) {
  if (!exists("Hsapiens"))
    library(BSgenome.Hsapiens.UCSC.hg19)
  if (!exists("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  vr <- predictCoding(vr, TxDb.Hsapiens.UCSC.hg19.knownGene, Hsapiens, varAllele = DNAStringSet(vr@alt))
  return(unique(vr)) ## I don't know why predictCoding output duplicated annotations??
}

# Annotate mutations against dbNSFP.
# vr: mutations in VRanges format.
# dbnsfp.fl: the tabix-indexed filename of dbNSFP database. Note that there are multiple annotation
#+for a single mutation in dbSNFP, e.g. there are 2 SIFT prediction scores available for a single
#+mutation, which was separated by a comma (i.e. ';'). Users have to preprocess dbNSFP to make sure
#+that there is only one annotation for a single mutation. This can be easily achieved with Python:

###########---Python script to preprocess dbNSFP----############################
#import gzip
#fn="/home/lixiangchun/.work/database/oncotator_v1_ds_Jan262014/dbNSFP_ds/hg19/dbNSFP2.4_variant.tabix_indexed.tsv.gz"
#f=gzip.open(fn)
#h=f.readline().split("\t")
#idxs = [0,1,2,3,4,5,7]
#Uniprot_aapos_i = h.index("Uniprot_aapos")
#SIFT_score_i = h.index("SIFT_score")
#Polyphen2_HDIV_score_i = h.index("Polyphen2_HDIV_score")
#CADD_raw_i = h.index("CADD_raw")
#h_idx = []
#h_idx.extend(idxs)
#h_idx.extend([Uniprot_aapos_i,SIFT_score_i,Polyphen2_HDIV_score_i,CADD_raw_i])
#print "\t".join([h[i] for i in h_idx])

#for s in f:
#  a = s.split("\t")
#  b = [a[i] for i in idxs]
#  Uniprot_aapos = a[Uniprot_aapos_i].split(";")[0]
#  SIFT_score = a[SIFT_score_i].split(';')[0]
#  Polyphen2_HDIV_score = a[Polyphen2_HDIV_score_i].split(';')[0]
#  CADD_raw = a[CADD_raw_i].split(';')[0]
#  b.extend([Uniprot_aapos, SIFT_score, Polyphen2_HDIV_score, CADD_raw])
#  print "\t".join(b)
#  f.close()

annotate.dbnsfp <- function(vr, dbnsfp.fl, tbx=NULL, verbose=FALSE) {
  if (is.null(tbx)) {
    tbx <- open(TabixFile(dbnsfp.fl))
  }
  h <- headerTabix(tbx)
  tbx.seqnames <- h$seqnames
  tbx.headers <- strsplit(h$header, "\t")[[1]]
  ref_i = which(tbx.headers=="ref")
  alt_i = which(tbx.headers=="alt")
  
  #SIFT_score_i = which(tbx.headers=="SIFT_score")
  #Polyphen2_HDIV_score_i = which(tbx.headers=="Polyphen2_HDIV_score")
  #CADD_raw_i = which(tbx.headers=="CADD_raw")
  
  gr <- GRanges(sub("^chr","",seqnames(vr)), ranges(vr))
  #gr <- unique(gr) # I don't want to remove duplicates because I desire to add dbNSFP annotation results to the original vr.
  message(sprintf("Querying dbNSFP database, number of GRanges = %d.", length(gr)))
  
  mc.cores <- getOption("mc.cores", 4L) ## 4 cores are used by default to query dbNSFP.
  
  # call scanTabix/seqminer::tabix.read to scan dbnsfp database for every mutation
  d <- mclapply(1:length(vr), function(k) {
    e <- vr[k]
    r <- rep(NA, length(tbx.headers))
    #a <- scanTabix(tbx, param = GRanges(sub("^chr","",seqnames(e)), ranges(e)))[[1]]
    # Or
    #a <- unlist(scanTabix(tbx, param = GRanges(sub("^chr","",seqnames(e)), ranges(e)))) ## cannot run in mclapply when mc.cores > 1 because workers cannot access tbx created in the main thread.
    a <- seqminer::tabix.read(tbx$path, sprintf("%s:%s-%s", sub("^chr","",seqnames(e)), start(e), end(e)))
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
    if (verbose && k >= 200 && k %% 200 == 0)
      message(sprintf(" Processed %d mutations.", k))
    return(r)
    }, mc.cores=ifelse(length(vr)>50, mc.cores, 1))
  message("Querying dbNSFP ended.")
  
  d <- do.call("rbind",d)
  d <- as.data.frame(d, stringsAsFactors=FALSE)
  colnames(d) <- tbx.headers
  
  if (!is.null(tbx))
    close(tbx)

  options(warn=-1) # ignore warnings in as.numeric when NA is present.
  #vr$SIFT_score <- as.numeric(as.character(d[, SIFT_score_i]))
  #vr$Polyphen2_HDIV_score <- as.numeric(as.character(d[, Polyphen2_HDIV_score_i]))
  #vr$CADD_raw <- as.numeric(as.character(d[, CADD_raw_i]))
  for (i in 9:ncol(d)) { # Process all scores in a for-loop
    mcols(vr)[,tbx.headers[i]] <- as.numeric(as.character(d[,i]))
  }
  return(vr)
}

get.bkgr.dbnsfp <- function(d, Entrez_Gene_Id, dbnsfp.fl, tbx=NULL, verbose=FALSE) {
  if (is.null(tbx)) {
    tbx <- open(TabixFile(dbnsfp.fl))
  }
  h <- headerTabix(tbx)
  tbx.seqnames <- h$seqnames
  tbx.headers <- strsplit(h$header, "\t")[[1]]
  ref_i = which(tbx.headers=="ref")
  alt_i = which(tbx.headers=="alt")
  
  #SIFT_score_i = which(tbx.headers=="SIFT_score")
  #Polyphen2_HDIV_score_i = which(tbx.headers=="Polyphen2_HDIV_score")
  #CADD_raw_i = which(tbx.headers=="CADD_raw")
  
  gr <- getLongestCCDS(d, Entrez_Gene_Id = Entrez_Gene_Id)
  gr <- GRanges(sub("^chr","",seqnames(gr)), ranges(gr))
  message(sprintf("Querying dbNSFP database for background information for gene %s, number of GRanges = %d.", Entrez_Gene_Id, length(gr)))

  # call scanTabix to scan dbnsfp database for every mutation
  d <- lapply(gr, function(e) {
    r <- rep(NA, length(tbx.headers))
    a <- scanTabix(tbx, param = GRanges(sub("^chr","",seqnames(e)), ranges(e)))[[1]]
    # 2nd way: 
    #a <- unlist(scanTabix(tbx, param = GRanges(sub("^chr","",seqnames(e)), ranges(e)))) ## cannot run in mclapply when mc.cores > 1
    # 3rd way:
    #region <- sprintf("%s:%s-%s", sub("^chr","",seqnames(e)), start(e), end(e))
    #a <- seqminer::tabix.read(tbx$path, region)
    if (length(a) < 1) {
      return(r)  
    }
    b <- strsplit(a,"\t")
    r <- do.call("rbind", b)
    return(r)
  })
  message("Querying dbNSFP for background information ended.")
  
  ##------- remove entries with all NA values. Don't do this in annotate.dbnsfp since the input for annotate.dbnsfp is a VRanges. There is no
  #+error when attach a new feature (even though it's NA) to an existing VRanges. ---------------##
  ## This is definitely required. When there is NA in d, error will occur in following code: vr <- VRanges(...)
  idxs <- which(sapply(d, function(a) any(is.na(a[1:4])) == FALSE)) # a[1:4] are chrom, position, ref and alt, which cannot be NA in VRanges. So I removed it in advance.
  d <- lapply(idxs, function(i) d[[i]])
  ##--------------------------
  d <- do.call("rbind",d)
  d <- as.data.frame(d, stringsAsFactors=FALSE)
  
  ##------- remove entries with all NA values ---------------##
  #d <- na.omit(d)  # Never do this here since some scores are NAs why the others are not. na.omit will remove rows with one or more NA(s). I done this above.
  ##--------------------------
  
  colnames(d) <- tbx.headers

  if (!is.null(tbx))
    close(tbx)

  pos <- as.numeric(d[,2])
  chr <- as.character(seqnames(gr)[1])
  if (!grepl("^chr",chr))
    chr <- sprintf("chr%s", chr)
  
  vr <- VRanges(chr, IRanges(pos, pos), ref = d$ref, alt = d$alt)
  options(warn=-1) # ignore warnings in as.numeric when NA is present.
  vr$aaref <- d$aaref
  vr$aaalt <- d$aaalt
  vr$Hugo_Symbol <- d$genename
  vr$Uniprot_aapos <- d$Uniprot_aapos
  
  #vr$SIFT_score <- as.numeric(d[, SIFT_score_i])
  #vr$Polyphen2_HDIV_score <- as.numeric(d[, Polyphen2_HDIV_score_i])
  #vr$CADD_raw <- as.numeric(d[, CADD_raw_i])
  for (i in 9:ncol(d)) {
    mcols(vr)[, tbx.headers[i]] <- as.numeric(d[,i])
  }
  
  vr <- vr[ref(vr) %in% c("A","C","G","T") & alt(vr) %in% c("A","C","G","T")] # keep only single nucleotide substitutions
  vr <- collapseMutationContext(vr=vr, removeSilent = FALSE)
  return(vr)
}

## Assign 96 mutation types to mutational signatures identified from TCGA dataset (not published yet)
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
  hyper.msig19 <- c("CA:T.T")
  prominent.msigs <- c(msig1,msig2,CtoTinCpG,APOBEC,Aristolochic_Acid_Signature,hyper.msig19)
  others <- setdiff(unique(context), prominent.msigs)

  collapseContext <- rep("others", length(vr))
  collapseContext[context %in% msig1] <- "Signature1.C[C>A]A"
  collapseContext[context %in% msig2] <- "Signature2.T[C>T]C"
  collapseContext[context %in% CtoTinCpG] <- "*CpG"
  collapseContext[context %in% APOBEC] <- "APOBEC"
  collapseContext[context %in% Aristolochic_Acid_Signature] <- "Aristolochic_Acid.C[T>A]G"
  collapseContext[context %in% hyper.msig19] <- "hyper.Signature19.T[C>A]T"
  
  a <- substring(context, 1, 2)
  flag <- !(context %in% prominent.msigs) # If contexts have not been assigned to above signatures.
  
  collapseContext[a == "CA" & flag] <- "others.C>A"
  collapseContext[a == "CG" & flag] <- "others.C>G"
  collapseContext[a == "CT" & flag] <- "others.C>T"
  collapseContext[a == "TA" & flag] <- "others.T>A"
  collapseContext[a == "TC" & flag] <- "others.T>C"
  collapseContext[a == "TG" & flag] <- "others.T>G"
  
  vr$collapseContext <- collapseContext
  return(vr)
}

## vr: object returned by collapseMutationContext
## HugoSymbol: Hugo_Symbol
## EntrezGeneId: Entrez_Gene_Id
## At least must either HugoSymbol or EntrezGeneId be present.

## An example:
## maf <- read.maf(fn)
## vr <- collapseMutationContext(maf)
## obs.vr <- getObsGeneSigature(vr,"VHL")

getObsGeneSigature <- function(vr, HugoSymbol=NULL, EntrezGeneId=NULL) {
  if (!"collapseContext" %in% colnames(mcols(vr)))
    vr <- collapseMutationContext(vr=vr)
  if (!is.null(EntrezGeneId)) {
    gene.vr <- subset(vr, Entrez_Gene_Id == EntrezGeneId)
  }
  else if (!is.null(HugoSymbol)) {
    gene.vr <- subset(vr, Hugo_Symbol == HugoSymbol)
  } else {
    stop("Missing Hugo_Symbol or Entrez_Gene_Id.")
  }
  #d <- data.frame(start=start(gene.vr), collapseContext=gene.vr$collapseContext)
  #return(d)
  return(gene.vr)
}

## Function called by getLongestCCDS(...) and getLongestIsoformPerPosNucleotieGRanges(...).  
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
## Get the CCDS genomic regions of longest isoform corresponding to Entrez_Gene_Id (or Hugo_Symbol).
# d: A data.frame refers to NCBI CCDS file. E.g. /Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2
# Entrez_Gene_Id or Hugo_Symbol: gene identifier, at least provide one.
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
  #d <- subset(d, ccds_status == "Public") # This has been done in read.ncbi.ccds.file(...)
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

## Identify the CCDS genomic regions of longest isoform and return per position per nucleotide GRanges.
# d: A data.frame refers to NCBI CCDS file. E.g. /Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2
# Entrez_Gene_Id or Hugo_Symbol: gene identifier, at least provide one.
getLongestIsoformPerPosNucleotieGRanges <- function(d, Entrez_Gene_Id=NULL, Hugo_Symbol=NULL) 
{
  #d <- subset(d, ccds_status == "Public") # This has been done in read.ncbi.ccds.file(...)
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

## Get background gene signature information.
# d: a data.frame refers to NCBI CCDS file with colname genes, 
#+ e.g.CCDS_current_file="/Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2"
#+ d <- read.table(CCDS_current_file, comment.char = "", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# EntrezGeneId: Entrez_Gene_Id.
# removeSilent: If TRUE call predictCoding(...) to annotate mutations and remove silent ones.
getBkgrGeneSignature <- function(d, HugoSymbol="VHL", EntrezGeneId=7428, removeSilent=FALSE) {
  mutTbl <- list(A=c("C","G","T"), C=c("A","G","T"), G=c("A","C","T"), T=c("A","C","G"))
  gr <- getLongestIsoformPerPosNucleotieGRanges(d, Hugo_Symbol = HugoSymbol, Entrez_Gene_Id = EntrezGeneId)
  if (is.null(gr)) {
    #return(data.frame(start=NA, collapseContext=NA))
    return(VRanges())
  }
  if (exists("Hsapiens")) {
    library("BSgenome.Hsapiens.UCSC.hg19")
  }
  seq <- getSeq(Hsapiens, gr)
  seq <- do.call("c", seq)
  s <- as.character(seq)
  bases <- sapply(1:nchar(s), function(i) substring(s,i,i))
  ref <- rep(bases, each=3)
  alt <- do.call("c",sapply(bases, function(b) mutTbl[b]))
  positions <- rep(start(gr), each=3)
  vr <- VRanges(seqnames(gr)[1], IRanges(positions, positions), ref = ref, alt = alt, sampleNames = "sampleId")
  cmc <- collapseMutationContext(vr=mutationContext(vr, Hsapiens), removeSilent = removeSilent)
  #bkgrSignature <- data.frame(start=start(cmc), collapseContext=cmc$collapseContext)
  #return(bkgrSignature)
  return(cmc)
}

## I failed to implement getBkgrGeneSignature(...) by using TxDb.Hsapiens.UCSC.hg19.knownGene.
## This legacy code is reserved for future work.
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

## Get hotspot statistic.
# y: the positions of mutations in a gene
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

## The core engine to perform MutSigCL.
# bkgr.vr: background mutation list in VRanges format returned by getBkgrGeneSignature(...).
# obs.vr: the observed mutation list in VRanges format returned by getObsGeneSigature(...).
# global: If TRUE perform MutSigCL without taking into account mutation contexts.
# nsim: number of simulations.
mutsigcl_core <- function(bkgr.vr, obs.vr, global=FALSE, nsim=1000, mc.cores=1) {
  
  mutation_sampling <- function(obs.vr, bkgr.vr, global) {
    if (global) { ## Sampling without taking into account mutation context
      x <- start(bkgr.vr)
      k <- length(obs.vr)
      #y <- sample(x, k, replace = TRUE)
      y <- faster.sample(x, k, replace = TRUE)
    } else { ## Sampling according to mutation context
      contexts <- table(obs.vr$collapseContext)
      y <- unlist(lapply(names(contexts), function(context) {
        x <- start(bkgr.vr)[bkgr.vr$collapseContext == context]
        k <- contexts[context]
        #sample(x, k, replace = TRUE) 
        faster.sample(x, k, replace = TRUE) 
      }))
    }
    r <- getHotspotStatistic(y)
    return(r)
  }
  y0 <- start(obs.vr)
  r <- getHotspotStatistic(y0)
  statistic0 <- r[1]
  hotspot.num0 <- r[2]
 
  r <- mclapply(1:nsim, function(i) {
    mutation_sampling(obs.vr, bkgr.vr, global)
  }, mc.cores = mc.cores)
  dat <- do.call("rbind", r)
  k <- sum(dat[, 1] >= statistic0)
  p.value <- (k + 1) / (nsim + 1)
  
  return(c(statistic0, hotspot.num0, p.value))
}

## An example to run mutsigcl_core2:
#-----------------------------------------------
#| ccds.file="/Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2"
#| d <- read.ncbi.ccds.file(ccds.file)
#| maf <- read.maf(maf.file)
#| maf.vr <- collapseMutationContext(maf)
#| Entrez_Gene_Id = "7428"
#| mutsigcl_core2(d, maf.vr, Entrez_Gene_Id)
#------------------------------------------------
mutsigcl_core2 <- function(d, maf.vr, Entrez_Gene_Id, nsim=1000, removeSilent=TRUE, mc.cores=1) {
  bkgr.vr <- getBkgrGeneSignature(d = d, EntrezGeneId = Entrez_Gene_Id, removeSilent = removeSilent)
  obs.vr <- getObsGeneSigature(vr = maf.vr, EntrezGeneId = Entrez_Gene_Id)
  cl <- mutsigcl_core(bkgr.vr, obs.vr, global = FALSE, nsim = nsim, mc.cores = mc.cores)
  global.p.value <- mutsigcl_core(bkgr.vr, obs.vr, global = FALSE, nsim = nsim, mc.cores = mc.cores)[3]
  message(sprintf("Finished running MutSigCL on gene %s, statistic = %g, hotspot.num = %d, p.value = %g, global.p.value = %g.", Entrez_Gene_Id, cl[1],cl[2],cl[3], global.p.value))
  out <- list(statistic=cl[1], hotspot.num=cl[2], p.value=cl[3], global.p.value=global.p.value)
  return(out)
}

## The core engine to perform MutSigFN. By default, MutSigCL is performed by taking
#+into account mutation context and WITHOUT taking into account mutation contexts.
# bkgr.vr: background mutation list in VRanges format returned by get.bkgr.dbnsfp(...).
# obs.vr: the observed mutation list in VRanges format returned by getObsGeneSigature(...).
# dbnsfp.fl: the filename of preprocessed dbNSFP database, see annotate.dbnsfp(...) description for detail.
# nsim: number of simulations.
mutsigfn_core <- function(bkgr.vr, obs.vr, dbnsfp.fl, nsim=1000, mc.cores=1) {
  FN.score.sampling <- function(bkgr.vr, obs.vr, global, score_name) {
    obs.scores <- na.omit(mcols(obs.vr)[,score_name])
    if (global) { # do not take into account mutation context
      bkgr.scores <- na.omit(mcols(bkgr.vr)[,score_name])
      #y <- sample(bkgr.scores, length(obs.scores), replace = TRUE)
      y <- faster.sample(bkgr.scores, length(obs.scores), replace = TRUE)
    } else { # take into account mutation context
      contexts <- table(obs.vr$collapseContext[!is.na(mcols(obs.vr)[,score_name])])
      y <- unlist(lapply(names(contexts), function(context) {
        x <- na.omit(mcols(bkgr.vr)[,score_name][bkgr.vr$collapseContext == context])
        k <- contexts[context]
        #sample(x, k, replace = TRUE)  
        faster.sample(x, k, replace = TRUE)
      }))
    }
    return(sum(y) >= sum(obs.scores))
  }
  FN.score.fast.sampling <- function(bkgr.l, obs.l, global, score_idx) {
    obs.scores <- obs.l[[score_idx]]$FN.score
    if (global) {
      bkgr.scores <- bkgr.l[[score_idx]]$FN.score
      #y <- sample(bkgr.scores, length(obs.scores), replace = TRUE)
      y <- faster.sample(bkgr.scores, length(obs.scores), replace = TRUE)
    } else {
      contexts <- table(obs.l[[score_idx]]$collapseContext)
      y <- unlist(lapply(names(contexts), function(context) {
        x <- bkgr.l[[score_idx]]$FN.score[bkgr.l[[score_idx]]$collapseContext == context]
        k <- contexts[context]
        #sample(x, k, replace = TRUE)  
        faster.sample(x, k, replace = TRUE) 
      }))
    }
    return(sum(y) >= sum(obs.scores))
  }
  
  ##--------------- functions to calculate p-value ----------------------------##
  calc.p.value <- function(bkgr.vr, obs.vr, global, score_name, nsim) {
    k <- sapply(1:nsim, function(i) FN.score.sampling(bkgr.vr, obs.vr, global, score_name))
    return((sum(k) + 1) / (nsim + 1))
  }
  fast.calc.p.value <- function(bkgr.l, obs.l, global, score_idx, nsim) {
    if (nrow(obs.l[[score_idx]]) == 0 || all(is.na(obs.l[[score_idx]]))) {
      # dbNSFP has no annotation for some mutations even they are annotated to be exonic by oncotator.
      # This `if` clause is able to eliminate ERROR message: Error in sum(k) : invalid 'type' (character) of argument, Calls: MutSigCLFN ... mutsigfn_core2 -> mutsigfn_core -> fast.calc.p.value
      #return(1) # Would be more appropriate to return NA since this has no contribution to combine.p.values(..., na.rm=TRUE)
      return(NA)
    }
    #k <- mclapply(1:nsim, function(i) FN.score.fast.sampling(bkgr.l, obs.l, global, score_idx), mc.cores = mc.cores)
    k <- lapply(1:nsim, function(i) FN.score.fast.sampling(bkgr.l, obs.l, global, score_idx))
    k <- unlist(k)
    return((sum(k) + 1) / (nsim + 1))
  }
  ## ---------End---------------------------------------------------------------------##
  
  #bkgr.vr <- annotate.dbnsfp(bkgr.vr, dbnsfp.fl, verbose = TRUE) # bkgr.vr returned from get.bkgr.dbnsfp(...) already contains dbNSFP annotation.
  obs.vr <- annotate.dbnsfp(obs.vr, dbnsfp.fl, verbose = TRUE)
  
  ##------------ transform SIFT scores of bkgr and obs to make it the larger the more damaging --------------##
  #bkgr.vr$SIFT_score <- 1.0 - bkgr.vr$SIFT_score # I have done this in preprocessing the dbNSFP with a Python script
  #obs.vr$SIFT_score <- 1.0 - obs.vr$SIFT_score # I have done this in preprocessing the dbNSFP with a Python script
  
  ## Make sure that the column names listed here are all presented in dbNSFP.
  #score_names <- c("SIFT_score","Polyphen2_HDIV_score","CADD_raw")
  #score_names <- c('SIFT_score','Polyphen2_HDIV_score','LRT_score','MutationAssessor_score','FATHMM_score','RadialSVM_score','LR_score','CADD_raw','GERP++_RS','phyloP46way_primate','phyloP46way_placental','phyloP100way_vertebrate','phastCons46way_primate','phastCons46way_placental','phastCons100way_vertebrate','SiPhy_29way_logOdds','mcap_score')
  #score_names <- c('SIFT_converted_rankscore','Polyphen2_HDIV_rankscore','Polyphen2_HVAR_rankscore','LRT_converted_rankscore','MutationTaster_converted_rankscore','MutationAssessor_rankscore','FATHMM_rankscore','RadialSVM_rankscore','LR_rankscore','CADD_raw_rankscore','GERP++_RS_rankscore','phyloP46way_primate_rankscore','phyloP46way_placental_rankscore','phyloP100way_vertebrate_rankscore','phastCons46way_primate_rankscore','phastCons46way_placental_rankscore','phastCons100way_vertebrate_rankscore','SiPhy_29way_logOdds_rankscore','mcap_score')
  #score_names <- intersect(score_names, colnames(mcols(bkgr.vr))) # in case some score names are not available
  
  ## Make sure that the functional impact and conservative scores are ended with "score". The other column names should not ended with "score".
  score_names <- colnames(mcols(bkgr.vr))[grep('score$',colnames(mcols(bkgr.vr)))]
  
  obs.l <- lapply(score_names, function(score_name) {
    na.omit(data.frame(FN.score=mcols(obs.vr)[,score_name], collapseContext=obs.vr$collapseContext))
  })
  bkgr.l <- lapply(score_names, function(score_name) {
    na.omit(data.frame(FN.score=mcols(bkgr.vr)[,score_name], collapseContext=bkgr.vr$collapseContext))
  })
  
  message("Performing in silico permutation to compute p-values.")
  
  ##----- calculate p-value by assumming Gaussian distribution of FN scores ----------##
  #SIFT.norm.pvalue <- pnorm(median(obs.l[[1]]$FN.score), mean(bkgr.l[[1]]$FN.score), mad(bkgr.l[[1]]$FN.score), lower.tail = FALSE)
  #Polyphen2.norm.pvalue <- pnorm(median(obs.l[[2]]$FN.score), mean(bkgr.l[[2]]$FN.score), mad(bkgr.l[[2]]$FN.score), lower.tail = FALSE)
  #CADD.norm.pvalue <- pnorm(median(obs.l[[3]]$FN.score), mean(bkgr.l[[3]]$FN.score), mad(bkgr.l[[3]]$FN.score), lower.tail = FALSE)
  
  ##-------------------- calcuate p-value with permutation ----------------------##
  global.pvalues <- mclapply(1:length(score_names), function(i) {
    fast.calc.p.value(bkgr.l, obs.l, global = TRUE, i, nsim)
  }, mc.cores = mc.cores)
  global.pvalues <- unlist(global.pvalues)
  combined.global.pvalue <- combine.p.values(global.pvalues)
  
  pvalues <- mclapply(1:length(score_names), function(i) {   ## pvalues by considering contexts.
    fast.calc.p.value(bkgr.l, obs.l, global = FALSE, i, nsim)
  }, mc.cores = mc.cores)
  pvalues <- unlist(pvalues)
  combined.pvalue <- combine.p.values(pvalues)
  ##--------------------End--------------------------------------------------------##
  
  ## --------------- get the median SIFT, Polyphen2_HDIV and mcap scores -----------##
  SIFT_score_name <- score_names[grep("SIFT", score_names)][1]
  Polyphen2_score_name <- score_names[grep("Polyphen2_HDIV", score_names)][1] # Make sure that length(Polyphen2_score_name) == 1, or ERROR will occur in: median(mcols(obs.vr)[,Polyphen2_score_name], na.rm = TRUE)
  mcap_score_name <- score_names[grep("mcap_score", score_names)][1]
  
  SIFT_score <- median(mcols(obs.vr)[,SIFT_score_name], na.rm = TRUE)
  Polyphen2_score <- median(mcols(obs.vr)[,Polyphen2_score_name], na.rm = TRUE)
  mcap_score <- median(mcols(obs.vr)[,mcap_score_name], na.rm = TRUE)
  ## ----------------End---------------------------------------------------------------##
  
  out <- c(SIFT_score, Polyphen2_score, mcap_score, global.pvalues, pvalues, combined.global.pvalue, combined.pvalue)
  headers <- c("SIFT_rankscore", "Polyphen2_rankscore", "mcap_score", paste(score_names, "global.pvalue", sep = "."), paste(score_names, "pvalue", sep = "."), "combined.global.pvalue", "combined.pvalue")
  names(out) <- headers
  gene_id <- obs.vr$Hugo_Symbol[1]
  message(sprintf("Finished running MutSigFN on %s, SIFT_rankscore = %g, Polyphen2_HDIV_rankscore = %g, mcap_score = %g, combined.global.pvalue = %g, combined.pvalue = %g.", gene_id, SIFT_score, Polyphen2_score, mcap_score, combined.global.pvalue, combined.pvalue))
  return(out)
}

# d: A data.frame refers to NCBI CCDS file. E.g. /Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2
# maf.vr: the VRanges object of observed mutations annotated by oncotator.
# e.g.:
#     maf = read.maf(...)
#     maf.vr = collapseMutationContext(maf)
# Entrez_Gene_Id: gene id
# dbnsfp.fl: the filename of dbNSFP.
# nsim: number of simutations.

## An example to run mutsigfn_core2:
#-----------------------------------------------
#| ccds.file="/Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2"
#| d <- read.ncbi.ccds.file(ccds.file)
#| maf <- read.maf(maf.file)
#| maf.vr <- collapseMutationContext(maf)
#| dbnsfp.fl = "dbsnfp2.4.gz"
#| Entrez_Gene_Id = "7428"
#| mutsigfn_core2(d, maf.vr, Entrez_Gene_Id, dbnsfp.fl)
#------------------------------------------------
mutsigfn_core2 <- function(d, maf.vr, Entrez_Gene_Id, dbnsfp.fl, nsim=1000, mc.cores=1) {
  tic0 <- proc.time()[3] # elapsed time
  ## It's not necessary to remove silent mutations here since dbNSFP does not contain silent
  #+mutations. Removing silent mutations may be able to reduce the burden of querying dbNSFP.
  #bkgr.vr <- getBkgrGeneSignature(d = d, EntrezGeneId = Entrez_Gene_Id, removeSilent = FALSE)
  bkgr.vr <- get.bkgr.dbnsfp(d, Entrez_Gene_Id = Entrez_Gene_Id, dbnsfp.fl = dbnsfp.fl)
  obs.vr <- getObsGeneSigature(vr = maf.vr, EntrezGeneId = Entrez_Gene_Id)
  out <- mutsigfn_core(bkgr.vr, obs.vr, dbnsfp.fl, nsim, mc.cores)
  message(sprintf("Elapsed time is %g secs.", proc.time()[3] - tic0))
  return(out)
}

# maf.file: the filename of mutation file annotated by oncotator
# ccds.file: NCBI CCDS file
# maf: an object returned by read.maf(...)
# d: a data.frame referred to ccds.file
# k: the minimum mutation frequency of genes to be analyzed.
# removeSilent: If TRUE, silent mutations are excluded from analysis.
# nsim: number of simulations.
# mc.cores: number of cpu cores used. Currently, mc.cores > 1 does not work.
mutsigclOldVersion <- function(maf.file, ccds.file, out.file=NULL, maf=NULL, d=NULL, removeSilent=TRUE, nsim=1000, mc.cores=1) {
  ## loading required packages
  library(maftools)
  library(VariantAnnotation)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(SomaticSignatures)
  library(fastcluster)
  library(seqminer)
  
  if (is.null(maf)) {
    #maf <- read.maf(maf.file, removeSilent=removeSilent)
    vr <- makeVRangesFromOncotator(maf.file = maf.file, removeIndel = TRUE, removeSilent = removeSilent)
    message(sprintf("Finished reading %s.", maf.file))
  }
  if (is.null(d)) {
    d <- read.ncbi.ccds.file(ccds.file)
  }
  
  ## Since read.maf(...) has already removed silent mutations, so there is no need to do
  ##+it here. But if you are serious about consistent annotation results, I recommend to
  ##+set removeSilent to TRUE and mutations will be re-annotated by VariantAnnotation
  ##+package which is used in background in silico mutation annotation.
  vr <- collapseMutationContext(maf, vr = vr, removeSilent = FALSE)
  
  x <- by(start(vr),vr$Entrez_Gene_Id, getHotspotStatistic)
  hotspot.num <- do.call("rbind", x)[,2]
  hotspotMutatedGeneIds <- names(hotspot.num)[!is.na(hotspot.num) & hotspot.num > 0]

  geneToMutatedSamples <- by(vr@sampleNames, vr$Entrez_Gene_Id, function(a) length(unique(a)))
  recurrentlyMutatedGeneIds <- names(geneToMutatedSamples)[geneToMutatedSamples>=3]
  hotspotMutatedGeneIds <- intersect(hotspotMutatedGeneIds, recurrentlyMutatedGeneIds)

  ###############--MUC16---TTN--
  long_genes <- c("94025","7273")
  hotspotMutatedGeneIds <- setdiff(hotspotMutatedGeneIds, long_genes)
  hotspotMutatedGeneIds <- setdiff(hotspotMutatedGeneIds, as.character(d$gene_id))
  hotspotMutatedGeneIds <- c("7157","7428")
  
  message(sprintf("Number of genes with hotspot mutations = %d.", length(hotspotMutatedGeneIds)))
  
  dat <- data.frame(Hugo_Symbol=vr$Hugo_Symbol, Entrez_Gene_Id=vr$Entrez_Gene_Id)
  dat <- apply(dat, 2, function(a) sub("\\s+","",as.character(a)))
  
  r <- mclapply(hotspotMutatedGeneIds, function(geneId, vr, d) {
    #library(BSgenome.Hsapiens.UCSC.hg19) ## For running in parallel, loading BSgenome in each worker is required. NOT WORK!!!!!!
    bkgr.vr <- getBkgrGeneSignature(d, EntrezGeneId = geneId, removeSilent = removeSilent)
    if (length(bkgr.vr) == 0)
      return(c(NA,NA,NA,NA))
    obs.vr <- getObsGeneSigature(vr, EntrezGeneId = geneId)
    res <- mutsigcl_core(bkgr.vr, obs.vr, global = FALSE, nsim = nsim)
    global.p.value <- mutsigcl_core(bkgr.vr, obs.vr, global = TRUE, nsim = nsim)[3]
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
  my.df$q.value <- p.adjust(as.numeric(my.df$p.value), method="BY")
  my.df$global.q.value <- p.adjust(as.numeric(my.df$global.p.value), method="BY")
  my.df$MutatedSamples <- sapply(hotspotMutatedGeneIds, function(e) geneToMutatedSamples[[e]])
  my.df <- my.df[order(my.df$q.value, decreasing = FALSE), ]
  if (!is.null(out.file)) {
    write.table(my.df, file = out.file, quote=FALSE, sep = "\t", row.names = FALSE)
  }
  invisible(my.df)
}

# mc.cores: number of cpu cores used in querying dbNSFP database.
MutSigFN <- function(maf.file, ccds.file, dbnsfp.file, out.file=NULL, nsim=10000, minSampleNum=3, mc.cores=NULL, hyperCutoff=NULL) {
  options(mc.cores = mc.cores)
  
  maf.vr <- makeVRangesFromOncotator(maf.file, removeIndel = TRUE, removeSilent = TRUE, hyperCutoff = hyperCutoff)
  d <- read.ncbi.ccds.file(ccds.file = ccds.file)
  
  a = unique(data.frame(maf.vr$Entrez_Gene_Id, maf.vr$Hugo_Symbol))
  Hugo_Symbols <- as.character(a[,2])
  names(Hugo_Symbols) <- as.character(a[,1])
  
  genesToMutatedSampleNum <- by(sampleNames(maf.vr), maf.vr$Entrez_Gene_Id, function(a) length(unique(a))) 
  recurrentlyMutatedGenes <- names(genesToMutatedSampleNum)[genesToMutatedSampleNum >= minSampleNum]
  recurrentlyMutatedGenes <- intersect(as.character(d$gene_id), recurrentlyMutatedGenes)
  
  ###############--MUC16---TTN--
  long_genes <- c("94025","7273")
  recurrentlyMutatedGenes <- setdiff(recurrentlyMutatedGenes, long_genes)
  message(sprintf("Number of recurrently mutated genes = %d.", length(recurrentlyMutatedGenes)))
  
  r = lapply(recurrentlyMutatedGenes, function(gene_id) {
    out <- mutsigfn_core2(d=d, maf.vr = maf.vr, Entrez_Gene_Id = gene_id, dbnsfp.fl = dbnsfp.file, nsim = nsim, mc.cores = mc.cores)
    return(unlist(out))
  })
  r = as.data.frame(do.call("rbind", r))
  mutatedSampleNums = sapply(recurrentlyMutatedGenes, function(gene_id) {
    genesToMutatedSampleNum[[gene_id]]
  })
  Hugo_Symbols = sapply(recurrentlyMutatedGenes, function(gene_id) {
    Hugo_Symbols[[gene_id]]
  })
  r <- data.frame(Hugo_Symbol=Hugo_Symbols, Entrez_Gene_Id=recurrentlyMutatedGenes, mutatedSampleNum=mutatedSampleNums, r)
  ##------------ correct for multiple hypothesis test -----------##
  r$BH.global.q.value <- p.adjust(r$combined.global.pvalue, method = "BH")
  r$BY.global.q.value <- p.adjust(r$combined.global.pvalue, method = "BY")
  r$BH.q.value <- p.adjust(r$combined.pvalue, method = "BH")
  r$BY.q.value <- p.adjust(r$combined.pvalue, method = "BY")
  ##------------------------------------------------------------##
  r <- r[order(r$BY.q.value, decreasing = FALSE),]
  if (!is.null(out.file))
    write.table(r, file=out.file, quote = FALSE, sep = "\t", row.names = FALSE)
  invisible(r)
}
MutSigCL <- function(maf.file, ccds.file, out.file=NULL, removeSilent=TRUE, nsim=10000, minSampleNum=3, mc.cores=1, hyperCutoff=NULL) {

  vr <- makeVRangesFromOncotator(maf.file = maf.file, removeIndel = TRUE, removeSilent = removeSilent, hyperCutoff = hyperCutoff)
  message(sprintf("Finished reading %s.", maf.file))
  
  d <- read.ncbi.ccds.file(ccds.file)
  
  ## Since makeVRangesFromOncotator(...) has already removed silent mutations if removeSilentis 
  ##+set to TRUE, so there is no need to do
  ##+it here. But if you are serious about consistent annotation results, I recommend to
  ##+set removeSilent to TRUE and mutations will be re-annotated by VariantAnnotation
  ##+package which is used in background in silico mutation annotation.
  vr <- collapseMutationContext(vr = vr, removeSilent = FALSE)
  
  x <- by(start(vr),vr$Entrez_Gene_Id, getHotspotStatistic)
  hotspot.num <- do.call("rbind", x)[,2]
  hotspotMutatedGeneIds <- names(hotspot.num)[!is.na(hotspot.num) & hotspot.num > 0]
  
  geneToMutatedSampleNum <- by(vr@sampleNames, vr$Entrez_Gene_Id, function(a) length(unique(a)))
  recurrentlyMutatedGeneIds <- names(geneToMutatedSampleNum)[geneToMutatedSampleNum >= minSampleNum]
  hotspotMutatedGeneIds <- intersect(hotspotMutatedGeneIds, recurrentlyMutatedGeneIds)
  
  ###############--MUC16---TTN--
  long_genes <- c("94025","7273")
  hotspotMutatedGeneIds <- setdiff(hotspotMutatedGeneIds, long_genes)
  hotspotMutatedGeneIds <- intersect(hotspotMutatedGeneIds, as.character(d$gene_id))
  ##hotspotMutatedGeneIds <- c("7157","7428") ## TP53, VHL
  
  message(sprintf("Number of genes with hotspot mutations = %d.", length(hotspotMutatedGeneIds)))
  
  r <- lapply(hotspotMutatedGeneIds, function(geneId, vr, d) {
    res <- mutsigcl_core2(d = d, maf.vr = vr, Entrez_Gene_Id = geneId, nsim = nsim, removeSilent = removeSilent, mc.cores = mc.cores)
    return(do.call("cbind", res))
  }, vr, d)
  
  dat <- data.frame(Hugo_Symbol=vr$Hugo_Symbol, Entrez_Gene_Id=vr$Entrez_Gene_Id, stringsAsFactors = FALSE)
  Hugo_Symbols <- sapply( hotspotMutatedGeneIds, function(geneId) return(dat[,1][dat[,2] == geneId][1]) )
  
  r <- do.call("rbind", r)
  r <- cbind(Hugo_Symbol=Hugo_Symbols, Entrez_Gene_Id=names(Hugo_Symbols), r)
  r <- as.data.frame(r, stringsAsFactors = FALSE)
  r$q.value <- p.adjust(as.numeric(r$p.value), method="BY")
  r$global.q.value <- p.adjust(as.numeric(r$global.p.value), method="BY")
  r$MutatedSampleNum <- sapply(hotspotMutatedGeneIds, function(e) geneToMutatedSampleNum[[e]])
  r <- r[order(r$q.value, decreasing = FALSE), ]
  if (!is.null(out.file)) {
    write.table(r, file = out.file, quote=FALSE, sep = "\t", row.names = FALSE)
  }
  invisible(r)
}

## The main engine to run MutSigCL and MutSigFN.
# maf.file: oncotator annotated file
# ccds.file: NCBI CCDS file
# dbnsfp.file: pre-processed dbNSFP database when running MutSigFN.
MutSigCLFN <- function(maf.file, ccds.file, dbnsfp.file=NULL, out.file=NULL, removeSilent=TRUE, method=c("CL","FN"), nsim=10000, minSampleNum=3, mc.cores=NULL, hyperCutoff=NULL) {
  method_selected <- match.arg(method)
  if (method_selected == "FN" && is.null(dbnsfp.file)) {
    stop("You must provide the filename of dbNSFP for dbnsfp.file to run MutSigFN.")
  }
  tic0 <- proc.time()[1]
  if (method_selected == "CL") {
    r <- MutSigCL(maf.file, ccds.file, out.file = out.file, removeSilent = removeSilent, nsim = nsim, minSampleNum = minSampleNum, mc.cores = mc.cores, hyperCutoff = hyperCutoff)
  } else if (method_selected == "FN") {
    r <- MutSigFN(maf.file, ccds.file, dbnsfp.file = dbnsfp.file, out.file = out.file, nsim = nsim, minSampleNum = minSampleNum, mc.cores = mc.cores, hyperCutoff = hyperCutoff)
  }
  tic <- proc.time()[1]
  message(sprintf("You have successfully finished running MutSig%s. It takes %.3f secs.", method_selected, tic - tic0))
  invisible(r)
}

#ccds.file="/Users/lixiangchun/BaiduYunPan/BaiduYunPan/Database/ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt.bz2"

