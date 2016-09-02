
## source code from https://seqqc.wordpress.com/2015/03/10/calculate-nucelotide-frequency-with-rsamtools-pileup/
pileup_freqs <- function(plp) {
    nucleotides <- levels(plp$nucleotide)
    res <- split(plp, plp$seqnames)
    res <- lapply(res, function (x) {split(x, x$pos)})
    res <- lapply(res, function (positionsplit) {
        nuctab <- lapply(positionsplit, function(each) {
                        chr = as.character(unique(each$seqnames))
                        pos = as.character(unique(each$pos))
                        tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
                        c(chr,pos, tablecounts)
                    })
        nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
        rownames(nuctab) <- NULL
        nuctab
    })
    res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
    rownames(res) <- NULL
    colnames(res) <- c("seqnames","start",levels(plp$nucleotide))
    res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
    res
}


