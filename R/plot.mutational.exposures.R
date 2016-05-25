

plot.mutational.exposures <-
  function(d = NULL,
           mutational.exposure.file = NULL,
           clustid = NULL,
           cols =  c('#6C4F4C','#8EA37A','#843818','#AEDD99','#F89867','#FF0000','#A79D90'),
           sorting = FALSE,
           use.colname.as.signature.name = FALSE,
           legend.ncol = 4,
           legend.ypos = 0,
           sorting.by.sample.name = FALSE,
           add.sample.name = FALSE)
  {
    orderDescending <- function(...)
      order(..., decreasing = TRUE)
    mutational.exposure.file.set <- FALSE
    if (is.null(d) && file.exists(mutational.exposure.file)) {
      d <- read.table(mutational.exposure.file,header=TRUE)
      d <- d[,c(2:ncol(d))]
      d <- as.data.frame(t(d))
      sample_names <- sub('_VS_.*','',rownames(d))
      sample_names <- sub('-VS-.*','',sample_names)
      rownames(d) <- sample_names
      if (sorting.by.sample.name) {
        sample_names <- sort(sample_names)
        d <- d[sample_names,]
      }
      mutational.exposure.file.set <- TRUE
    }
	if (!is.null(d) && !is.data.frame(d)) {
		d <- as.data.frame(d)
	}
    if (is.null(clustid)) {
      clustid <- rep(1,nrow(d))
    }
    if (!is.factor(clustid)) {
      clustid <- as.factor(clustid)
    }
    group <- as.integer(clustid)
    mutationalSignatureNames <- colnames(d)
    d$group <- group
    col.names <- c('group', mutationalSignatureNames)
    ## Sort by columns sequentially
    if (sorting) {
      d <- d[do.call('orderDescending', d[, col.names]),]
    }
    d$group <- NULL
    a <-
      colSums(d[, mutationalSignatureNames]) / sum(d[, mutationalSignatureNames])
    d <- prop.table(as.matrix(d), 1)
    legend.text <-
      paste('Signature', 1:length(mutationalSignatureNames))
    colnames(d) <- NULL
    rownames(d) <- NULL
    if (use.colname.as.signature.name) {
      legend.text <-
        sprintf("%s (%.3g%s)", mutationalSignatureNames, a * 100, '%')
    } else {
      legend.text <-
        sprintf("Signature %d (%.3g%s)",
                1:length(mutationalSignatureNames),
                a * 100,
                '%')
    }
    #pdf('MutationalSignature.pdf', height=3.5, width=9)
    #cols <- c('#6C4F4C', '#8EA37A', '#843818', '#AEDD99', '#F89867', 'black', '#A79D90')
    p <-
      barplot(
        t(d),
        ylab = 'Percentage of mutations',
        space = 0,
        las = 1,
        border = NA,
        col = cols,
        legend.text = legend.text,
        args.legend = list(
          x = max(dim(d)),
          y = legend.ypos,
          ncol = legend.ncol,
          xpd = TRUE,
          border = NA,
          box.col = NA,
          x.intersp = 0.1,
          cex = 0.8
        )
      )
    if (add.sample.name && mutational.exposure.file.set) {
      mtext(sample_names, 1, at=1:nrow(d),las=2,adj=1,padj=-0.8)
      if (legend.ypos < 1) {
        print('INFO - You may need to set legend.ypos = 1.15 when add.sample.name is set to TRUE.')
      }
    }
    #print(idx)
    idx <- table(clustid)
    idx <- cumsum(idx)
    if (length(idx) > 1) {
      for (i in idx[1:(length(idx) - 1)]) {
        abline(
          v = p[i],
          lwd = 0.7,
          lty = 'dashed',
          col = 'black'
        )
      }
    }
    #dev.off()
  }


#library(hustlxc)
#data('icgc.hcc')
#d <- icgc.hcc
#mutationalSignatureNames <- paste('mutationalSignature',1:7,sep="")
#clustid <- d$nmf_clustid
#clustid <- rep(0, nrow(d))
#plot.mutational.exposures(d[,mutationalSignatureNames],clustid)
