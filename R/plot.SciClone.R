addHighlightLegend2 <- function (data, positionsToHighlight, scale) 
{
    if ((is.null(positionsToHighlight))) {
        return()
    }
    names(positionsToHighlight) = c("chr", "st", "name")
    addpts = merge(data, positionsToHighlight, by.x = c("chr", 
        "st"), by.y = c("chr", "st"))
    if (length(addpts[, 1]) == 0) {
        return()
    }
    non.trivial.names <- addpts$name[addpts$name != ""]
    if (length(non.trivial.names) == 0) {
        return()
    }
    par(new=TRUE)

    #plot.default(x = -10000, y = 1, type = "p", pch = 19, cex = 0.4, 
     #   col = "#00000000", xlim = c(0, 1000), ylim = c(0, 1000), 
      #  axes = FALSE, ann = FALSE)

    plot.default(x = -10000, y = 1, type = "p", pch = 19, cex = 0.4, 
        col = "#00000000", xlim = c(0, 100), ylim = c(0, 1000), 
        axes = FALSE, ann = FALSE, font=3)

    #ypos = rev(seq(0, 900, (900/13)))[1:13]
	#step.w = 12
    #ypos = rev(seq(0, 900, (800/step.w)))[1:step.w]

	step.w = nrow(na.omit(positionsToHighlight))
    ypos = rev(seq(20, 900, length.out=step.w))
	if (step.w == 1) { ypos = c(100) }
	if (step.w == 2) { ypos = rev(seq(100, 200, 100)) }
	if (step.w == 3) { ypos = rev(seq(100, 300, 100)) }
	if (step.w == 4) { ypos = rev(seq(100, 400, 100)) }
	if (step.w == 5) { ypos = rev(seq(100, 500, 100)) }
	if (step.w == 6) { ypos = rev(seq(100, 600, 100)) }
	if (step.w <= 6) {ypos = ypos - 50}

    non.trivial.indices <- (1:length(addpts$name))[addpts$name != ""]
    ncol = ceiling(length(non.trivial.names)/step.w)
    xpos = 40
    offset = 1
    nxt <- 1
    for (n in 1:ncol) {
        names = non.trivial.names[offset:(offset + 12)]
        names = as.character(names[!(is.na(names))])
        num = length(names)
        for (i in 1:num) {
            text(xpos, ypos[i], paste(nxt, ". ", names[i], sep = ""), 
                cex = 0.5, pos = 4, font=3, col='grey40')
            nxt <- nxt + 1
        }
        xpos = xpos + 250
        offset = offset + step.w
    }
}

plot.SciClone <- function (sco, outputFile, cnToPlot = c(1, 2, 3, 4), showCopyNumberScatterPlots = TRUE, 
    highlightSexChrs = TRUE, positionsToHighlight = NULL, highlightsHaveNames = FALSE, 
    overlayClusters = TRUE, overlayIndividualModels = TRUE, showHistogram = FALSE, 
    showTitle = TRUE, biggerText = FALSE, highlightsOnHistogram = FALSE, 
    highlightCnPoints = FALSE, smp=NULL) 
{
    densityData = sco@densities
    vafs.merged = sco@vafs.merged
    sampleNames = sco@sampleNames
    dimensions = sco@dimensions
	#nclust = length(unique(sc@clust$cluster.assignments))
	nclust = dim(sc@clust$cluster.lower)[2]
    if (max(cnToPlot) > 4 | min(cnToPlot) < 1) {
        print("sciClone supports plotting of copy numbers between 1 and 4 at this time")
    }
    minimumLabelledPeakHeight = 999
    onlyLabelHighestPeak = TRUE
    add.legend <- FALSE
    addpts <- NULL
    if (!is.null(positionsToHighlight)) {
        names(positionsToHighlight) = c("chr", "st", "name")
        addpts = merge(vafs.merged, positionsToHighlight, by.x = c("chr", 
            "st"), by.y = c("chr", "st"))
        if ((dim(addpts)[1] > 0) & (any(addpts$name != ""))) {
            if (showCopyNumberScatterPlots & (length(cnToPlot) < 
                2) & highlightsHaveNames) {
                add.legend <- TRUE
            }
        }
    }
    if (highlightsHaveNames & add.legend) {
        if (!(is.null(positionsToHighlight))) {
            if (length(positionsToHighlight) < 3) {
                print("ERROR: named plot requires names in the third column of the positionsToHighlight data frame")
                return(0)
            }
            cnToPlot = c(2)
        }
        else {
            print("ERROR: highlightsHaveNames requires a 3-column dataframe of positions and names (chr, pos, name)")
            return(0)
        }
    }
    clust = NULL
    if (overlayClusters) {
        if (is.null(sco@clust[1])) {
            print("ERROR: can't overlay clusters when clustering was not done on the input data")
            return(0)
        }
        else {
            clust = sco@clust
        }
    }
    num.rows <- length(cnToPlot) + 1
    if (add.legend) {
    #    num.rows <- num.rows + 1
    }
    textScale = 1
    axisPosScale = 1
    if (biggerText) {
        textScale = 1.4
        axisPosScale = 0.9
    }
    height <- 8.5 * (num.rows/5)
    width <- 3.7
    spacing = 1
    scale = 1
    if (num.rows == 2) {
        spacing = 1.5
        scale = 1.5
    }
    if (num.rows == 3) {
        spacing = 1
        scale = 1
    }
    pdf(file = outputFile, width = width, height = height, bg = "white")
    numClusters = 0
    if (!is.null(clust)) {
        numClusters = max(clust$cluster.assignments)
    }
    for (d in 1:dimensions) {
        name = sampleNames[d]
        vafs = getOneSampleVafs(vafs.merged, name, numClusters)
        par(mfcol = c(num.rows, 1), mar = c(0.5, 3/spacing, 1, 
            1.5/spacing), oma = c(3/spacing, 0.5, 4/spacing, 
            0.5), mgp = c(3, 1, 0))
        densities = densityData[[d]]$densities
        factors = densityData[[d]]$factors
        peakPos = densityData[[d]]$peakPos
        peakHeights = densityData[[d]]$peakHeights
        maxDepth = densityData[[d]]$maxDepth
        maxDensity = densityData[[d]]$maxDensity
        scalingFactor = 1/maxDensity
        plot.default(x = c(1:10), y = c(1:10), ylim = c(0, 1.1), 
            xlim = c(0, 100), axes = FALSE, ann = FALSE, col = "#00000000", 
            xaxs = "i", yaxs = "i")
        rect(0, 0, 100, 1.1, col = "#00000011", border = NA)
        axis(side = 2, at = c(0, 2), labels = c("", ""), las = 1, 
            cex.axis = 0.6 * textScale, hadj = 0.6 * textScale, 
            lwd = 0.5, lwd.ticks = 0.5, tck = -0.01)
        colors = c("#1C366099", "#67B32E99", "#F4981999", "#E5242099")
        density.curve.width <- 4
        for (i in cnToPlot) {
            if (!(is.null(densities[[i]])) & (!showHistogram | 
                (i != 2))) {
                lines(densities[[i]]$x, scalingFactor * factors[[i]], 
                  col = colors[i], lwd = density.curve.width)
                if (length(peakHeights[[i]]) > 0) {
                  ppos = c()
                  if (onlyLabelHighestPeak) {
                    ppos = which(peakHeights[[i]] == max(peakHeights[[i]]))
                  }
                  else {
                    ppos = which((peakHeights[[i]] == max(peakHeights[[i]])) & 
                      (peakHeights[[i]] > minimumLabelledPeakHeight))
                  }
                  if (length(ppos) > 0) {
                    text(x = peakPos[[i]][ppos], y = (scalingFactor * 
                      peakHeights[[i]][ppos]) + 1.7, labels = signif(peakPos[[i]][ppos], 
                      3), cex = 0.7, srt = 0, col = colors[[i]])
                  }
                }
            }
            else if (showHistogram & (i == 2)) {
                v = vafs[which(vafs$cleancn == 2 & vafs$adequateDepth == 
                  1), ]
                frequencies <- data.frame(x = v$vaf, row.names = NULL, 
                  stringsAsFactors = NULL)
                bin.width <- 2.5
                num.breaks <- ceiling(100/bin.width) + 1
                breaks <- unlist(lapply(0:(num.breaks - 1), function(x) 100 * 
                  x/(num.breaks - 1)))
                h <- hist(v$vaf, breaks = breaks, plot = FALSE)
                h$density <- h$density/max(h$density)
                plot(h, add = TRUE, freq = FALSE, col = "white", 
                  border = "black")
            }
        }
        model.style <- 4
        model.style <- 1
        model.width <- density.curve.width/2
        individual.model.style <- 3
        individual.model.width <- density.curve.width/2
        if (!(is.null(clust))) {
            maxFitDensity <- max(clust$fit.y[d, ])
            lines(clust$fit.x, clust$fit.y[d, ]/maxFitDensity, 
                type = "l", col = "grey50", lty = model.style, 
                lwd = model.width)
            if (overlayIndividualModels == TRUE) {
                for (i in 1:numClusters) {
                  lines(clust$fit.x, clust$individual.fits.y[[i]][d, 
                    ]/maxFitDensity, type = "l", col = "grey50", 
                    lty = individual.model.style, lwd = individual.model.width)
                }
            }
            if (highlightsOnHistogram) {
                if (!is.null(positionsToHighlight)) {
                  addpts = merge(vafs, positionsToHighlight, 
                    by.x = c("chr", "st"), by.y = c("chr", "st"))
                  for (i in 1:length(addpts$vaf)) {
                    if (addpts$name[i] != "") {
                      vaf <- addpts$vaf[i]
                      nearest.indx <- which(unlist(lapply(clust$fit.x, 
                        function(x) abs(x - vaf))) == min(abs(clust$fit.x - 
                        vaf)))[1]
                      vaf.y <- clust$fit.y[d, nearest.indx]/maxFitDensity
                      label <- as.character(addpts$name[i])
                      cex <- 1
                      text(x = vaf, y = vaf.y, label = "*", cex = cex)
                      text(x = vaf, y = vaf.y + 0.1, label = label, 
                        cex = cex)
                    }
                  }
                }
            }
        }
        lcol = colors[cnToPlot]
        lty = c(1, 1, 1, 1)
        lwd = c(2, 2, 2, 2)
        pchs = c(NA, NA, NA, NA)
        pt.bgs = lcol
        leg = c("1 Copy", "2 Copies", "3 Copies", "4 Copies")
        leg = leg[cnToPlot]
        if ((length(cnToPlot) == 1) & (cnToPlot[1] == 2)) {
            if (showHistogram == FALSE) {
                lty = c(1)
                lwd = c(2)
                pchs = c(NA)
                pt.bgs = lcol
            }
            else {
                lcol = "black"
                lty = c(0)
                lwd = c(0)
                pchs = c(22)
                pt.bgs = "white"
            }
        }
        if (!(is.null(clust))) {
            leg = c(leg, "Model Fit")
            lcol = c(lcol, "grey50")
            lty = c(lty, model.style)
            lwd = c(lwd, model.width)
            pt.bgs = c(pt.bgs, "grey50")
            pchs = c(pchs, NA)
            if (overlayIndividualModels == TRUE) {
                leg = c(leg, "Component Fits")
                lcol = c(lcol, "grey50")
                lty = c(lty, individual.model.style)
                lwd = c(lwd, 2)
                pt.bgs = c(pt.bgs, "grey50")
                pchs = c(pchs, NA)
            }
        }
        legend(x = "topright", lwd = lwd, lty = lty, legend = leg, 
            col = lcol, bty = "n", cex = ((0.6/scale) * textScale), 
            y.intersp = 1.25, pch = pchs, pt.bg = pt.bgs)
        axis(side = 3, at = c(0, 20, 40, 60, 80, 100), labels = c(0, 
            20, 40, 60, 80, 100), cex.axis = ((0.6/scale) * textScale), 
            lwd = 0.5, lwd.ticks = 0.5, padj = (((scale * 3.5) - 
                3.5) + 1.4) * (1/textScale), tck = -0.05)

        #mtext("Variant Allele Frequency", adj = 0.5, padj = -3.1 * (1/textScale), cex = 0.6 * (textScale * 0.8), side = 3)
        mtext(paste("Variant Allele Frequency (", smp, ", n=", nclust, ")", sep=""), adj = 0.5, padj = -3.1 * (1/textScale), cex = 0.6 * (textScale * 0.8), side = 3)


        mtext("Density (a.u.)", side = 2, cex = 0.6 * (textScale * 
            0.8), padj = -4.2 * (axisPosScale))
        if (showTitle) {
            title = ""
            if (is.null(sampleNames[d])) {
                title = "Clonality Plot"
            }
            else {
                title = paste(sampleNames[d], "Clonality Plot", sep = " ")
            }
            mtext(title, adj = 0.5, padj = -5 * (1/textScale), 
                cex = 0.65 * textScale, side = 3)
        }
        if (showCopyNumberScatterPlots) {
            for (i in cnToPlot) {
                v = vafs[which(vafs$cleancn == i & vafs$adequateDepth == 
                  1), ]
                drawScatterPlot(v, highlightSexChrs, positionsToHighlight, 
                  colors, i, maxDepth, highlightsHaveNames, overlayClusters, 
                  scale, textScale, axisPosScale, labelOnPlot = FALSE, 
                  highlightCnPoints = highlightCnPoints)
                axis(side = 1, at = c(0, 20, 40, 60, 80, 100), 
                  labels = c(0, 20, 40, 60, 80, 100), cex.axis = (0.6/scale) * 
                    textScale, lwd = 0.5, lwd.ticks = 0.5, padj = ((-scale * 
                    5) + 5 - 1.4) * (1/textScale), tck = -0.05)
                if (length(cnToPlot) < 2 & highlightsHaveNames) {
                  addHighlightLegend2(v, positionsToHighlight,scale)
                }
                else {
                  if (highlightsHaveNames) {
                    print("WARNING: highlighted point naming is only supported when plotting only CN2 regions (cnToPlot=c(2))")
                    print("Instead labeling directly on plot")
                  }
                }
            }
        }
    }
    devoff <- dev.off()
}

#library(sciClone)
#files = list.files(pattern='.Rdata')

#for (fn in files) {
#	load(fn)
#	fig=paste(smp,".pdf",sep="")
#	plot1d(sc,fig, highlightSexChrs=FALSE, positionsToHighlight=anno.smp, highlightsHaveNames=TRUE, overlayClusters=TRUE, showTitle=FALSE, cnToPlot=c(2), biggerText=FALSE)
#	cat(smp, fn, fig, "\n")
#}
