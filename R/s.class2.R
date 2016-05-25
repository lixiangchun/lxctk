"dudi.acm" <- function (df, row.w = rep(1, nrow(df)), scannf = TRUE, nf = 2) {
    if (!all(unlist(lapply(df, is.factor)))) 
        stop("All variables must be factors")
    df <- as.data.frame(df)
    X <- acm.disjonctif(df)
    lig <- nrow(X)
    col <- ncol(X)
    var <- ncol(df)
    if (length(row.w) != lig) 
        stop("Non convenient row weights")
    if (any(row.w < 0)) 
        stop("row weight < 0")
    row.w <- row.w/sum(row.w)
    col.w <- apply(X, 2, function(x) sum(x*row.w))
    if (any(col.w == 0)) 
        stop("One category with null weight")
    X <- t(t(X)/col.w) - 1
    col.w <- col.w/var
    X <- as.dudi(data.frame(X), col.w, row.w, scannf = scannf, 
        nf = nf, call = match.call(), type = "acm")
    rcor <- matrix(0, ncol(df), X$nf)
    rcor <- row(rcor) + 0 + (0+1i) * col(rcor)
    floc <- function(x) {
        i <- Re(x)
        j <- Im(x)
        x <- X$l1[, j] * X$lw
        qual <- df[, i]
        poicla <- unlist(tapply(X$lw, qual, sum))
        z <- unlist(tapply(x, qual, sum))/poicla
        return(sum(poicla * z * z))
    }
    rcor <- apply(rcor, c(1, 2), floc)
    rcor <- data.frame(rcor)
    row.names(rcor) <- names(df)
    names(rcor) <- names(X$l1)
    X$cr <- rcor
    return(X)
}

"boxplot.acm" <- function (x, xax = 1, ...) {
    # correction d'un bug par P. Cornillon 29/10/2004
    if (!inherits(x, "acm")) 
        stop("Object of class 'acm' expected")
    if ((xax < 1) || (xax > x$nf)) 
        stop("non convenient axe number")
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    oritab <- eval.parent(as.list(x$call)[[2]])
    nvar <- ncol(oritab)
    if (nvar <= 7) 
        sco.boxplot(x$l1[, xax], oritab[, 1:nvar], clabel = 1)
    else if (nvar <= 14) {
        par(mfrow = c(1, 2))
        sco.boxplot(x$l1[, xax], oritab[, 1:(nvar%/%2)], clabel = 1.3)
        sco.boxplot(x$l1[, xax], oritab[, (nvar%/%2 + 1):nvar], 
            clabel = 1.3)
    }
    else {
        par(mfrow = c(1, 3))
        if ((a0 <- nvar%/%3) < nvar/3) 
            a0 <- a0 + 1
        sco.boxplot(x$l1[, xax], oritab[, 1:a0], clabel = 1.6)
        sco.boxplot(x$l1[, xax], oritab[, (a0 + 1):(2 * a0)], 
            clabel = 1.6)
        sco.boxplot(x$l1[, xax], oritab[, (2 * a0 + 1):nvar], 
            clabel = 1.6)
    }
}

"acm.burt" <- function (df1, df2, counts = rep(1, nrow(df1))) {
    if (!all(unlist(lapply(df1, is.factor)))) 
        stop("All variables must be factors")
    if (!all(unlist(lapply(df2, is.factor)))) 
        stop("All variables must be factors")
    if (nrow(df1) != nrow(df2)) 
        stop("non convenient row numbers")
    if (length(counts) != nrow(df2)) 
        stop("non convenient row numbers")
    g1 <- acm.disjonctif(df1)
    g1 <- g1 * counts
    g2 <- acm.disjonctif(df2)
    burt <- as.matrix(t(g1)) %*% as.matrix(g2)
    burt <- data.frame(burt)
    names(burt) <- names(g2)
    row.names(burt) <- names(g1)
    return(burt)
} 

"acm.disjonctif" <- function (df) {
    acm.util.df <- function(i) {
        cl <- df[,i]
        cha <- names(df)[i] 
        n <- length(cl)
        cl <- as.factor(cl)
        x <- matrix(0, n, length(levels(cl)))
        x[(1:n) + n * (unclass(cl) - 1)] <- 1
        dimnames(x) <- list(row.names(df), paste(cha,levels(cl),sep="."))
        return(x)
    }
    G <- lapply(1:ncol(df), acm.util.df)
    G <- data.frame (G, check.names = FALSE)
    return(G)
}


fac2disj<- function(fac, drop = FALSE) {
  ## Returns the disjunctive table corrseponding to a factor
  n <- length(fac)
  fac <- as.factor(fac)
  if(drop)
    fac <- factor(fac)
  x <- matrix(0, n, nlevels(fac))
  x[(1:n) + n * (unclass(fac) - 1)] <- 1
  dimnames(x) <- list(names(fac), as.character(levels(fac)))
  return(data.frame(x, check.names = FALSE))
}


scatterutil.eti <- function (x, y, label, clabel, boxes = TRUE, coul = rep(1, length(x)), 
    horizontal = TRUE, bg = "white") 
{
    if (length(label) == 0) 
        return(invisible())
    if (is.null(label)) 
        return(invisible())
    if (any(label == "")) 
        return(invisible())
    cex0 <- par("cex") * clabel
    for (i in 1:(length(x))) {
        cha <- as.character(label[i])
        cha <- paste(" ", cha, " ", sep = "")
        x1 <- x[i]
        y1 <- y[i]
        xh <- strwidth(cha, cex = cex0)
        yh <- strheight(cha, cex = cex0) * 5/3
        if (!horizontal) {
            tmp <- scatterutil.convrot90(xh, yh)
            xh <- tmp[1]
            yh <- tmp[2]
        }
        if (boxes) {
            rect(x1 - xh/2, y1 - yh/2, x1 + xh/2, y1 + yh/2, 
                col = bg, border = coul[i])
        }
        if (horizontal) {
            text(x1, y1, cha, cex = cex0, col = coul[i])
        }
        else {
            text(x1, y1, cha, cex = cex0, col = coul[i], srt = 90)
        }
    }
}

scatterutil.sub <- function (cha, csub, possub = "bottomleft") 
{
    cha <- as.character(cha)
    if (length(cha) == 0) 
        return(invisible())
    if (is.null(cha)) 
        return(invisible())
    if (is.na(cha)) 
        return(invisible())
    if (any(cha == "")) 
        return(invisible())
    if (csub == 0) 
        return(invisible())
    cex0 <- par("cex") * csub
    cha <- paste(" ", cha, " ", sep = "")
    xh <- strwidth(cha, cex = cex0)
    yh <- strheight(cha, cex = cex0) * 5/3
    if (possub == "bottomleft") {
        x1 <- par("usr")[1]
        y1 <- par("usr")[3]
        rect(x1, y1, x1 + xh, y1 + yh, col = "white", border = 0)
        text(x1 + xh/2, y1 + yh/2, cha, cex = cex0)
    }
    else if (possub == "topleft") {
        x1 <- par("usr")[1]
        y1 <- par("usr")[4]
        rect(x1, y1, x1 + xh, y1 - yh, col = "white", border = 0)
        text(x1 + xh/2, y1 - yh/2, cha, cex = cex0)
    }
    else if (possub == "bottomright") {
        x1 <- par("usr")[2]
        y1 <- par("usr")[3]
        rect(x1, y1, x1 - xh, y1 + yh, col = "white", border = 0)
        text(x1 - xh/2, y1 + yh/2, cha, cex = cex0)
    }
    else if (possub == "topright") {
        x1 <- par("usr")[2]
        y1 <- par("usr")[4]
        rect(x1, y1, x1 - xh, y1 - yh, col = "white", border = 0)
        text(x1 - xh/2, y1 - yh/2, cha, cex = cex0)
    }
}

scatterutil.grid <- function (cgrid) 
{
    col <- "lightgray"
    lty <- 1
    xaxp <- par("xaxp")
    ax <- (xaxp[2] - xaxp[1])/xaxp[3]
    yaxp <- par("yaxp")
    ay <- (yaxp[2] - yaxp[1])/yaxp[3]
    a <- min(ax, ay)
    v0 <- seq(xaxp[1], xaxp[2], by = a)
    h0 <- seq(yaxp[1], yaxp[2], by = a)
    abline(v = v0, col = col, lty = lty)
    abline(h = h0, col = col, lty = lty)
    if (cgrid <= 0) 
        return(invisible())
    cha <- paste(" d = ", a, " ", sep = "")
    cex0 <- par("cex") * cgrid
    xh <- strwidth(cha, cex = cex0)
    yh <- strheight(cha, cex = cex0) * 5/3
    x1 <- par("usr")[2]
    y1 <- par("usr")[4]
    rect(x1 - xh, y1 - yh, x1 + xh, y1 + yh, col = "white", border = 0)
    text(x1 - xh/2, y1 - yh/2, cha, cex = cex0)
}

scatterutil.base <- function (dfxy, xax, yax, xlim, ylim, grid, addaxes, cgrid, include.origin, 
    origin, sub, csub, possub, pixmap, contour, area, add.plot) 
{
    df <- data.frame(dfxy)
    if (!is.data.frame(df)) 
        stop("Non convenient selection for df")
    if ((xax < 1) || (xax > ncol(df))) 
        stop("Non convenient selection for xax")
    if ((yax < 1) || (yax > ncol(df))) 
        stop("Non convenient selection for yax")
    x <- df[, xax]
    y <- df[, yax]
    if (is.null(xlim)) {
        x1 <- x
        if (include.origin) 
            x1 <- c(x1, origin[1])
        x1 <- c(x1 - diff(range(x1)/10), x1 + diff(range(x1))/10)
        xlim <- range(x1)
    }
    if (is.null(ylim)) {
        y1 <- y
        if (include.origin) 
            y1 <- c(y1, origin[2])
        y1 <- c(y1 - diff(range(y1)/10), y1 + diff(range(y1))/10)
        ylim <- range(y1)
    }
    if (!is.null(pixmap)) {
        if (is.null(class(pixmap))) 
            pixmap <- NULL
        if (is.na(charmatch("pixmap", class(pixmap)))) 
            pixmap <- NULL
    }
    if (!is.null(contour)) {
        if (!is.data.frame(contour)) 
            contour <- NULL
        if (ncol(contour) != 4) 
            contour <- NULL
    }
    if (!is.null(area)) {
        if (!is.data.frame(area)) 
            area <- NULL
        if (!is.factor(area[, 1])) 
            area <- NULL
        if (ncol(area) < 3) 
            area <- NULL
    }
    if (!add.plot) 
        plot.default(0, 0, type = "n", asp = 1, xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
            xaxs = "i", yaxs = "i", frame.plot = FALSE)
    if (!is.null(pixmap)) {
        plot(pixmap, add = TRUE)
    }
    if (!is.null(contour)) {
        apply(contour, 1, function(x) segments(x[1], x[2], x[3], 
            x[4], lwd = 1))
    }
    if (grid & !add.plot) 
        scatterutil.grid(cgrid)
    if (addaxes & !add.plot) 
        abline(h = 0, v = 0, lty = 1, lwd=0.4)
    if (!is.null(area)) {
        nlev <- nlevels(area[, 1])
        x1 <- area[, 2]
        x2 <- area[, 3]
        for (i in 1:nlev) {
            lev <- levels(area[, 1])[i]
            a1 <- x1[area[, 1] == lev]
            a2 <- x2[area[, 1] == lev]
            polygon(a1, a2)
        }
    }
    if (csub > 0) 
        scatterutil.sub(sub, csub, possub)
    return(list(x = x, y = y))
}

scatterutil.star <- function (x, y, z, cstar, coul = rep(1, length(x)), star.lwd=0.5) 
{
    z <- z/sum(z)
    x1 <- sum(x * z)
    y1 <- sum(y * z)
    for (i in which(z > 0)) {
        hx <- cstar * (x[i] - x1)
        hy <- cstar * (y[i] - y1)
        segments(x1, y1, x1 + hx, y1 + hy, col = coul, lwd=star.lwd)
    }
}

scatterutil.ellipse <- function (x, y, z, cellipse, axesell, coul = rep(1, length(x)), ellipse.lwd=2, axesell.lwd=1, star.lwd=0.5) 
{
    if (any(is.na(z))) 
        return(invisible())
    if (sum(z * z) == 0) 
        return(invisible())
    util.ellipse <- function(mx, my, vx, cxy, vy, coeff) {
        lig <- 100
        epsi <- 1e-10
        x <- 0
        y <- 0
        if (vx < 0) 
            vx <- 0
        if (vy < 0) 
            vy <- 0
        if (vx == 0 && vy == 0) 
            return(NULL)
        delta <- (vx - vy) * (vx - vy) + 4 * cxy * cxy
        delta <- sqrt(delta)
        l1 <- (vx + vy + delta)/2
        l2 <- vx + vy - l1
        if (l1 < 0) 
            l1 <- 0
        if (l2 < 0) 
            l2 <- 0
        l1 <- sqrt(l1)
        l2 <- sqrt(l2)
        test <- 0
        if (vx == 0) {
            a0 <- 0
            b0 <- 1
            test <- 1
        }
        if ((vy == 0) && (test == 0)) {
            a0 <- 1
            b0 <- 0
            test <- 1
        }
        if (((abs(cxy)) < epsi) && (test == 0)) {
            a0 <- 1
            b0 <- 0
            test <- 1
        }
        if (test == 0) {
            a0 <- 1
            b0 <- (l1 * l1 - vx)/cxy
            norm <- sqrt(a0 * a0 + b0 * b0)
            a0 <- a0/norm
            b0 <- b0/norm
        }
        a1 <- 2 * pi/lig
        c11 <- coeff * a0 * l1
        c12 <- (-coeff) * b0 * l2
        c21 <- coeff * b0 * l1
        c22 <- coeff * a0 * l2
        angle <- 0
        for (i in 1:lig) {
            cosinus <- cos(angle)
            sinus <- sin(angle)
            x[i] <- mx + c11 * cosinus + c12 * sinus
            y[i] <- my + c21 * cosinus + c22 * sinus
            angle <- angle + a1
        }
        return(list(x = x, y = y, seg1 = c(mx + c11, my + c21, 
            mx - c11, my - c21), seg2 = c(mx + c12, my + c22, 
            mx - c12, my - c22)))
    }
    z <- z/sum(z)
    m1 <- sum(x * z)
    m2 <- sum(y * z)
    v1 <- sum((x - m1) * (x - m1) * z)
    v2 <- sum((y - m2) * (y - m2) * z)
    cxy <- sum((x - m1) * (y - m2) * z)
    ell <- util.ellipse(m1, m2, v1, cxy, v2, cellipse)
    if (is.null(ell)) 
        return(invisible())
    polygon(ell$x, ell$y, border = coul, lwd=ellipse.lwd)
    if (axesell) 
        segments(ell$seg1[1], ell$seg1[2], ell$seg1[3], ell$seg1[4], 
            lty = 2, col = coul, lwd=axesell.lwd)
    if (axesell) 
        segments(ell$seg2[1], ell$seg2[2], ell$seg2[3], ell$seg2[4], 
            lty = 2, col = coul, lwd=axesell.lwd)
}

s.class2 <- function (dfxy, fac, wt = rep(1, length(fac)), xax = 1, yax = 2, 
    cstar = 1, cellipse = 1.5, axesell = TRUE, label = levels(fac), 
    clabel = 1, cpoint = 1, pch = 20, col = rep(1, length(levels(fac))), 
    xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, origin = c(0, 
        0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", 
    cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE,
	ellipse.lwd=2, axesell.lwd=1, star.lwd=0.5, axes=TRUE, xlab='NMDS1', ylab='NMDS2')
{
    opar <- par(mar = par("mar"))
	if (axes) {
		par(mar = c(0.1, 0.1, 0.1, 0.1) * 40)
	} else {
		par(mar = c(0.1, 0.1, 0.1, 0.1))
	}
    on.exit(par(opar))
    dfxy <- data.frame(dfxy)
    if (!is.data.frame(dfxy)) 
        stop("Non convenient selection for dfxy")
    if (any(is.na(dfxy))) 
        stop("NA non implemented")
    if (!is.factor(fac)) 
        stop("factor expected for fac")
    dfdistri <- fac2disj(fac) * wt
    coul <- col
    w1 <- unlist(lapply(dfdistri, sum))
    dfdistri <- t(t(dfdistri)/w1)
    coox <- as.matrix(t(dfdistri)) %*% dfxy[, xax]
    cooy <- as.matrix(t(dfdistri)) %*% dfxy[, yax]
    if (nrow(dfxy) != nrow(dfdistri)) 
        stop(paste("Non equal row numbers", nrow(dfxy), nrow(dfdistri)))
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
        contour = contour, area = area, add.plot = add.plot)
    if (cpoint > 0) 
        for (i in 1:ncol(dfdistri)) {
            pch <- rep(pch, length = nrow(dfxy))
            points(coo$x[dfdistri[, i] > 0], coo$y[dfdistri[, 
                i] > 0], pch = pch[dfdistri[, i] > 0], cex = par("cex") * 
                cpoint, col = coul[i])
        }
    if (cstar > 0) 
        for (i in 1:ncol(dfdistri)) {
            scatterutil.star(coo$x, coo$y, dfdistri[, i], cstar = cstar, 
                coul[i], star.lwd=star.lwd)
        }
    if (cellipse > 0) 
        for (i in 1:ncol(dfdistri)) {
            scatterutil.ellipse(coo$x, coo$y, dfdistri[, i], 
                cellipse = cellipse, axesell = axesell, coul[i], ellipse.lwd=ellipse.lwd, axesell.lwd=axesell.lwd)
        }
    if (clabel > 0) 
        scatterutil.eti(coox, cooy, label, clabel, coul = col)
    box(lwd=0.4)
	if (axes) {
		title(xlab=xlab, ylab=ylab)
		axis(side=1, lwd=0.4)	
		axis(side=2, las=1, lwd=0.4)
	}
    invisible(match.call())
}

