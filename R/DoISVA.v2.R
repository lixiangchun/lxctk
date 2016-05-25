fastICA.v2 <- function (X, n.comp, alg.typ = c("parallel", "deflation"), fun = c("logcosh", 
    "exp"), alpha = 1, method = c("R", "C"), row.norm = FALSE, 
    maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL) 
{
    dd <- dim(X)
    d <- dd[dd != 1L]
    if (length(d) != 2L) 
        stop("data must be matrix-conformal")
    X <- if (length(d) != length(dd)) 
        matrix(X, d[1L], d[2L])
    else as.matrix(X)
    if (alpha < 1 || alpha > 2) 
        stop("alpha must be in range [1,2]")
    method <- match.arg(method)
    alg.typ <- match.arg(alg.typ)
    fun <- match.arg(fun)
    n <- nrow(X)
    p <- ncol(X)
    if (n.comp > min(n, p)) {
        message("'n.comp' is too large: reset to ", min(n, p))
        n.comp <- min(n, p)
    }
    if (is.null(w.init)) 
        w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
    else {
        if (!is.matrix(w.init) || length(w.init) != (n.comp^2)) 
            stop("w.init is not a matrix or is the wrong size")
    }
    if (method == "R") {
        if (verbose) 
            message("Centering")
        X <- scale(X, scale = FALSE)
        X <- if (row.norm) 
            t(scale(X, scale = row.norm))
        else t(X)
        if (verbose) 
            message("Whitening")

		#BEGIN: Added by Xiangchun Li
		X.colname=colnames(X)
		if (any(is.na(X.colname))) {
			message("Removing NA value from matrix X, script modified by Li Xiangchun.")
			I=which(!is.na(X.colname))
			X=X[,I]	
		}
		#END: Added by Xiangchun Li

        V <- X %*% t(X)/n
        s <- La.svd(V)
        D <- diag(c(1/sqrt(s$d)))
        K <- D %*% t(s$u)
        K <- matrix(K[1:n.comp, ], n.comp, p)
        X1 <- K %*% X
        a <- if (alg.typ == "deflation") 
            ica.R.def(X1, n.comp, tol = tol, fun = fun, alpha = alpha, 
                maxit = maxit, verbose = verbose, w.init = w.init)
        else if (alg.typ == "parallel") 
            ica.R.par(X1, n.comp, tol = tol, fun = fun, alpha = alpha, 
                maxit = maxit, verbose = verbose, w.init = w.init)
        w <- a %*% K
        S <- w %*% X
        A <- t(w) %*% solve(w %*% t(w))
        return(list(X = t(X), K = t(K), W = t(a), A = t(A), S = t(S)))
    }
    else if (method == "C") {
        a <- .C(icainc_JM, as.double(X), as.double(w.init), as.integer(p), 
            as.integer(n), as.integer(n.comp), as.double(alpha), 
            as.integer(1), as.integer(row.norm), as.integer(1L + 
                (fun == "exp")), as.integer(maxit), as.double(tol), 
            as.integer(alg.typ != "parallel"), as.integer(verbose), 
            X = double(p * n), K = double(n.comp * p), W = double(n.comp * 
                n.comp), A = double(p * n.comp), S = double(n.comp * 
                n))
        X1 <- matrix(a$X, n, p)
        K <- matrix(a$K, p, n.comp)
        W <- matrix(a$W, n.comp, n.comp)
        A <- matrix(a$A, n.comp, p)
        S <- matrix(a$S, n, n.comp)
        list(X = X1, K = K, W = W, A = A, S = S)
    }
}


isvaFn.v2 <- function (data.m, pheno.v, ncomp = NULL) 
{
    lm.o <- lm(t(data.m) ~ pheno.v)
    res.m <- t(lm.o$res)
    model <- model.matrix(~1 + pheno.v)
    if (is.null(ncomp)) {
        rmt.o <- EstDimRMT(res.m)
        ncomp <- rmt.o$dim
        print(paste("Number of candidate ISVs = ", ncomp, sep = ""))
    }
    else {
        print("no need to estimate dimensionality")
    }
    fICA.o <- fastICA.v2(res.m, n.comp = ncomp)
    tmp.m <- t(fICA.o$A)
    isv.m <- tmp.m
    sd <- 1/sqrt(ncol(data.m) - 3)
    for (k in 1:ncol(tmp.m)) {
        cor.v <- as.vector(cor(t(data.m), tmp.m[, k]))
        z.v <- 0.5 * log((1 + cor.v)/(1 - cor.v))
        pv.v <- 2 * pnorm(abs(z.v), 0, sd, lower.tail = FALSE)
        tmp.s <- sort(pv.v, decreasing = FALSE, index.return = TRUE)
        #qv.o <- qvalue(pv.v)             # NA value will induce ERROR, commented by Li Xiangchun

		## BEGIN: Added by Li Xiangchun
		if (any(is.na(pv.v))) {
			message('NA p-value found, remove them to get script running; modified by Li Xiangchun.')
			qv.o <- qvalue(na.omit(pv.v))
		}
		## END: Added by Li Xiangchun		

        nsig <- length(which(qv.o$qvalues < 0.05))
        if (nsig < 500) {
            nsig <- 500
        }
        red.m <- data.m[tmp.s$ix[1:nsig], ]
        fICA.o <- fastICA.v2(red.m, n.comp = ncomp)
        cor.v <- abs(cor(tmp.m[, k], t(fICA.o$A)))
        kmax <- which.max(cor.v)
        isv.m[, k] <- t(fICA.o$A)[, kmax]
        print(paste("Built ISV ", k, sep = ""))
    }
    return(list(n.isv = ncomp, isv = isv.m))
}


DoISVA.v2 = function (data.m, pheno.v, cf.m = NULL, factor.log, pvthCF = 0.01, 
    th = 0.05, ncomp = NULL) 
{
    isva.o <- isvaFn.v2(data.m, pheno.v, ncomp)
    if (is.null(cf.m) == FALSE) {
        tmp.m <- cbind(pheno.v, cf.m)
        treatfactor <- c(FALSE, factor.log)
        pv.m <- matrix(nrow = ncol(isva.o$isv), ncol = 1 + ncol(cf.m))
        colnames(pv.m) <- c("POI", colnames(cf.m))
        for (c in 1:ncol(tmp.m)) {
            if (treatfactor[c] == FALSE) {
                for (sv in 1:ncol(isva.o$isv)) {
                  lm.o <- lm(isva.o$isv[, sv] ~ as.numeric(tmp.m[, 
                    c]))
                  pv.m[sv, c] <- summary(lm.o)$coeff[2, 4]
                }
            }
            else {
                for (sv in 1:ncol(isva.o$isv)) {
                  lm.o <- lm(isva.o$isv[, sv] ~ as.factor(tmp.m[, 
                    c]))
                  pv.m[sv, c] <- pf(summary(lm.o)$fstat[1], summary(lm.o)$fstat[2], 
                    summary(lm.o)$fstat[3], lower.tail = FALSE)
                }
            }
        }
        print("Selecting ISVs")
        selisv.idx <- vector()
        for (sv in 1:nrow(pv.m)) {
            ncf <- length(which(pv.m[sv, 2:ncol(pv.m)] < pvthCF))
            minpv <- min(pv.m[sv, 2:ncol(pv.m)])
            phpv <- pv.m[sv, 1]
            if (ncf > 0) {
                if (minpv < phpv) {
                  selisv.idx <- c(selisv.idx, sv)
                }
            }
        }
        if (length(selisv.idx) == 0) {
            print("No ISVs selected because none correlated with the given confounders. Rerun ISVA with cf.m=NULL option")
            stop
        }
    }
    else {
        selisv.idx <- 1:ncol(isva.o$isv)
        pv.m <- NULL
    }
    print("Running final multivariate regressions with selected ISVs")
    selisv.m <- matrix(isva.o$isv[, selisv.idx], ncol = length(selisv.idx))
    mod <- model.matrix(~pheno.v + selisv.m)
    modNULL <- model.matrix(~selisv.m)
    df1 <- dim(mod)[2]
    df0 <- dim(modNULL)[2]
    pv.v <- rep(0, nrow(data.m))
    Id <- diag(ncol(data.m))
    resid <- data.m %*% (Id - mod %*% solve(t(mod) %*% mod) %*% 
        t(mod))
    rss1 <- rowSums(resid * resid)
    rm(resid)
    residNULL <- data.m %*% (Id - modNULL %*% solve(t(modNULL) %*% 
        modNULL) %*% t(modNULL))
    rssNULL <- rowSums(residNULL * residNULL)
    rm(residNULL)
    fstats <- ((rssNULL - rss1)/(df1 - df0))/(rss1/(ncol(data.m) - 
        df1))
    pv.v <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (ncol(data.m) - 
        df1))
    pv.s <- sort(pv.v, decreasing = FALSE, index.return = TRUE)
    qv.v <- qvalue(pv.s$x)$qvalue
    ntop <- length(which(qv.v < th))
    print(paste("Number of DEGs after ISV adjustment = ", ntop, 
        sep = ""))
    if (ntop > 0) {
        pred.idx <- pv.s$ix[1:ntop]
        lm.o <- lm(t(data.m[pred.idx, ]) ~ pheno.v + selisv.m)
        tstats.v <- unlist(lapply(summary(lm.o), function(x) {
            x$coeff[2, 3]
        }))
        lm.m <- cbind(tstats.v, pv.s$x[1:ntop], qv.v[1:ntop])
        colnames(lm.m) <- c("t-stat", "P-value", "q-value")
    }
    else {
        pred.idx <- NULL
        lm.m <- NULL
    }
    return(list(spv = pv.s$x, qv = qv.v, rk = pv.s$ix, ndeg = ntop, 
        deg = pred.idx, lm = lm.m, isv = selisv.m, nsv = length(selisv.idx), 
        pvCF = pv.m, selisv = selisv.idx))
}

