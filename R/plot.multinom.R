	
cat.RR <- function(multinomObj)
{
	data <- summary(multinomObj)$coefficients
	labels <- rownames(data)
	x0 <- data[1,]
	out <- list()
	out[[labels[1]]] <- x0 
	n <- nrow(data)
	m <- ncol(data)
	for (i in 2:n) {
		out[[labels[[i]]]] <- data[i, 2:m]	
	}
	coeff <- do.call('c', out)

	data <- summary(multinomObj)$standard.errors
	labels <- rownames(data)
	x0 <- data[1,]
	out <- list()
	out[[labels[1]]] <- x0 
	n <- nrow(data)
	m <- ncol(data)
	for (i in 2:n) {
		out[[labels[[i]]]] <- data[i, 2:m]	
	}
	std.err <- do.call('c', out)
	z <- coeff / std.err	
	p.values <- (1 - pnorm(abs(z), 0, 1)) * 2
	list(coeff=coeff, std.err=std.err, Ward.ratio=z, p.values=p.values)	
}

cat.CI <- function(multinomObj)
{
	data <- confint(multinomObj)
	data <- as.data.frame(data)
	labels <- rownames(summary(multinomObj)$coefficients)
	out <- list()
	n <- nrow(data)
	m <- ncol(data)
	end <- seq(2, m, by=2)
	beg <- end - 1
	k <- m / 2
	for (i in 1:k) {
		if (i > 1) {
			li <- data[2:n, c(beg[i], end[i])]	
		} else {
			li <- data[1:n, c(beg[i], end[i])]	
		}
		colnames(li) <- c('2.5%(CI)','97.5%(CI)')
		rownames(li) <- paste(labels[i], rownames(li), sep=":")
		out[[i]] <- li
	}
	do.call('rbind', out)
}



plot.multinom <- function(multinomObj, plot=TRUE, remove.intercept=TRUE, coef.names=NULL,main=NULL, main.xycoord=c(0.4, 1),...)
{
	r <- summary(multinomObj)	
	## http://www.cognopod.com/RADLIBS/radlib.php?r=40
	## calculate p values using Wald tests
	coeff <- r$coefficients
	if (is.null(nrow(coeff))) { # dependent variable has 2 choices, i.e. dummy variable
		variable <- multinomObj$coefnames
		coeff <- r$coefficients
		std.err <- r$standard.errors
		z <- coeff / std.err
		p.values <- (1 - pnorm(abs(z), 0, 1)) * 2
		x <- list(coeff=coeff, std.err=std.err, Ward.ratio=z, p.values=p.values)
		x <- do.call('cbind', x)
		RR <- exp(coeff)
		CI <- as.data.frame(exp(confint(multinomObj)))
	} else {
		obj <- cat.RR(multinomObj)
		z <- obj$Ward.ratio
		p.values <- obj$p.values
		x <- list(coeff=obj$coeff, std.err=obj$std.err, Ward.ratio=obj$Ward.ratio, p.values=p.values)
		x <- do.call('cbind', x)
		RR <- exp(obj$coeff)
		CI <- exp(cat.CI(multinomObj))
		variable <- rownames(CI)
		remove.intercept <- TRUE
	}
	##colnames(CI) <- c('2.5%(CI)', '97.5%(CI)')

	if (plot) {
		if (!is.null(coef.names)) {
			variable <- c('',coef.names)
		}
		HR <- RR
        lower.95 <- CI[,1]
        upper.95 <- CI[,2]
		##HR.95.CI <- paste(paste("(",paste(lower.95, upper.95, sep=" ~ "), sep=""),")",sep="")  
		##HR.95.CI <- paste(HR, HR.95.CI, sep=" ")
		HR.95.CI <- sprintf("%.3g (%.3g ~ %.3g)", HR, lower.95, upper.95)
		labeltext <-  cbind(variable, HR.95.CI, p.values=sprintf('%.3g', p.values)) 
		labeltext <- rbind(c('Variable', "RR (95% CI)", "P-value"), labeltext)
        is.summary <- c(1, rep(0, 200))   
		XLAB="Risk ratio (95% CI)"	
		if (remove.intercept) {
			forestplot(labeltext[-2,],mean=c(NA,as.numeric(HR[-1])),lower=c(NA,as.numeric(lower.95[-1])),upper=c(NA, as.numeric(upper.95[-1])), zero=TRUE, col=meta.colors(box="black",line="black", summary="black", zero='black'), boxsize=0.5, is.summary=is.summary, align=c('l','l','l'), xlab=XLAB, ...)
        } else {
			forestplot(labeltext,mean=c(NA,as.numeric(HR)),lower=c(NA,as.numeric(lower.95)),upper=c(NA, as.numeric(upper.95)), zero=TRUE, col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"), boxsize=0.5, is.summary=is.summary, align=c('l','l','l'), xlab=XLAB, ...)
        }  
		if (!is.null(main)) {
			text(main.xycoord[1], main.xycoord[2], main, xpd=TRUE)
		}
	}		
	res <- cbind(x, RR, CI)
	invisible(res)
}

