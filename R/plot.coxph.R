
plot.coxph <- function(coxTable=NULL, coxObj=NULL, clip=c(0,8), variable.display=NULL, ...)
{
	require("grid") || stop("`grid' package not found")
    require("rmeta") || stop("`rmeta' package not found")
	XLAB="Hazard Ratio (95% CI)"
	if (!is.null(coxTable)) {
		tabletext=coxTable[,1:3]		
		m=suppressWarnings( as.numeric( coxTable[,4] ) )
		l=suppressWarnings( as.numeric( coxTable[,5] ) )
		u=suppressWarnings( as.numeric( coxTable[,6] ) )
		is.summary <- try(as.integer(coxTable[,7]), silent=TRUE)
		if (class(is.summary) == 'try-error') {
			is.summary <- rep(0, 20)
			is.summary[1] <- 1
		}
		#forestplot(tabletext,as.numeric(m),as.numeric(l),as.numeric(u), align=c('l','c','c','c'), clip=c(0,8), zero=TRUE, xlab=XLAB, is.summary=is.summary,col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"), boxsize=0.5)
		forestplot(tabletext,as.numeric(m),as.numeric(l),as.numeric(u), align=c('l','c','c','c'), clip=clip, zero=TRUE, xlab=XLAB, is.summary=is.summary,col=meta.colors(all.elements='black'), boxsize=0.5, ...)
	}
	else if (!is.null(coxObj)) {
		r=summary(coxObj)
		if (is.null(variable.display)) {
			variable.display=rownames(r$conf.int)
		}
		HR <- r$conf.int[,1]
		lower.95 <- r$conf.int[,3]
		upper.95 <- r$conf.int[,4]
		p.values <- r$coefficients[,5]

		#HR=sprintf("%.3g", HR)
		#lower.95=sprintf("%.3g",lower.95)
		#upper.95=sprintf("%.3g",upper.95)

		#HR.95.CI=paste(paste("(",paste(lower.95, upper.95, sep=" ~ "), sep=""),")",sep="")
		#HR.95.CI=paste(HR, HR.95.CI, sep=" ")
		HR.95.CI <- sprintf("%.3g (%.3g ~ %.3g)", HR, lower.95, upper.95)
		p.values <- sprintf("%.3g", p.values)

		labeltext=cbind(variable.display, HR.95.CI, p.values)	
		labeltext=rbind(c('Variable', "HR (95% CI)", "P-value"), labeltext)

		is.summary = c(1, rep(0, 20))		

		#forestplot(labeltext,mean=c(NA,as.numeric(HR)),lower=c(NA,as.numeric(lower.95)),upper=c(NA, as.numeric(upper.95)), zero=TRUE, col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"), boxsize=0.5, is.summary=is.summary, align=c('l','l','l'), clip=clip, xlab=XLAB)
		forestplot(labeltext,mean=c(NA,as.numeric(HR)),lower=c(NA,as.numeric(lower.95)),upper=c(NA, as.numeric(upper.95)), zero=TRUE, col=meta.colors(all.elements='black'), boxsize=0.5, is.summary=is.summary, align=c('l','l','l'), clip=clip, xlab=XLAB, ...)
	}
}
