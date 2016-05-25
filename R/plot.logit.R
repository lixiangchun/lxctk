
plot.logit <- function(logitObj=NULL, remove.intercept=TRUE, coef.names=NULL, main=NULL, main.xycoord=c(0.3,0.7),...)
{
	#require("grid") || stop("`grid' package not found")
    #require("rmeta") || stop("`rmeta' package not found")
	XLAB="Risk ratio (95% CI)"
		r <- summary(logitObj)
		if (!is.null(coef.names)) {
			variable <- coef.names
		} else {
			variable <- rownames(r$coefficients)
		}
		CI <- exp(confint(logitObj))
		lower.95 <- CI[,1]
		upper.95 <- CI[,2]
		HR <- exp(r$coefficients[,1])
		p.values <- r$coefficients[,4]

		p.values <- sprintf("%.3g", p.values)
		HR.95.CI <- sprintf("%.3g (%.3g ~ %.3g)", HR, lower.95, upper.95)

		labeltext <- cbind(variable, HR.95.CI, p.values)	
		labeltext <- rbind(c('Variable', "RR (95% CI)", "P-value"), labeltext)

		is.summary = c(1, rep(0, 20))		
		
		if (remove.intercept) {
			#forestplot(labeltext[-2,],mean=c(NA,as.numeric(HR[-1])),lower=c(NA,as.numeric(lower.95[-1])),upper=c(NA, as.numeric(upper.95[-1])), zero=TRUE, col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"), boxsize=0.5, is.summary=is.summary, align=c('l','l','l'), clip=clip, xlab=XLAB)
			forestplot(labeltext[-2,],mean=c(NA,as.numeric(HR[-1])),lower=c(NA,as.numeric(lower.95[-1])),upper=c(NA, as.numeric(upper.95[-1])), zero=TRUE, col=meta.colors(box="black",line="black", summary="black", zero='black'), boxsize=0.5, is.summary=is.summary, align=c('l','l','l'), xlab=XLAB,...)
		} else {
			#forestplot(labeltext,mean=c(NA,as.numeric(HR)),lower=c(NA,as.numeric(lower.95)),upper=c(NA, as.numeric(upper.95)), zero=TRUE, col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"), boxsize=0.5, is.summary=is.summary, align=c('l','l','l'), xlab=XLAB,...)
			forestplot(labeltext,mean=c(NA,as.numeric(HR)),lower=c(NA,as.numeric(lower.95)),upper=c(NA, as.numeric(upper.95)), zero=TRUE, col=meta.colors(all.elements='black'), boxsize=0.5, is.summary=is.summary, align=c('l','l','l'), xlab=XLAB,...)
		}	
		if (!is.null(main)) {
			x <- main.xycoord[1]
			y <- main.xycoord[2]
			text(x, y, main, xpd=TRUE)
		}
}


