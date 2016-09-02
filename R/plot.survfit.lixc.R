plot.survfit.lixc <- function(survdiff.formula, data, legend.labels=c('Wild-type','Mutant'), ci.lab.xycoord=c(70,0.2), legend.xycoord=c(40,0.4), col=c("cyan4", "brown"),no.ci.lab=FALSE,pngfig=NULL,...)
{
	diff <- try(survdiff(survdiff.formula, data=data))
	if (class(diff) == 'try-error') next;
	surv.p.value <- pchisq(diff$chisq, 1, lower.tail=F)

	#x <- data[, gene]
	#x <- x[!is.na(x)]
	#mutant_count = sum(x > 0)
	#wild_type_count = sum(x == 0)
	mutant_count <- diff$n[2]
	wild_type_count <- diff$n[1]

	###################CoxPH##########################
	r <- coxph(survdiff.formula, data=data)
	ci <- exp(confint(r))
	HR <- exp(r$coefficients)
	cox.p.value <- summary(r)$coefficients[5]
	lower <- ci[1]
	upper <- ci[2]

	if (is.infinite(ci[2]))
		texts <- sprintf("HR=%.3g (%s CI: NA); p=%.3g", HR, "95%", surv.p.value)
	else
		texts <- sprintf("HR=%.3g (%s CI: %.3g ~ %.3g); p=%.3g", HR, "95%", lower, upper, surv.p.value)
	if (no.ci.lab)
		texts <- sprintf("Log-rank test, p=%.2g", surv.p.value)

	####################KM-Curve#######################
	if (!is.null(pngfig)) {
		png(pngfig, width=2000, height=2000, res=300, family="Arials")
	}
	fit <- survfit(survdiff.formula, data=data)
	#plot(fit, col=c("cyan4", "brown"), xlab="", axes=FALSE, lwd=line.lwd, lty=c(1, 1), cex.main=2, font.main=2, xaxs='i',yaxs='i', main=main)
	plot(fit, col=col, xlab="", axes=FALSE, lty=c(1, 1), xaxs='i',yaxs='i',...)
	axis(side=1, padj=0)
	axis(side=2, at=seq(0,1, 0.1), las=2)
	mtext("Time since operation (months)", side=1, padj=3, cex=1.5, font=1)
	mtext("Overall survival (%)", side=2, padj=-3, cex=1.5)
	text(ci.lab.xycoord[1], ci.lab.xycoord[2], label=texts, font=1, cex=1.2)

	#legend.label1=paste(legend.labels[1], " (events ", mutant_count, ")", sep="")
	#legend.label2=paste(legend.labels[2], " (events ", wild_type_count, ")", sep="")
	#legend(legend.xycoord[1], legend.xycoord[2], legend=c(legend.label1, legend.label2), fill=c("brown", "cyan4"), box.lwd=0, box.col='white', cex=1.2, border='white',x.intersp=0.2)
	
	if (length(diff$n)>length(legend.labels) || length(diff$n)>length(col)) {
		warning('Event number is not equal to number of user provided legend labels or colors for KM curve.')
	}
	legend.text <- sprintf("%s (n=%g)", legend.labels, diff$n)
	legend(legend.xycoord[1], legend.xycoord[2], legend=legend.text, fill=col, box.lwd=0, box.col=NA, cex=1.2, border=col,x.intersp=0.2)
	if (!is.null(pngfig)) {
		dev.off()
	}
	invisible(list(surv.diff=diff, coxph.res=c))
}

