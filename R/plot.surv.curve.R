

plot.surv.curve <- function(gene, survdiff.formula, coxph.formula, data, legend.labels=c('Mutant', 'Wild-type'), ci.lab.xycoord=c(70,0.2), legend.xycoord=c(40,0.4), pngfig=NULL, main=gene, line.lwd=par("lwd"))
{
	diff = try(survdiff(survdiff.formula, data=data))
	if (class(diff) == 'try-error') next;
	p = pchisq(diff$chisq, 1, lower.tail=F)

	#mutant_count = sum(data[,gene] > 0)
	#wild_type_count = sum(data[,gene] == 0)

	x <- data[, gene]
	x <- x[!is.na(x)]
	mutant_count = sum(x > 0)
	wild_type_count = sum(x == 0)

	###################CoxPH##########################
	c = coxph(coxph.formula, data=data)
	r = summary(c)
	r.NR = nrow(r)
	r.NC = ncol(r)
	cox.pvalue = r$coefficients[r.NR, r.NC]

	ci = exp(confint(c))

	ci.NR = nrow(ci)
	signature.ci = ci[ci.NR,]
	HR = format( r$coefficients[ci.NR,][2], digit=2 )
	lower = format( signature.ci[1], digit=2 )
	upper = format( signature.ci[2], digit=2 )
	p = format(p, digit=3, scientific=TRUE)

	if (is.infinite(signature.ci[2])) {
		texts = paste('HR=', HR, " (95% CI: NA); p=", p, sep="")
	}
	else {
		texts = paste('HR=', HR, " (95% CI: ", lower, "-", upper, "); p=", p, sep="")
	}

	####################KM-Curve#######################
	if (!is.null(pngfig)) {
		png(pngfig, width=2000, height=2000, res=300, family="Arials")
	}
	fit = survfit(survdiff.formula, data=data)
	plot(fit, col=c("cyan4", "brown"), xlab="", axes=F, lwd=line.lwd, lty=c(1, 1), cex.main=2, font.main=2, xaxs='i',yaxs='i', main=main)
	#axis(side=1, at=seq(0, 120, 20), padj=0)
	axis(side=1, padj=0)
	axis(side=2, at=seq(0,1, 0.1), las=2)
	mtext("Time since operation (months)", side=1, padj=3, cex=1.5, font=1)
	mtext("Overall survival (%)", side=2, padj=-3, cex=1.5)
	text(ci.lab.xycoord[1], ci.lab.xycoord[2], label=texts, font=1, cex=1.2)
	#legend(legend.xycoord[1], legend.xycoord[2], legend=c(paste("mutant (events ", mutant_count, ")", sep=""), paste("wild-type (events ", wild_type_count, ")", sep="")), fill=c("brown", "cyan4"), box.lwd=0, box.col='white', cex=1.2, border='white')
	#legend(legend.xycoord[1], legend.xycoord[2], legend=c(paste("Hypermutated (events ", mutant_count, ")", sep=""), paste("Non-hypermutated (events ", wild_type_count, ")", sep="")), fill=c("brown", "cyan4"), box.lwd=0, box.col='white', cex=1.2, border='white',x.intersp=0.2)

	legend.label1=paste(legend.labels[1], " (events ", mutant_count, ")", sep="")
	legend.label2=paste(legend.labels[2], " (events ", wild_type_count, ")", sep="")
	legend(legend.xycoord[1], legend.xycoord[2], legend=c(legend.label1, legend.label2), fill=c("brown", "cyan4"), box.lwd=0, box.col='white', cex=1.2, border='white',x.intersp=0.2)
	if (!is.null(pngfig)) {
		dev.off()
	}
	invisible(list(surv.diff=diff, coxph.res=c))
}
