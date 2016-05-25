

unicoxph <- function(cox.f,data,variable.display=NULL, plot=TRUE, ...)
{
	variables <- all.vars(cox.f)
	n <- length(variables)
	lci <- rep(NA, n)
	uci <- rep(NA, n)
	OR <- rep(NA, n)
	p.value <- rep(NA, n)
	
	ci.l <- list()	
	OR.l <- list()
	p.value.l <- list()

	for (i in 3:n) {
		e <- variables[i]
		#surv.f <- eval( do.call('substitute', list(surv.f.template, list(marker=as.name(e)))) )
		surv.f <- as.formula( sprintf("Surv(%s,%s) ~ %s", variables[1], variables[2], e) )
		r <- coxph(surv.f, data=data)
		ci <- exp(confint(r))
		odd.ratio <- exp(r$coefficients)
		p <- summary(r)$coefficients[,5]
		ci.l[[e]] <- ci
		OR.l[[e]] <- odd.ratio
		p.value.l[[e]] <- p
	}
	
	ci.x <- do.call('rbind', ci.l)
	OR.x <- do.call('c', OR.l)
	p.value.x <- do.call('c', p.value.l)
	p.value.x <- sprintf("%.3g", p.value.x)
	lci <- ci.x[,1]
	uci <- ci.x[,2]
	CI <- sprintf('%.3g (%.3g ~ %.3g)', OR.x, lci,uci)

	if (is.null(variable.display)) variable.display <- names(OR.x)
	coxTable <- cbind(variable.display, CI, p.value.x, OR.x, lci, uci)
	coxTable <- rbind(c('Variable','HR (95% CI)','P-value', 'OR', 'lci', 'uci'), coxTable)
	if (plot)
		plot.coxph(coxTable=coxTable, ...)
	invisible(coxTable)
}

