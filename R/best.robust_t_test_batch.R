###!/ifshk4/BC_PUB/biosoft/PIPE_RD/Package/R-3.1.1/bin/Rscript

normality_test <- function(vals1, vals2)
{
	r1 <- try(shapiro.test(vals1), silent=TRUE)
	r2 <- try(shapiro.test(vals2), silent=TRUE)
	r3 <- try(shapiro.test(c(vals1, vals2)), silent=TRUE)
	if (class(r1) == 'try-error') {
		p1 <- NA
	} else {
		p1 <- r1$p.value
	}
	if (class(r2) == 'try-error') {
		p2 <- NA
	} else {
		p2 <- r2$p.value
	}
	if (class(r3) == 'try-error') {
		p3 <- NA
	} else {
		p3 <- r3$p.value
	}
	p <- c(p1, p2, p3)
	p
}

best.robust_t_test_batch <- 
	function(tumor_expression, normal_expression, best.robust_t_test_model_stanDso, outfile, ...)
{

tumor_expression_rownames <- rownames(tumor_expression)
normal_expression_rownames <- rownames(normal_expression)
if (any(tumor_expression_rownames != normal_expression_rownames)) {
	print('The column names of 2 input data frame does not match!', file=stderr())
	quit('no')
}

#if (!is.null(infile)) {
#	genes <- read.table(infile,header=FALSE,stringsAsFactors=FALSE)$V1
#} else {
#	genes <- rownames(normal_expression)
#}

genes <- normal_expression_rownames
print(length(genes))

normality_test_pvals <- list()
mu1 <- list()
mu2 <- list()
sd1 <- list()
sd2 <- list()
nu1 <- list()
nu2 <- list()
mu_diff <- list()
sigma_diff <- list()
nu_diff <- list()

num <- list()
N <- length(genes)
t.test.pvals <- rep(1, N)
t.test.grp <- list()
wilcox.test.grp <- list()
for (i in 1:N)
{
	gene <- genes[i]

	vals1 <- as.numeric(normal_expression[gene,]) ## normal gene expression
	vals2 <- as.numeric(tumor_expression[gene,])  ## tumor gene expression

	p <- normality_test(vals1, vals2)
	normality_test_pvals[[gene]] <- p

	vals <- c(vals2, vals1)
	grp <- c(rep(1, length(vals2)), rep(2, length(vals1)))
	r <- try(t.test(vals2, vals1))
	if (class(r) == 'try-error') {
		t.test.grp[[gene]] <- c(NA, NA, NA)
		t.pval <- 1
	} else {
		t.test.grp[[gene]] <- c(r$estimate, r$p.value)
		t.pval <- r$p.value
	}
	wt <- try(wilcox.test(vals2, vals1))
	if (class(wt) == 'try-error') {
		wilcox.test.grp[[gene]] <- c(NA, NA, NA)
		wt.pval <- 1
	} else {
		m1 <- median(vals2)
		m2 <- median(vals1) 
		wilcox.test.grp[[gene]] <- c(m1, m2, wt$p.value)
		wt.pval <- wt$p.value
	}
#if (t.pval < 0.05 || wt.pval < 0.05) {
	fit <- best.robust_t_test(vals,grp,best.robust_t_test_model_stanDso,parameters=c('mu','sigma','nu','mu_diff','sigma_diff', 'nu_diff'),...)
	if (!is.array(fit)) {
		## MCMC often get stuck if the region [unifLo, unifHi] is too large. Just use this if the above failed!
		fit <- best.robust_t_test(vals, grp, best.robust_t_test_model_stanDso, unifLo=0,unifHi=.Machine$double.xmax,parameters = c('mu','sigma','nu','mu_diff','sigma_diff', 'nu_diff'), ...)
	}
#} else {
#	fit <- NA
#}
	if (!is.array(fit)) {
		print('is.array(fit) is FALSE.')
		NA_Vals <- c(NA, NA, NA, NA, NA)
		mu1[[gene]] <- NA_Vals
		mu2[[gene]] <- NA_Vals
		sd1[[gene]] <- NA_Vals
		sd2[[gene]] <- NA_Vals
		nu1[[gene]] <- NA_Vals
		nu2[[gene]] <- NA_Vals
		mu_diff[[gene]] <- NA_Vals
		sigma_diff[[gene]] <- NA_Vals
		nu_diff[[gene]] <- NA_Vals
	} else {	
		##if (class(fit) == 'try-error') next;
		x <- summary(fit)$summary
		##mu1[[gene]] <- x[,c('mean','2.5%','97.5%','Rhat')][1,]
		mu1[[gene]] <- x[,c('mean','2.5%','97.5%','n_eff','Rhat')][1,]
		mu2[[gene]] <- x[,c('mean','2.5%','97.5%','n_eff','Rhat')][2,]
		sd1[[gene]] <- x[,c('mean','2.5%','97.5%','n_eff','Rhat')][3,]
		sd2[[gene]] <- x[,c('mean','2.5%','97.5%','n_eff','Rhat')][4,]
		nu1[[gene]] <- x[,c('mean','2.5%','97.5%','n_eff','Rhat')][5,]
		nu2[[gene]] <- x[,c('mean','2.5%','97.5%','n_eff','Rhat')][6,]
		mu_diff[[gene]] <- x[,c('mean','2.5%','97.5%','n_eff','Rhat')][7,]
		sigma_diff[[gene]] <- x[,c('mean','2.5%','97.5%','n_eff','Rhat')][8,]
		nu_diff[[gene]] <- x[,c('mean','2.5%','97.5%','n_eff','Rhat')][9,]
		print(x[,c('mean','2.5%','97.5%','n_eff','Rhat')])
	}
	print(r)
	print(wt)
	cat("##GENE:", gene, "###idx:", i, "##Amp. Tumor:",length(vals2), ", ##Normal:", length(vals1), "##\n")
	num[[gene]] <- c(length(vals2), length(vals1))
	#print(showConnections(all=TRUE))
	##t.test.pvals[i] <- r$p.value
	### Close by hand, or there will be ERROR:
	### all connections are in use Stan model  does not contain samples.
	closeAllConnections()
}

x1 <- do.call('rbind', mu1)
x2 <- do.call('rbind', mu2)
y1 <- do.call('rbind', sd1)
y2 <- do.call('rbind', sd2)

normality_test_pvals <- do.call('rbind', normality_test_pvals)
colnames(normality_test_pvals) <- c('tumor_normality_pval', 'normal_normality_pval', 'combined_normality_pval')

mu_diff <- do.call('rbind', mu_diff)
sigma_diff <- do.call('rbind', sigma_diff)
nu_diff <- do.call('rbind', nu_diff)

normality1 <- do.call('rbind', nu1)
normality2 <- do.call('rbind', nu2)
number <- do.call('rbind', num)

t.test.grp <- do.call('rbind',t.test.grp)
wilcox.test.grp <- do.call('rbind', wilcox.test.grp)

colnames(x1) <- c('posterior.mean.tumor','posterior.mean.tumor.CI.2.5%','posterior.mean.tumor.CI.97.5%','mean.tumor.n_eff','mean.tumor.Rhat')
colnames(x2) <- c('posterior.mean.normal','posterior.mean.normal.CI.2.5%','posterior.mean.normal.CI.97.5%','mean.normal.n_eff','mean.normal.Rhat')
colnames(y1) <- c('posterior.sd.tumor','posterior.sd.tumor.CI.2.5%','posterior.sd.tumor.CI.97.5%','sd.tumor.n_eff','sd.tumor.Rhat')
colnames(y2) <- c('posterior.sd.normal','posterior.sd.normal.CI.2.5%','posterior.sd.normal.CI.97.5%','sd.normal.n_eff','sd.normal.Rhat')

colnames(t.test.grp) <- c('sample.mean.tumor', 'sample.mean.normal', 't.test.pval')
colnames(wilcox.test.grp) <- c('sample.median.tumor', 'sample.median.normal', 'wilcox.test.pval')

colnames(normality1) <- c('posterior.nu.tumor','posterior.nu.tumor.CI.2.5%','posterior.nu.tumor.CI.97.5%','nu.tumor.n_eff','nu.tumor.Rhat')
colnames(normality2) <- c('posterior.nu.normal','posterior.nu.normal.CI.2.5%','posterior.nu.normal.CI.97.5%','nu.normal.n_eff','nu.normal.Rhat')

colnames(mu_diff) <- c('posterior.mu_diff','posterior.mu_diff.CI.2.5%','posterior.mu_diff.CI.97.5%','mu_diff.n_eff','mu_diff.Rhat')
colnames(sigma_diff) <- c('posterior.sigma_diff','posterior.sigma_diff.CI.2.5%','posterior.sigma_diff.CI.97.5%','sigma_diff.n_eff','sigma_diff.Rhat')
colnames(nu_diff) <- c('posterior.nu_diff','posterior.nu_diff.CI.2.5%','posterior.nu_diff.CI.97.5%','nu_diff.n_eff','nu_diff.Rhat')


colnames(number) <- c('tumor.sample.num', 'normal.sample.num')

dat <- cbind(normality_test_pvals, t.test.grp, wilcox.test.grp, x1, y1, x2, y2, normality1, normality2, mu_diff, sigma_diff, nu_diff, number)
##write.table(dat, file='differentially_expressed_gene.bayes.txt',quote=FALSE,sep="\t")
write.table(dat, file=outfile,quote=FALSE,sep="\t")

print("Still_waters_run_deep")

}

#args <- commandArgs(TRUE)
#infile <- args[1]
#outfile <- args[2]

#print(args)
#normalizePath(path=c('/ifshk1/BC_CANCER/03user/lixiangchun/Software/INSTALL/Rpackages','/home/lixiangchun/R/x86_64-unknown-linux-gnu-library/3.1'))

#library(lxctk)
#library(rstan)
#library(Rcpp)
#library(methods)
#library(utils)
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

##load('../TCGA_Expression/LIHC_Expression_Matrix.RData')
# [1] "normal_expression" "tumor_expression" 
##load('../CNV/CNV_vs_Expression/best.robust_t_test_model_stanDso.RData')

#load('../../pre_processed_data/gdac.broadinstitute.org_COAD.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0.normal.txt.RData')
#print('finish loading normal data file.', file=stderr())
#load('../../pre_processed_data/gdac.broadinstitute.org_COAD.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0.tumor.txt.RData')
#print('finish loading tumor data file.', file=stderr())

#load('/ifshk5/PC_HUMAN_AP/PMO/F13TSHNCKF0797_HUMdjpX/lixc/TCGA/Colon_Cancer/Methylation/differential_probe_methyl_analysis/tumor_vs_normal/best.robust_t_test_model2_stanDso.RData')


#target_genes <- read.table(infile,header=FALSE,stringsAsFactors=FALSE)$V1
#tumor <- tumor[target_genes,]
#normal <- normal[target_genes,]
#best.robust_t_test_batch_running(tumor, normal, best.robust_t_test_model_stanDso, outfile=outfile, iter=4000)

