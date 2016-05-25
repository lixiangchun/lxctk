

landscape.data <- function(genes=NULL, mutsigcv.sig.gene.file, gene.mut.freq.file, mut.per.smp.file, gene.mut.mat.file, use.q.value=T, top.k.gene=12, discard.genes=c("TTN"))
{
	if (is.null(genes)) {
		d=read.table(gene.mut.freq.file, stringsAsFactors=F)	
		genes=rownames(d)[1:top.k.gene]
		if (!is.null(discard.genes)) {
			genes=setdiff(genes, discard.genes)
		}
	}
	d=read.table(mut.per.smp.file, stringsAsFactors=F)
	upper.panel=data.frame(sample=rownames(d), silence=d$Syn, nonsilence=d$TotalNum-d$Syn)
	
	d=read.table(gene.mut.mat.file, stringsAsFactors=F)
	middle.panel=d[, genes]

	d=read.table(gene.mut.freq.file, stringsAsFactors=F)
	gene.freq=d$NonsilentFreq
	names(gene.freq)=rownames(d)	
	nonsilence_sample_num=d$MutatedSampleNum
	names(nonsilence_sample_num)=rownames(d)
#	syn=d$Syn
#	names(syn)=rownames(d)
	nonsilent_sample_percent=d$NonsilentFreq * 100
	nonsilent_sample_percent=paste(format(nonsilent_sample_percent, digits=2), "%", sep="")
	names(nonsilent_sample_percent)=rownames(d)
	left.panel=data.frame(gene=genes,nonsilence_sample_num=nonsilence_sample_num[genes], nonsilent_sample_percent=nonsilent_sample_percent[genes])	
#	left.panel=data.frame(gene=genes, syn=syn[genes], nonsilence_sample_num=nonsilence_sample_num[genes], nonsilent_sample_percent=nonsilent_sample_percent[genes])	

	d=read.table(mutsigcv.sig.gene.file, header=T, stringsAsFactors=F, sep="\t")
	if (use.q.value)
	{
		q=d$q
		names(q)=d$gene
	}
	if (!use.q.value)
	{
		q=d$p
		names(q)=d$gene
	}
	q[q==0] <- min(q[q>0]) / 10
	right.panel=try(data.frame(gene=genes, q=q[genes]))
	if (class(right.panel) == 'try-error') {right.panel=NULL}

	list(upper.panel=upper.panel, left.panel=left.panel, middle.panel=t(middle.panel), right.panel=right.panel)
}

#d=read.table('../basic_stats/WES.exonic.gene_mut_freq.txt', stringsAsFactors=F)
#genes=rownames(d)[d$MutatedSampleNum>=10]

#r=landscape.data(genes=NULL, "../MutSigCV/GC.WES.sig_genes.txt", "../basic_stats/WES.exonic.gene_mut_freq.txt", "../basic_stats/WES.exonic.mut_per_sample.txt", "../basic_stats/WES.exonic.gene_mut_mat.txt", F, top.k.gene=13, discard.genes="TTN")
#gc.subtype =read.table("subtype", stringsAsFactors=F)
#colnames(gc.subtype)=c('sample', 'subtype')
#gc.subtype$sample=paste(gc.subtype$sample, "-WES", sep="")
#gc.subtype$subtype[gc.subtype$subtype=='Diffuse']=1
#gc.subtype$subtype[gc.subtype$subtype=='Intestinal']=2
#gc.subtype$subtype=as.integer(gc.subtype$subtype)

#save(r, gc.subtype, file="r.RData")
