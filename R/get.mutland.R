get.mutland <- function(genes=NULL, mutsigcv.sig.gene.file, gene.mut.freq.file, mut.per.smp.file, gene.mut.mat.file, use.q.value=TRUE, discard.genes=c("TTN"))
{
	dGeneMutFreq=read.table(gene.mut.freq.file, stringsAsFactors=FALSE)	
	dMutPerSmp=read.table(mut.per.smp.file, stringsAsFactors=FALSE)
	dGeneMutMat=read.table(gene.mut.mat.file, stringsAsFactors=FALSE)
	dMutSig=read.table(mutsigcv.sig.gene.file, header=TRUE, stringsAsFactors=FALSE, sep="\t")

	if (is.null(genes)) {
		q.vals=c(0.1, 0.2, 0.25)
		genes_lt_q0.1=dMutSig$gene[dMutSig$q<q.vals[1]]
		genes_lt_q0.2=dMutSig$gene[dMutSig$q<q.vals[2]]
		genes_lt_q0.25=dMutSig$gene[dMutSig$q<q.vals[3]]
		genelist=list(genes_lt_q0.1, genes_lt_q0.2, genes_lt_q0.25)	
		x=c(length(genes_lt_q0.1), length(genes_lt_q0.2), length(genes_lt_q0.25)) - 10
		i=which.min(x)
		genes=genelist[[i]]
		cat("Automaticall selected q-value cutoff = ", q.vals[i], ".\n", sep="")
		if (!is.null(discard.genes)) {
			genes=setdiff(genes, discard.genes)
		}
	}
	
	if (length(genes) == 0) {
		genes=dMutSig$gene[1:12]
		cat("CAUTION: Number of MutSig genes is zero, use top 10 genes!\n", file=stderr())
	}
	genes=intersect(colnames(dGeneMutMat), genes)

	upper.panel=data.frame(sample=rownames(dMutPerSmp), silence=dMutPerSmp$Syn, nonsilence=dMutPerSmp$TotalNum-dMutPerSmp$Syn)
	middle.panel=dGeneMutMat[, genes]

	gene.freq=dGeneMutFreq$NonsilentFreq
	names(gene.freq)=rownames(dGeneMutFreq)	
	nonsilence_sample_num=dGeneMutFreq$MutatedSampleNum
	names(nonsilence_sample_num)=rownames(dGeneMutFreq)

#	syn=dGeneMutFreq$Syn
#	names(syn)=rownames(dGeneMutFreq)

	nonsilent_sample_percent=dGeneMutFreq$NonsilentFreq * 100
	nonsilent_sample_percent=paste(format(nonsilent_sample_percent, digits=2), "%", sep="")
	names(nonsilent_sample_percent)=rownames(dGeneMutFreq)
	left.panel=data.frame(gene=genes,nonsilence_sample_num=nonsilence_sample_num[genes], nonsilent_sample_percent=nonsilent_sample_percent[genes])	

#	left.panel=data.frame(gene=genes, syn=syn[genes], nonsilence_sample_num=nonsilence_sample_num[genes], nonsilent_sample_percent=nonsilent_sample_percent[genes])	

	if (use.q.value)
	{
		q=dMutSig$q
		names(q)=dMutSig$gene
	}
	if (!use.q.value)
	{
		q=dMutSig$p
		names(q)=dMutSig$gene
	}
	q[q==0] <- min(q[q>0]) / 10
	right.panel=try(data.frame(gene=genes, q=q[genes]))
	if (class(right.panel) == 'try-error') {right.panel=NULL}

	list(upper.panel=upper.panel, left.panel=left.panel, middle.panel=t(middle.panel), right.panel=right.panel)
}

#r=landscape.data2(NULL,'mutsig.out.mutsig.gene.txt', '../iCGA/BasicAna/GastricCancer.gene_mut_freq.txt','../iCGA/BasicAna/GastricCancer.mut_per_smp.txt','../iCGA/BasicAna/GastricCancer.gene_mut_mat.txt')

