

COSMIC_signatures_annotation <- function(x, cutoff=0.7, annot_outfile)
{
	Signatures <- list(Signature.1 = "Cancer types: Signature 1 has been found in all cancer types and in most cancer samples.\nProposed aetiology: Signature 1 is the result of an endogenous mutational process initiated by spontaneous deamination of 5-methylcytosine.\nAdditional mutational features: Signature 1 is associated with small numbers of small insertions and deletions in most tissue types.\nComments: The number of Signature 1 mutations correlates with age of cancer diagnosis.",
	Signature.2 = "Cancer types: Signature 2 has been found in 22 cancer types, but most commonly in cervical and bladder cancers. In most of these 22 cancer types, Signature 2 is present in at least 10% of samples.\nProposed aetiology: Signature 2 has been attributed to activity of the AID/APOBEC family of cytidine deaminases. On the basis of similarities in the sequence context of cytosine mutations caused by APOBEC enzymes in experimental systems, a role for APOBEC1, APOBEC3A and/or APOBEC3B in human cancer appears more likely than for other members of the family.\nAdditional mutational features: Transcriptional strand bias of mutations has been observed in exons, but is not present or is weaker in introns.\nComments: Signature 2 is usually found in the same samples as Signature 13. It has been proposed that activation of AID/APOBEC cytidine deaminases is due to viral infection, retrotransposon jumping or to tissue inflammation. Currently, there is limited evidence to support these hypotheses. A germline deletion polymorphism involving APOBEC3A and APOBEC3B is associated with the presence of large numbers of Signature 2 and 13 mutations and with predisposition to breast cancer. Mutations of similar patterns to Signatures 2 and 13 are commonly found in the phenomenon of local hypermutation present in some cancers, known as kataegis, potentially implicating AID/APOBEC enzymes in this process as well.",
	Signature.3 = "Cancer types:Signature 3 has been found in breast, ovarian, and pancreatic cancers.\nProposed aetiology: Signature 3 is associated with failure of DNA double-strand break-repair by homologous recombination.\nAdditional mutational features: Signature 3 associates strongly with elevated numbers of large (longer than 3bp) insertions and deletions with overlapping microhomology at breakpoint junctions.\nComments: Signature 3 is strongly associated with germline and somatic BRCA1 and BRCA2 mutations in breast, pancreatic, and ovarian cancers. In pancreatic cancer, responders to platinum therapy usually exhibit Signature 3 mutations.",
	Signature.4="Cancer types: Signature 4 has been found in head and neck cancer, liver cancer, lung adenocarcinoma, lung squamous carcinoma, small cell lung carcinoma, and oesophageal cancer.\nProposed aetiology: Signature 4 is associated with smoking and its profile is similar to the mutational pattern observed in experimental systems exposed to tobacco carcinogens (e.g., benzo[a]pyrene). Signature 4 is likely due to tobacco mutagens.\nAdditional mutational features: Signature 4 exhibits transcriptional strand bias for C>A mutations, compatible with the notion that damage to guanine is repaired by transcription-coupled nucleotide excision repair. Signature 4 is also associated with CC>AA dinucleotide substitutions.\nComments: Signature 29 is found in cancers associated with tobacco chewing and appears different from Signature 4.",
	Signature.5="Cancer types: Signature 5 has been found in all cancer types and most cancer samples.\nProposed aetiology: The aetiology of Signature 5 is unknown.\nAdditional mutational features: Signature 5 exhibits transcriptional strand bias for T>C substitutions at ApTpN context.\nComments: N/A",
	Signature.6="Cancer types: Signature 6 has been found in 17 cancer types and is most common in colorectal and uterine cancers. In most other cancer types, Signature 6 is found in less than 3% of examined samples.\nProposed aetiology: Signature 6 is associated with defective DNA mismatch repair and is found in microsatellite unstable tumours.\nAdditional mutational features: Signature 6 is associated with high numbers of small (shorter than 3bp) insertions and deletions at mono/polynucleotide repeats.\nComments: Signature 6 is one of four mutational signatures associated with defective DNA mismatch repair and is often found in the same samples as Signatures 15, 20, and 26.",
	Signature.7="Cancer types: Signature 7 has been found predominantly in skin cancers and in cancers of the lip categorized as head and neck or oral squamous cancers.\nProposed aetiology: Based on its prevalence in ultraviolet exposed areas and the similarity of the mutational pattern to that observed in experimental systems exposed to ultraviolet light Signature 7 is likely due to ultraviolet light exposure.\nAdditional mutational features: Signature 7 is associated with large numbers of CC>TT dinucleotide mutations at dipyrimidines. Additionally, Signature 7 exhibits a strong transcriptional strand-bias indicating that mutations occur at pyrimidines (viz., by formation of pyrimidine-pyrimidine photodimers) and these mutations are being repaired by transcription-coupled nucleotide excision repair.\nComments: N/A",
	Signature.8="Cancer types: Signature 8 has been found in breast cancer and medulloblastoma.\nProposed aetiology: The aetiology of Signature 8 remains unknown.\nAdditional mutational features: Signature 8 exhibits weak strand bias for C>A substitutions and is associated with double nucleotide substitutions, notably CC>AA.\nComments:N/A",
	Signature.9="Cancer types: Signature 9 has been found in chronic lymphocytic leukaemias and malignant B-cell lymphomas.\nProposed aetiology: Signature 9 is characterized by a pattern of mutations that has been attributed to polymerase η, which is implicated with the activity of AID during somatic hypermutation.\nAdditional mutational features: N/A\nComments: Chronic lymphocytic leukaemias that possess immunoglobulin gene hypermutation (IGHV-mutated) have elevated numbers of mutations attributed to Signature 9 compared to those that do not have immunoglobulin gene hypermutation.",
	Signature.10="Cancer types: Signature 10 has been found in six cancer types, notably colorectal and uterine cancer, usually generating huge numbers of mutations in small subsets of samples.\nProposed aetiology: It has been proposed that the mutational process underlying this signature is altered activity of the error-prone polymerase POLE. The presence of large numbers of Signature 10 mutations is associated with recurrent POLE somatic mutations, viz., Pro286Arg and Val411Leu.\nAdditional mutational features: Signature 10 exhibits strand bias for C>A mutations at TpCpT context and T>G mutations at TpTpT context.\nComments: Signature 10 is associated with some of most mutated cancer samples. Samples exhibiting this mutational signature have been termed ultra-hypermutators.",
	Signature.11="Cancer types: Signature 11 has been found in melanoma and glioblastoma.\nProposed aetiology: Signature 11 exhibits a mutational pattern resembling that of alkylating agents. Patient histories have revealed an association between treatments with the alkylating agent temozolomide and Signature 11 mutations.\nAdditional mutational features: Signature 11 exhibits a strong transcriptional strand-bias for C>T substitutions indicating that mutations occur on guanine and that these mutations are effectively repaired by transcription-coupled nucleotide excision repair.\nComments: N/A",
	Signature.12="Cancer types: Signature 12 has been found in liver cancer.\nProposed aetiology: The aetiology of Signature 12 remains unknown.\nAdditional mutational features: Signature 12 exhibits a strong transcriptional strand-bias for T>C substitutions.\nComments: Signature 12 usually contributes a small percentage (<20%) of the mutations observed in a liver cancer sample.",
	Signature.13="Cancer types: Signature 13 has been found in 22 cancer types and seems to be commonest in cervical and bladder cancers. In most of these 22 cancer types, Signature 13 is present in at least 10% of samples.\nProposed aetiology: Signature 13 has been attributed to activity of the AID/APOBEC family of cytidine deaminases converting cytosine to uracil. On the basis of similarities in the sequence context of cytosine mutations caused by APOBEC enzymes in experimental systems, a role for APOBEC1, APOBEC3A and/or APOBEC3B in human cancer appears more likely than for other members of the family. Signature 13 causes predominantly C>G mutations. This may be due to generation of abasic sites after removal of uracil by base excision repair and replication over these abasic sites by REV1.\nAdditional mutational features: Transcriptional strand bias of mutations has been observed in exons, but is not present or is weaker in introns.\nComments: Signature 2 is usually found in the same samples as Signature 13. It has been proposed that activation of AID/APOBEC cytidine deaminases is due to viral infection, retrotransposon jumping or to tissue inflammation. Currently, there is limited evidence to support these hypotheses. A germline deletion polymorphism involving APOBEC3A and APOBEC3B is associated with the presence of large numbers of Signature 2 and 13 mutations and with predisposition to breast cancer. Mutations of similar patterns to Signatures 2 and 13 are commonly found in the phenomenon of local hypermutation present in some cancers, known as kataegis, potentially implicating AID/APOBEC enzymes in this process as well.",
	Signature.14="Cancer types: Signature 14 has been observed in four uterine cancers and a single adult low-grade glioma sample.\nProposed aetiology: The aetiology of Signature 14 remains unknown.\nAdditional mutational features: N/A\nComments: Signature 14 generates very high numbers of somatic mutations (>200 mutations per MB) in all samples in which it has been observed.",
	Signature.15="Cancer types: Signature 15 has been found in several stomach cancers and a single small cell lung carcinoma.\nProposed aetiology: Signature 15 is associated with defective DNA mismatch repair.\nAdditional mutational features: Signature 15 is associated with high numbers of small (shorter than 3bp) insertions and deletions at mono/polynucleotide repeats.\nComments: Signature 15 is one of four mutational signatures associated with defective DNA mismatch repair and is often found in the same samples as Signatures 6, 20, and 26.",
	Signature.16="Cancer types: Signature 16 has been found in liver cancer.\nProposed aetiology: The aetiology of Signature 16 remains unknown.\nAdditional mutational features: Signature 16 exhibits an extremely strong transcriptional strand bias for T>C mutations at ApTpN context, with T>C mutations occurring almost exclusively on the transcribed strand.\nComments: N/A",
	Signature.17="Cancer types: Signature 17 has been found in oesophagus cancer, breast cancer, liver cancer, lung adenocarcinoma, B-cell lymphoma, stomach cancer and melanoma.\nProposed aetiology: The aetiology of Signature 17 remains unknown.\nAdditional mutational features: N/A\nComments: N/A",
	Signature.18="Cancer types: Signature 18 has been found commonly in neuroblastoma. Additionally, Signature 18 has been also observed in breast and stomach carcinomas.\nProposed aetiology: The aetiology of Signature 18 remains unknown.\nAdditional mutational features: N/A\nComments:N/A",
	Signature.19="Cancer types: Signature 19 has been found only in pilocytic astrocytoma.\nProposed aetiology: The aetiology of Signature 19 remains unknown.\nAdditional mutational features: N/A\nComments: N/A",
	Signature.20="Cancer types: Signature 20 has been found in stomach and breast cancers.\nProposed aetiology: Signature 20 is believed to be associated with defective DNA mismatch repair.\nAdditional mutational features: Signature 20 is associated with high numbers of small (shorter than 3bp) insertions and deletions at mono/polynucleotide repeats.\nComments: Signature 20 is one of four mutational signatures associated with defective DNA mismatch repair and is often found in the same samples as Signatures 6, 15, and 26.",
	Signature.21="Cancer types: Signature 21 has been found only in stomach cancer.\nProposed aetiology: The aetiology of Signature 21 remains unknown.\nAdditional mutational features: N/A\nComments: Signature 21 is found only in four samples all generated by the same sequencing centre. The mutational pattern of Signature 21 is somewhat similar to the one of Signature 26. Additionally, Signature 21 is found only in samples that also have Signatures 15 and 20. As such, Signature 21 is probably also related to microsatellite unstable tumours.",
	Signature.22="Cancer types: Signature 22 has been found in urothelial (renal pelvis) carcinoma and liver cancers.\nProposed aetiology: Signature 22 has been found in cancer samples with known exposures to aristolochic acid. Additionally, the pattern of mutations exhibited by the signature is consistent with the one previous observed in experimental systems exposed to aristolochic acid.\nAdditional mutational features: Signature 22 exhibits a very strong transcriptional strand bias for T>A mutations indicating adenine damage that is being repaired by transcription-coupled nucleotide excision repair.\nComments: Signature 22 has a very high mutational burden in urothelial carcinoma; however, its mutational burden is much lower in liver cancers.",
	Signature.23="Cancer types: Signature 23 has been found only in a single liver cancer sample.\nProposed aetiology: The aetiology of Signature 23 remains unknown.\nAdditional mutational features: Signature 23 exhibits very strong transcriptional strand bias for C>T mutations.\nComments:N/A",
	Signature.24="Cancer types: Signature 24 has been observed in a subset of liver cancers.\nProposed aetiology: Signature 24 has been found in cancer samples with known exposures to aflatoxin. Additionally, the pattern of mutations exhibited by the signature is consistent with that previous observed in experimental systems exposed to aflatoxin.\nAdditional mutational features: Signature 24 exhibits a very strong transcriptional strand bias for C>A mutations indicating guanine damage that is being repaired by transcription-coupled nucleotide excision repair.\nComments: N/A",
	Signature.25="Cancer types: Signature 25 has been observed in Hodgkin lymphomas.\nProposed aetiology: The aetiology of Signature 25 remains unknown.\nAdditional mutational features: Signature 25 exhibits transcriptional strand bias for T>A mutations.\nComments: This signature has only been identified in Hodgkin’s cell lines. Data is not available from primary Hodgkin lymphomas.",
	Signature.26="Cancer types: Signature 26 has been found in breast cancer, cervical cancer, stomach cancer and uterine carcinoma.\nProposed aetiology: Signature 26 is believed to be associated with defective DNA mismatch repair.\nAdditional mutational features: Signature 26 is associated with high numbers of small (shorter than 3bp) insertions and deletions at mono/polynucleotide repeats.\nComments: Signature 26 is one of four mutational signatures associated with defective DNA mismatch repair and is often found in the same samples as Signatures 6, 15 and 20.",
	Signature.27="Cancer types: Signature 27 has been observed in a subset of kidney clear cell carcinomas.\nProposed aetiology: The aetiology of Signature 27 remains unknown.\nAdditional mutational features: Signature 27 exhibits very strong transcriptional strand bias for T>A mutations. Signature 27 is associated with high numbers of small (shorter than 3bp) insertions and deletions at mono/polynucleotide repeats.\nComments: N/A",
	Signature.28="Cancer types: Signature 28 has been observed in a subset of stomach cancers.\nProposed aetiology: The aetiology of Signature 28 remains unknown.\nAdditional mutational features: N/A\nComments: N/A",
	Signature.29="Cancer types: Signature 29 has been observed only in gingivo-buccal oral squamous cell carcinoma.\nProposed aetiology: Signature 29 has been found in cancer samples from individuals with a tobacco chewing habit.\nAdditional mutational features: Signature 29 exhibits transcriptional strand bias for C>A mutations indicating guanine damage that is most likely repaired by transcription-coupled nucleotide excision repair. Signature 29 is also associated with CC>AA dinucleotide substitutions.\nComments: The Signature 29 pattern of C>A mutations due to tobacco chewing appears different from the pattern of mutations due to tobacco smoking reflected by Signature 4.",
	Signature.30="Cancer types: Signature 30 has been observed in a small subset of breast cancers.\nProposed aetiology: The aetiology of Signature 30 remains unknown.\nAdditional mutational features: N/A\nComments: N/A")
	f = file(annot_outfile, 'w')
	nc = ncol(x)
	Signature_Names = colnames(x)
	My_Signature_Names = rownames(x)
	for (i in 1:nc) {
		y <- as.numeric(x[,i])  ## cosine similarities
		if (any(y > cutoff)) {
			annot <- try(Signatures[[Signature_Names[i]]], silent=TRUE)
			if (class(annot) == 'try-error') { 
				annot = 'null'
			}
			Signature_Name = paste(c( "Your input:",My_Signature_Names[which(y>cutoff)],"\nCOSMIC signature:",Signature_Names[[i]],  "\nCosine similarity:", paste(signif(y[y>=cutoff]*100,4),"%",sep="")  ))
			cat(Signature_Name, file=f)
			cat("\n",file=f)
			cat(annot,file=f)
			cat("\n",file=f)
			cat("\n",file=f)
		}
	}	
	close(f)
}

## source code from "lsa" package
cosine <- function (x, y = NULL) 
{
    if (is.matrix(x) && is.null(y)) {
        co = array(0, c(ncol(x), ncol(x)))
        f = colnames(x)
        dimnames(co) = list(f, f)
        for (i in 2:ncol(x)) {
            for (j in 1:(i - 1)) {
                co[i, j] = cosine(x[, i], x[, j])
            }
        }
        co = co + t(co)
        diag(co) = 1
        return(as.matrix(co))
    }
    else if (is.vector(x) && is.vector(y)) {
        return(crossprod(x, y)/sqrt(crossprod(x) * crossprod(y)))
    }
    else if (is.vector(x) && is.matrix(y)) {
        co = vector(mode = "numeric", length = ncol(y))
        names(co) = colnames(y)
        for (i in 1:ncol(y)) {
            co[i] = cosine(x, y[, i])
        }
        return(co)
    }
    else {
        stop("argument mismatch. Either one matrix or two vectors needed as input.")
    }
}

COSMIC_signatures_corr <- function(mutational_processes_file, plot=FALSE, outfile, cutoff=0.7, annot_outfile, ...)
{
	d <- read.table(mutational_processes_file, header=TRUE, stringsAsFactors=FALSE)	
	data('COSMIC_signatures')
	Alexandrov <- COSMIC_signatures[,4:ncol(COSMIC_signatures)]
	lixc <- d[,3:ncol(d)]

	x <- matrix(nr=ncol(lixc), nc=ncol(Alexandrov), NA) 
	for (i in 1:ncol(lixc)) {
	 for (j in 1:ncol(Alexandrov)) {
	     x[i, j] <- cosine(lixc[,i], Alexandrov[,j])
	 }   
	}
	colnames(x) <- colnames(Alexandrov)
	rownames(x) <- paste('MSig',1:nrow(x),sep="")
	if (!missing(outfile)) {
		write.table(x, file=outfile, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
	}
	if (plot) {
		txt <- round(x,2)
		txt[txt<cutoff] <- NA
		NMF::aheatmap(x, 'Reds', txt=txt, ...)
	}
	# plot with plotly
	#Cosine <- x
	#p = plot_ly(z=Cosine,type='heatmap', colorscale='Reds', y=rownames(x), x=colnames(x)) %>% layout(yaxis=list(title='Li X.C. WGS et al'),xaxis=list(title='Alexandrov et al. 21 mutational signatures',ticklen=0, showticklabels=FALSE))
	#htmlwidgets::saveWidget(as.widget(p), "Cosine_similarity_with_Alexandrov_etal.html")

	if (!missing(annot_outfile)) {
		COSMIC_signatures_annotation(x, cutoff, annot_outfile)
	}
	invisible(x)
}
