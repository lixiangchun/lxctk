

run.cn.farms.from.Celfiles <- function(filenames, outfile="segmentation.RData", annotDir=NULL, 
						normMethod="SOR", cores=2, runtype="ff") {

createAnnotation(filenames = filenames, annotDir=annotDir)
print("Finished running createAnnotation(...)", file=stderr())

normData <- normalizeCels(filenames, method=normMethod, cores=cores, alleles=TRUE, runtype=runtype, annotDir=annotDir)
print("Finished running normalizeCels(...)", file=stderr())


summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(10)
callParam <- list(cores=cores, runtype=runtype)
slData <- slSummarization(normData,
	summaryMethod=summaryMethod,
	summaryParam=summaryParam,
	callParam=callParam,
	summaryWindow="std")

print("Finished running slSummarization(...)", file=stderr())
show(slData)

npData <- normalizeNpData(filenames, cores, runtype=runtype, annotDir=annotDir)
print("Finished running normalizeNpData(...)", file=stderr())

combData <- combineData(slData, npData, runtype=runtype)
print("Finished running combineData(...)", file=stderr())

show(combData)


windowMethod <- "std"
windowParam <- list()
windowParam$windowSize <- 5
windowParam$overlap <- TRUE
summaryMethod <- "Variational"
summaryParam <- list()
summaryParam$cyc <- c(20)
callParam <- list(cores=cores, runtype=runtype)
mlData <- mlSummarization(slData,
                 windowMethod =windowMethod,
                 windowParam  =windowParam,
                 summaryMethod=summaryMethod,
                 summaryParam =summaryParam,
                 callParam    =callParam)

print("Finished running mlSummarization(...)", file=stderr())

colnames(assayData(mlData)$L_z) <- sampleNames(mlData)
segments <- dnaCopySf(
	x = assayData(mlData)$L_z,
	chrom = fData(mlData)$chrom,
	maploc = fData(mlData)$start,
	cores = cores,
	smoothing=FALSE
)

print("Finished running dnaCopySf(...)", file=stderr())
save(segments, file=outfile)
print("Successfully finished running run.cn.farms.from.Celfiles(...)", file=stderr())

}


##filenames=list.files("/ifshk5/PC_HUMAN_AP/PMO/F13TSHNCKF0797_HUMdjpX/lixc/STAD_SNP6.0/GSE31168", "*.CEL", full.names=TRUE)
##run.cn.farms.from.Celfiles(filenames, outfile="segmentation.RData", annotDir=NULL, normMethod="SOR", cores=8, runtype="ff")

