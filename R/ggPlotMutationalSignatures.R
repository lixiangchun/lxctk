
theme_ss <- function() {
  
  t = theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(hjust = 0.5))
  
  return(t)
}
theme_small_axis <- function(x = TRUE, y = TRUE, size = 6, family = "mono") {
  ## decrease the x/y-axis label size
  template = element_text(size = size, family = family)
  t = theme_ss()
  if(x)
    t = t + theme(axis.text.x = template)
  if(y)
    t = t + theme(axis.text.y = template)
  
  return(t)
}  
#setwd('/Users/lixiangchun/Public/WorkSpace/Project/TCGA/Analysis/decipherMutationalProcesses/als/output')
#infile = 'Rank_eq_13.processes.txt'
infile="/Users/lixiangchun/Public/WorkSpace/Project/TJMUCH-Prof.ChenKexin-CRC/lixc/revision/decipherMutationalSignatures/output/Rank_eq_3.processes.txt"

## infile: mutational process file from decipherMutationalSignatures
## d: a data frame refers to infile
## color_pallete: ColorBrewer color name
## 
ggPlotMutationalSignatures <- function(mutational.processes.file, d=NULL, color_pallete="Dark2", xlabel="Mutational contexts", ylabel="Percentage of mutations", y.upper.value=NA)
{
  if (is.null(d)) {
    d <- read.table(mutational.processes.file, header = TRUE, stringsAsFactors = FALSE)
  }
  d <- lxctk::sortDataFrame(d, c('types', 'subtypes'))
  d <- reshape2::melt(d)
  d$alteration <- d$types
  d$context <- d$subtypes
  substr(d$context, 2, 2) <- '.'
  d$motif <- paste(d$types, d$context)
  d <- d[, c('motif', 'variable', 'value', 'alteration', 'context')]
  colnames(d) <- c('motif', 'signature', 'value', 'alteration', 'context')
  d$value <- d$value * 100
  if (!is.na(y.upper.value)) {
	d$value[d$value > y.upper.value] = y.upper.value
  }
  d$signature <- gsub("processes", "Signature", d$signature)
  
  p = ggplot(d)
  p = p + geom_bar(
    aes_string(x = "context", y = "value", fill = "alteration"),
    stat = "identity",
    position = "identity"
  )
  p = p + facet_grid(signature ~ alteration)
  p = p + theme_ss() + theme_small_axis()
  p = p + theme(legend.position = "none")
  p = p + scale_fill_brewer(palette = color_pallete)
  p = p + xlab(xlabel) + ylab(ylabel)
  p
}
