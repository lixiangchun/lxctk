## Draws a single "column" or "stack".
## X and Y coordinates determine the area of the column
## The Z coordinate determines the height of the column
## We include "lit=FALSE" arguments to remove the nasty shiny surfaces caused by lighting
stackplot.3d<-function(x,y,z,alpha=1,topcol="#078E53",sidecol="#aaaaaa"){
  
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to draw the sides and ends of the column separately  
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
  ## Determine the coordinates of each surface of the column and its edges
  x1=c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
  z1=c(rep(0,4),rep(c(0,0,z,z),4))
  y1=c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
  x2=c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2))
  z2=c(rep(c(0,z),4),rep(0,8),rep(z,8) )
  y2=c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2) )
  
  ## These lines create the sides of the column and its coloured top surface
  rgl.quads(x1,z1,y1,col=rep(sidecol,each=4),alpha=alpha,lit=FALSE)
  rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z,4),c(y[1],y[1],y[2],y[2]),
            col=rep(topcol,each=4),alpha=1,lit=FALSE) 
  ## This line adds black edges to the column
  rgl.lines(x2,z2,y2,col="#000000",lit=FALSE)
}
# Example:
# stackplot.3d(c(0,1),c(0,1),3,alpha=0.6)

## Calls stackplot.3d repeatedly to create a barplot
## z is the heights of the columns and must be an appropriately named vector
context3d<-function(z,alpha=1,scalexy=10,scalez=1,gap=0.2, 
		broadcolors=c("#805D3F","#72549A","#5EAFB2","#3F4F9D","#F2EC3C","#74B655")) {
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to add each column sequentially
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))

  ## Recreate Broad order
  types=c("C.G.G.C","T.A.A.T","C.A.G.T","T.G.A.C","C.T.G.A","T.C.A.G")
  contexts=c("TxT","CxT","AxT","GxT","TxC","CxC","AxC","GxC",
             "TxA","CxA","AxA","GxA","TxG","CxG","AxG","GxG")
  typeorder=c()
  for(type in types){
    typeorder=c(typeorder,paste(type,contexts,sep="_"))
  }
  z=z[typeorder]
  
  ## Reorder data into 6 regions
  set1=c(1:4,17:20,5:8,21:24,9:12,25:28,13:16,29:32)
  set2=set1+32
  set3=set1+64
  neworder=c(set1,set2,set3)

  ## Define dimensions of the plot 
  dimensions=c(12,8)
  
  ## Scale column area and the gap between columns 
  y=seq(1,dimensions[1])*scalexy
  x=seq(1,dimensions[2])*scalexy
  gap=gap*scalexy
  
  ## Scale z coordinate
  z=z*scalez
  
  ## Set up colour palette
  #broadcolors=c("#805D3F","#72549A","#5EAFB2","#3F4F9D","#F2EC3C","#74B655")
  colors=as.vector(sapply(broadcolors,rep,16))
  
  ## Plot each of the columns
  for(i in 1:dimensions[1]){
    for(j in 1:dimensions[2]){
      it=(i-1)*dimensions[2]+j # Variable to work out which column to plot; counts from 1:96
      stackplot.3d(c(gap+x[j],x[j]+scalexy),
                   c(-gap-y[i],-y[i]-scalexy),
                   z[neworder[it]],
                   alpha=alpha,
                   topcol=colors[neworder[it]],
                   sidecol=colors[neworder[it]])
    }
  }
  ## Set the viewpoint and add axes and labels
  rgl.viewpoint(theta=50,phi=40,fov=0)
  axes3d("y-+",labels=TRUE)
}

legoplot <- function(originalGenomes_file=system.file('data/stad_regular_GC_originalGenomes.txt', package='lxctk'), counts=NULL, alpha=0.8, scalez=0.5, gap=0.1, contexts_file=system.file('data/contexts.txt',package='lxctk'), snapshot_fig, piechart_fig)
{
	library(rgl)
	#contexts <- read.table('/Users/lixiangchun/Public/WorkSpace/Project/DigestiveSystemCancer/LiverCancer/decipherMutationalProcesses/mult/contexts.txt',header=FALSE,stringsAsFactors=FALSE)
	contexts <- read.table(contexts_file, header=FALSE, stringsAsFactors=FALSE)
	if (is.null(counts)) {
		x <- read.table(originalGenomes_file)
		counts <- rowSums(x)
	}
	names(counts) <- contexts$V3
	col1 <- c('#702E33','#BE5C81','#836E8D','#EC8C70','#596C22','#A7B762')
	C.A=col1[1]
	C.G=col1[2]
	C.T=col1[3]
	T.A=col1[4]
	T.C=col1[5]
	T.G=col1[6]

	col <- c(C.G,T.A,C.A,T.G,C.T,T.C)
	context3d(counts/32,alpha=alpha, scalez=scalez, gap=gap, broadcolors=col)
	#context3d(counts/32, alpha=0.8, scalez=0.5, gap=0.1, broadcolors=col)

	## Save your images to files if you wish
	if (!missing(snapshot_fig)) {
		rgl.snapshot(filename=snapshot_fig)
	}

	##broadcolors=c("#805D3F","#72549A","#5EAFB2","#3F4F9D","#F2EC3C","#74B655")

	col1 <- c(C.A,C.G,C.T,T.A,T.C,T.G)
	y <- by(counts, contexts$V2, sum)
	if (!missing(piechart_fig)) {
		pdf('pie.pdf',width=5, height=4)
		pie(y,col=col1)
		dev.off()
		print(y)
		print(y/sum(y))
	}
}
