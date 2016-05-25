diag.mcmc <- function( codaObject, parName=varnames(codaObject)[1], 
                       DBDAplColors = c("skyblue","black","royalblue","steelblue"),
                       figName=NULL) {
  ##DBDAplColors = c("skyblue","black","royalblue","steelblue")
  #openGraph(height=5,width=7)
  if (class(codaObject) == 'stanfit' && class(codaObject) != 'mcmc.list') {
	codaObject <- stanfit2mcmc.list(codaObject)
  }
  if (!is.null(figName)) {
	png(figName, width=7*200, height=5*200, res=200)
  }
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot and gelman.plot are from CODA package:
  require(coda)
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors ) 
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
                       col=DBDAplColors )
  )  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning: coda::gelman.plot fails for ",parName))
    }
  }
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  #DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot2(codaObject,parName,plColors=DBDAplColors, main='')
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(figName) ) {
    #saveGraph( file=paste0(saveName,"Diag",parName), type=saveType)
	dev.off()	
  }
}
