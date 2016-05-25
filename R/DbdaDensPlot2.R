DbdaDensPlot2 <- function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL, 
	xlab='Param. Value', main=paste('Posterior distribution of\n', xlab) ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject) # or nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  for ( cIdx in 1:nChain ) {
    densInfo = density(codaObject[,c(parName)][[cIdx]]) 
    xMat = cbind(xMat,densInfo$x)
    yMat = cbind(yMat,densInfo$y)
    hdiLims = cbind(hdiLims,HDIofMCMC(codaObject[,c(parName)][[cIdx]])) ## coda provideds HPDinterval(...)
  }
  matplot( xMat , yMat , type="l" , col=plColors , 
           main=main , xlab=xlab , ylab="Density")
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims) , 0 , "95% HDI" , adj=c(0.5,-0.2) )
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  MCSE = sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth) 
  text( max(xMat) , max(yMat) , adj=c(1.0,1.0) , cex=0.8 ,
        paste(paste("MCSE =",signif(MCSE,2)), sprintf("\nESS = %.1f\n",EffChnLngth)) )

  yPos <- seq(max(yMat),min(yMat), length.out=12)
  xPos <- seq(min(xMat), max(xMat), length.out=100)
  k <- 1
  text( xPos[17], yPos[1], cex=0.8, labels="95% HDI:")
  for (i in 2:5) {
    text( xPos[18], yPos[i], cex=0.6, labels=sprintf("(%.4g~%.4g)", 
          hdiLims[1,k], hdiLims[2,k]), col=plColors[k])
    k <- k + 1
  }

}

