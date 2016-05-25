library(sm)
vioplot2 <- function(datas,range=1.5,h=NULL,ylim=NULL,names=NULL, horizontal=FALSE, 
  col="magenta", border=NA, lty=1, lwd=1, rectCol="black", colMed="white", pchMed=19, at, add=FALSE, wex=0.8, 
  drawRect=TRUE, box=TRUE, axes=TRUE, color.brewer=TRUE)
{
	color.generator <- function(color.name, i) brewer.pal(8,color.name)[i]
	if (color.brewer) {
		col <- sapply(rownames(brewer.pal.info), color.generator, 1)
		#border <- sapply(rownames(brewer.pal.info), color.generator, 2)
		rectCol <- sapply(rownames(brewer.pal.info), color.generator, 3)
		colMed <- sapply(rownames(brewer.pal.info), color.generator, 5)
	}

	if (is.data.frame(datas)) {
		labeled.names <- colnames(datas)
		datas <- as.list(datas)
	} else if (is.list(datas)) {
		labeled.names <- names(datas)
	} 
	else {
		stop("Input object is NOT data frame!")
	}
    n <- length(datas)
   
    if(missing(at)) at <- 1:n
    
    # pass 1
    #
    # - calculate base range
    # - estimate density
    #
    
    # setup parameters for density estimation
    upper  <- vector(mode="numeric",length=n)
    lower  <- vector(mode="numeric",length=n) 
    q1     <- vector(mode="numeric",length=n)
    q3     <- vector(mode="numeric",length=n)
    med    <- vector(mode="numeric",length=n)
    base   <- vector(mode="list",length=n)
    height <- vector(mode="list",length=n)
    baserange <- c(Inf,-Inf)
    
    # global args for sm.density function-call   
    args <- list(display="none")
       
    if (!(is.null(h)))
        args <- c(args, h=h)
   
            
    for(i in 1:n) {

        data<-datas[[i]]
        
        # calculate plot parameters
        #   1- and 3-quantile, median, IQR, upper- and lower-adjacent
        
        data.min <- min(data)
        data.max <- max(data)
        q1[i]<-quantile(data,0.25)
        q3[i]<-quantile(data,0.75)
        med[i]<-median(data)
        iqd <- q3[i]-q1[i]
        upper[i] <- min( q3[i] + range*iqd, data.max )
        lower[i] <- max( q1[i] - range*iqd, data.min )
        
       
        #   strategy:
        #       xmin = min(lower, data.min))
        #       ymax = max(upper, data.max))
        #

        est.xlim <- c( min(lower[i], data.min), max(upper[i], data.max) ) 
        
        # estimate density curve
        
        smout <- do.call("sm.density", c( list(data, xlim=est.xlim), args ) )

        
        # calculate stretch factor
        #
        #  the plots density heights is defined in range 0.0 ... 0.5 
        #  we scale maximum estimated point to 0.4 per data
        #
        
        hscale <- 0.4/max(smout$estimate) * wex
        
        
        # add density curve x,y pair to lists
        
        base[[i]]   <- smout$eval.points
        height[[i]] <- smout$estimate * hscale
        
        
        # calculate min,max base ranges
        
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1],t[1])
        baserange[2] <- max(baserange[2],t[2])

    }
   
    # pass 2
    #
    # - plot graphics
    
    # setup parameters for plot

    if(!add){
      xlim <- if(n==1) 
               at + c(-.5, .5)
              else 
               range(at) + min(diff(at))/2 * c(-1,1)
    
      if (is.null(ylim)) {
         ylim <- baserange
      }
    }
    if (is.null(names)) {
        #label <- 1:n
		label <- labeled.names
    } else {
        label <- names
    }

    boxwidth <- 0.05 * wex
    
        
    # setup plot

    if(!add)
      plot.new()
    if(!horizontal) {
      if(!add){
        plot.window(xlim = xlim, ylim = ylim)
		if (axes) {
         axis(2, las=1)
         axis(1,at = at, label=label )
		}
      }  
      
      if (box) { box(); }
      
      for(i in 1:n) {
       
          # plot left/right density curve
          if (!color.brewer) {
          polygon( c(at[i]-height[[i]], rev(at[i]+height[[i]])), c(base[[i]], rev(base[[i]])), col = col, border=border, lty=lty, lwd=lwd)
		  } else {
          polygon( c(at[i]-height[[i]], rev(at[i]+height[[i]])), c(base[[i]], rev(base[[i]])), col = col[i], border=border, lty=lty, lwd=lwd)
          }
        
        
          if(drawRect){
            # plot IQR
            if (!color.brewer) { lines( at[c( i, i)], c(lower[i], upper[i]) ,lwd=lwd, lty=lty) }
            if (color.brewer) { lines( at[c( i, i)], c(lower[i], upper[i]) ,lwd=lwd, lty=lty, col=rectCol[i]) }
        
            # plot 50% KI box
        
            #rect( at[i]-boxwidth/2, q1[i], at[i]+boxwidth/2, q3[i], col=rectCol)
            if (!color.brewer) { rect( at[i]-boxwidth/2, q1[i], at[i]+boxwidth/2, q3[i], col=rectCol) }
            if (color.brewer) { rect( at[i]-boxwidth/2, q1[i], at[i]+boxwidth/2, q3[i], col=rectCol[i], border=rectCol[i]) }
        
            # plot median point
        
            if (!color.brewer) { points( at[i], med[i], pch=pchMed, col=colMed ) }
            if (color.brewer) { points( at[i], med[i], pch=pchMed, col=colMed[i] ) }
         }
      }

    }
    else {
      if(!add){
        plot.window(xlim = ylim, ylim = xlim)
        axis(1)
        axis(2,at = at, label=label )
      }
      
      box()
      for(i in 1:n) {
       
          # plot left/right density curve
        
          if (!color.brewer) { polygon( c(base[[i]], rev(base[[i]])), c(at[i]-height[[i]], rev(at[i]+height[[i]])), col = col, border=border, lty=lty, lwd=lwd) }
          if (color.brewer) { polygon( c(base[[i]], rev(base[[i]])), c(at[i]-height[[i]], rev(at[i]+height[[i]])), col = col[i], border=border, lty=lty, lwd=lwd) }
        
        
          if(drawRect){
            # plot IQR
            if (!color.brewer) lines( c(lower[i], upper[i]), at[c(i,i)] ,lwd=lwd, lty=lty)
            if (color.brewer) lines( c(lower[i], upper[i]), at[c(i,i)] ,lwd=lwd, lty=lty, col=rectCol[i])
        
            # plot 50% KI box
        
            if (!color.brewer) rect( q1[i], at[i]-boxwidth/2, q3[i], at[i]+boxwidth/2,  col=rectCol)
            if (color.brewer) rect( q1[i], at[i]-boxwidth/2, q3[i], at[i]+boxwidth/2,  col=rectCol[i], border=rectCol[i])
        
            # plot median point
            if (!color.brewer) points( med[i], at[i], pch=pchMed, col=colMed )
            if (color.brewer) points( med[i], at[i], pch=pchMed, col=colMed[i] )
          }
      }

      
    }

    
    invisible (list( upper=upper, lower=lower, median=med, q1=q1, q3=q3))
}

