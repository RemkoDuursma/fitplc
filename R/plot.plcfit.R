#' Plot a fitted vulnerability curve
#' @description Standard plots of fitted curves (objects returned by \code{\link{fitplc}}, \code{\link{fitplcs}}, \code{\link{fitcond}} or \code{\link{fitconds}}), with plenty of options for customization.
#' @param xlab,ylab Optionally, X and Y axis labels (if not provided, a default is used).
#' @param ylim Optionally, Y-axis limits.
#' @param pch Optionally, the plotting symbol (default = 19, filled circles)
#' @param selines Option for the confidence interval around Px, either 'parametric' (confidence interval computed with \code{\link{confint}}), 'bootstrap' (computed with non-parametric bootstrap) or 'none' (no plotting of the confidence interval).
#' @param plotrandom Logical. If TRUE (default is FALSE), plots the predictions for the random effects (only if random effects were included in the model fit).
#' @param multiplier Multiply the scaled data (for plotting).
#' @param x A fitted curve returned by \code{fitplc}
#' @param plotPx Logical (default TRUE), whether to plot a vertical line for the P50.
#' @param plotci Logical (default TRUE), whether to plot the confidence interval (if computed with bootci=TRUE).
#' @param plotdata Logical (default TRUE), whether to add the data to the plot.
#' @param add Logical (default FALSE), whether to add the plot to a current device. This is useful to overlay two plots or curves, for example.
#' @param citype Either 'polygon' (default), or 'lines', specifying formatting of the confidence interval in the plot.
#' @param linecol The color(s) of the fitted curve (or color of the random effects curves if plotrandom=TRUE).
#' @param pointcol The color(s) of the data points.
#' @param linecol2 The color of the fixed effects curve (if plotrandom=TRUE; otherwise ignored).
#' @param pxlinecol The color of the lines indicating Px and its confidence interval 
#' @param pxcex Character size for the Px label above the Y-axis.
#' @param what Either 'relk' or 'embol'; it will plot either relative conductivity or percent embolism.
#' @param \dots Further parameters passed to \code{plot}, or \code{points} (when \code{add=TRUE})
#' @export
#' @rdname plot.plcfit
#' @importFrom graphics abline mtext plot
#' @importFrom stats approx coef confint
plot.plcfit <- function(x, xlab=NULL, ylab=NULL, ylim=NULL, pch=19, 
                        plotPx=TRUE, plotci=TRUE, plotdata=TRUE, add=FALSE,
                        multiplier=NULL,
                        selines=c("parametric","bootstrap","none"),
                        plotrandom=FALSE,linecol="black",
                        linecol2="blue",
                        pxlinecol="red",
                        pxcex=0.7,
                        citype=c("polygon","lines"),
                        what=c("relk","embol"), ...){
  
  
  if(x$model == "sigmoidal")plotPx <- FALSE # temporary!!
  
  if(is.null(multiplier)){
    multiplier <- x$Kmax
  }
  
  selines <- match.arg(selines)
  citype <- match.arg(citype)
  
  if(is.null(xlab))xlab <- expression(Water~potential~~(-MPa))
  
  type <- ifelse(plotdata, 'p', 'n')
  what <- match.arg(what)
  
  # override
  if(x$condfit){
    what <- "relk"
  }
  
  if(plotrandom && !x$fitran)
    stop("To plot random effects predictions, refit with 'random' argument.")
  
  if(what == "relk"){
    if(is.null(ylab)){
      if(!x$condfit){
        ylab <- "Relative conductivity (0 - 1)"
      } else {
        ylab <- "Conductivity / conductance (in units provided)"
      }
    }
      
    x$data$Y <- x$data$relK
    if(is.null(ylim))ylim <- c(0,multiplier*max(x$data$Y))
  }
  
  
  if(what == "embol"){
    
    if(is.null(ylab))ylab <- "% Embolism"
    
    x$data$Y <- relk_to_plc(x$data$relK)
    if(x$bootci){
      x$pred$lwr <- toEmbol(x$pred$lwr)
      x$pred$upr <- toEmbol(x$pred$upr)
    }
    x$pred$pred <- toEmbol(x$pred$pred)
    if(is.null(ylim))ylim <- c(0,100)
    
    if(x$fitran && plotrandom){
      
      ng <- length(x$prednlme)
      for(i in 1:ng){
        x$prednlme[[i]]$y <- toEmbol(x$prednlme[[i]]$y)
      }
      x$prednlmefix$y <- toEmbol(x$prednlmefix$y)
    }
    
  } else {
    x$data$Y <- x$data$relK
  }
  
  if(!add){
    with(x, {
      plot(data$P, multiplier * data$Y, ylim=ylim, pch=pch,
           xlab=xlab,
           type=type,
           ylab=ylab, ...)
    })
  } else {
    with(x, {
      points(data$P, multiplier * data$Y, pch=pch, type=type,...)
    })
  }
  if(!plotrandom){
    if(x$bootci && plotci){
      if(citype == "lines"){
        with(x$pred,{
          lines(x, multiplier * lwr, type='l', lty=5, col=linecol)
          lines(x, multiplier * upr, type='l', lty=5, col=linecol)
        })
      }
      if(citype == "polygon"){
        with(x$pred, addpoly(x,multiplier * lwr,multiplier * upr))
        # replot points
        if(plotdata){
          with(x, {
            points(data$P, multiplier * data$Y, pch=pch, type=type,...)
          })
        }
      }
      
    }
    with(x$pred,{
      lines(x, multiplier * pred, type='l', lty=1, col=linecol)
    })
  }
  
  if(plotrandom){
    for(i in 1:length(x$prednlme)){
      with(x$prednlme[[i]], lines(x,multiplier * y,type='l',col=linecol))
    }  
    with(x$prednlmefix, lines(x,multiplier * y,type='l',lwd=2, col=linecol2))
  }
  
  if(plotPx){
    if(!x$fitran){
      px <- coef(x$fit)["PX"]
      
      if(selines == "bootstrap"){
        px_ci <- x$bootpars[2,2:3]
      } 
      if(selines == "parametric") {
        px_ci <- x$ci[2,]
      }
      
    } else {
      px <- fixef(x$nlmefit)["PX"]
      px_ci <- x$cinlme[2,]
    }
    
    abline(v=px, col=pxlinecol)
    if(selines != "none")abline(v=px_ci, col=pxlinecol, lty=5)
    
    mtext(side=3, at=px, text=bquote(P[.(x$x)]), 
          line=0, col=pxlinecol, cex=pxcex)
  }
  
}


#'@export
#'@param onepanel For plotting of many curve fits, plot all curves in one panel (TRUE) or in separate panels (FALSE)
#'@param legend Logical (default TRUE), whether to include a simple legend when plotting multiple fits
#'@param legendwhere As in \code{\link{legend}}, specification of where to place legend (e.g. 'bottomleft'; coordinates not accepted)
#'@rdname plot.plcfit
#'@importFrom grDevices rainbow
plot.manyplcfit <- function(x, what=c("relk","embol"), 
                            onepanel=FALSE, linecol=NULL, 
                            pointcol=NULL,
                            pch=19, 
                            legend=TRUE, legendwhere="topright", ...){
  
  what <- match.arg(what)
  np <- length(x)
  
  if(length(pch) < np)pch <- rep(pch,np)
  
  if(!onepanel){
    if(is.null(linecol) | length(linecol) < np)linecol <- rep("black",np)
    if(is.null(pointcol) | length(pointcol) < np)pointcol <- rep("black",np)
    
    for(i in 1:np)plot(x[[i]], linecol=linecol[i], pch=pch[i], what=what, ...)
  } else {
    if(is.null(linecol))linecol <- rainbow(np)
    if(is.null(pointcol))pointcol <- rainbow(np)
    
    plot(x[[1]], pch=pch[1], col=pointcol[1], what=what,...)
    if(np > 1){
      for(i in 2:np){
        plot(x[[i]], add=TRUE, linecol=linecol[i], col=pointcol[i], pch=pch[i], what=what, ...)
      }
    }
    # If plotting lines, plot them again to make sure they are on top
    for(i in 1:np){
      plot(x[[i]], add=TRUE, linecol=linecol[i], 
                       what=what,
                       plotPx=FALSE, 
                       plotdata=FALSE, plotci=FALSE)
    }
    if(legend){
      legend(legendwhere, names(x), lty=1, col=linecol)
    }
    
  }
}







