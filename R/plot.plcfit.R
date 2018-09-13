#' Plot a fitted vulnerability curve
#' @description Standard plots of fitted curves (objects returned by \code{\link{fitplc}}, \code{\link{fitplcs}}, \code{\link{fitcond}} or \code{\link{fitconds}}), with plenty of options for customization.
#' @param xlab,ylab Optionally, X and Y axis labels (if not provided, a default is used).
#' @param ylim Optionally, Y-axis limits.
#' @param pch Optionally, the plotting symbol (default = 19, filled circles)
#' @param px_ci Option for the confidence interval around Px, either 'parametric' (confidence interval computed with \code{\link{confint}}), 'bootstrap' (computed with non-parametric bootstrap) or 'none' (no plotting of the confidence interval) (formerly argument was called \code{selines})
#' @param px_ci_type Either 'vertical' (default), or 'horizontal', to plot confidence limits for Px.
#' @param px_ci_label Logical (default TRUE), whether to write a label next to the CI for Px.
#' @param plotrandom Logical. If TRUE (the default is FALSE), plots the predictions for the random effects (only if random effects were included in the model fit).
#' @param multiplier Multiply the scaled data (for plotting).
#' @param x A fitted curve returned by \code{fitplc}
#' @param plotPx Logical (default TRUE), whether to plot a vertical line for the P50.
#' @param plotci Logical (default TRUE), whether to plot the confidence interval (if computed with bootci=TRUE).
#' @param plotdata Logical (default TRUE), whether to add the data to the plot.
#' @param plotfit Logical (default TRUE), whether to add the fitted curve to the plot.
#' @param add Logical (default FALSE), whether to add the plot to a current device. This is useful to overlay two plots or curves, for example.
#' @param citype Either 'polygon' (default), or 'lines', specifying formatting of the confidence interval in the plot.
#' @param linecol The color(s) of the fitted curve (or color of the random effects curves if plotrandom=TRUE).
#' @param linetype Line type for fitted curve (see options for \code{lty} in \code{\link{par}}).
#' @param linelwd Width of the line (see options for \code{lwd} in \code{\link{par}}).
#' @param pointcol The color(s) of the data points.
#' @param linecol2 The color of the fixed effects curve (if plotrandom=TRUE; otherwise ignored).
#' @param cicol The color of the confidence interval band (if plotted).
#' @param pxlinecol The color of the lines indicating Px and its confidence interval 
#' @param pxcex Character size for the Px label above the Y-axis.
#' @param what Either 'relk' or 'PLC' (or synonym 'embol'); it will plot either relative conductivity or percent loss conductivity (percent embolism).
#' @param selines Obsolete; use \code{px_ci}
#' @param xaxis Either 'positive' (default), so that water potential is plotted as positive values, or 'negative', plotting negative-valued water potentials.
#' @param \dots Further parameters passed to \code{plot}, or \code{points} (when \code{add=TRUE})
#' @export
#' @rdname plot.plcfit
#' @importFrom graphics abline mtext plot box
#' @importFrom stats approx coef confint
plot.plcfit <- function(x, xlab=NULL, ylab=NULL, ylim=NULL, pch=19, 
                        plotPx=TRUE, plotci=TRUE, plotdata=TRUE, plotfit=TRUE, add=FALSE,
                        multiplier=NULL,
                        px_ci=c("bootstrap","parametric","none"),
                        px_ci_type=c("vertical","horizontal"),
                        px_ci_label=TRUE,
                        plotrandom=FALSE,
                        pointcol="black",
                        linecol="black",
                        linetype=1,
                        linelwd=1,
                        linecol2="blue",
                        pxlinecol="red",
                        pxcex=0.7,
                        citype=c("polygon","lines"),
                        cicol=alpha("lightgrey",0.8),
                        what=c("relk","PLC","embol"), 
                        selines=NULL,
                        xaxis=c("positive","negative"),
                        ...){
  
  
  if(missing(multiplier)){
    multiplier <- x$Kmax
  }
  
  if(!missing(selines)){
    warning("Argument 'selines' is now called 'px_ci'.")
    px_ci <- match.arg(selines, eval(formals(plot.plcfit)$px_ci))
  } else {
    px_ci <- match.arg(px_ci)
  }
  
  citype <- match.arg(citype)
  px_ci_type <- match.arg(px_ci_type)
  
  xaxis <- match.arg(xaxis)
  xsign <- if(xaxis == "negative")-1 else 1
  
  type <- ifelse(plotdata, 'p', 'n')
  
  what <- match.arg(what)
  if(what == "embol")what <- "PLC"
  
  # override
  if(x$condfit){
    what <- "relk"
  }
  
  if(plotrandom && !x$fitran){
    stop("To plot random effects predictions, refit with 'random' argument.")
  }
  
  # Set x-axis label
  if(missing(xlab)){
    if(xaxis == "positive"){
      xlab <- expression(Water~potential~~(-MPa))
    } else {
      xlab <- expression(Water~potential~~(MPa))
    }
  }
  
  # Set y-axis label
  if(missing(ylab)){
    if(what == "relk"){
      
      if(!x$condfit){
        ylab <- "Relative conductivity (0 - 1)"
      } else {
        ylab <- "Conductivity / conductance (in units provided)"
      }
    
    } else if(what == "PLC"){
      ylab <- "Percent loss conductivity (%)"
    }
  }
  
  # Set data
  if(what == "relk"){
    x$data$Y <- x$data$relK
  } else if(what == "PLC"){
    x$data$Y <- x$data$PLC
    x$pred$fit <- relk_to_plc(x$pred$fit)
    if(x$bootci){
      x$pred$lwr <- relk_to_plc(x$pred$lwr)
      x$pred$upr <- relk_to_plc(x$pred$upr)
    }
  }
  
  # Set y-axis limit
  if(is.null(ylim))ylim <- c(0,multiplier*max(x$data$Y))

  if(x$fitran && plotrandom){
    
    ng <- length(x$pred$ran)
    if(what=="PLC"){
      
      x$pred$ran <- lapply(x$pred$ran, 
                           function(f){
                             f$fit <- relk_to_plc(f$fit)
                             return(f)
                           })
      
      if(!x$fitran)x$pred$fit <- relk_to_plc(x$pred$fit)
      
    }
  }

  if(!add){
    with(x, {
      plot(xsign * data$P, multiplier * data$Y, ylim=ylim, pch=pch,
           col=pointcol,
           xlab=xlab,
           type=type,
           ylab=ylab, ...)
    })
  } else {
    with(x, {
      points(xsign * data$P, multiplier * data$Y, pch=pch, 
             type=type, col=pointcol, ...)
    })
  }
  
  if(!plotrandom){
    if(x$bootci && plotci){
      if(citype == "lines"){
        with(x$pred,{
          lines(xsign * x, multiplier * lwr, type='l', lty=5, col=linecol)
          lines(xsign * x, multiplier * upr, type='l', lty=5, col=linecol)
        })
      }
      if(citype == "polygon"){
        with(x$pred, addpoly(xsign * x, multiplier * lwr, multiplier * upr, col=cicol))
        
        # replot points
        if(plotdata){
          with(x, {
            points(xsign * data$P, multiplier * data$Y, pch=pch, type=type, col=pointcol, ...)
          })
        }
        
        # replot box
        box()
      }
      
    }
    if(plotfit){
      with(x$pred,{
        lines(xsign * x, multiplier * fit, type='l', lty=linetype, col=linecol, lwd=linelwd)
      })
    }
  }
  
  if(plotrandom){
    
    for(i in seq_along(x$pred$ran)){
      with(x$pred$ran[[i]], 
           lines(xsign * x, multiplier * fit, 
                 type='l', col=linecol)
           )
    } 
    with(x$pred, lines(xsign * x, multiplier * fit, type='l', 
                       lwd=2, col=linecol2))
  }
  
  if(plotPx){
    
    px <- coef(x)["PX","Estimate"]
    
    if(px_ci_type != "horizontal"){
      abline(v=xsign * px, col=pxlinecol)
      mtext(side=3, at=xsign * px, text=bquote(P[.(x$x)]), 
              line=0, col=pxlinecol, cex=pxcex)
    }
    
    haveboot <- any(grepl("Boot", colnames(coef(x))))
    havenorm <- any(grepl("Norm", colnames(coef(x))))
    
    if(px_ci == "bootstrap" && !haveboot)px_ci <- "parametric"
    
    # Confidence lines for the P50
    if(px_ci != "none"){
      if(px_ci == "bootstrap" && !haveboot)px_ci <- "parametric"
      if(px_ci == "parametric" && !havenorm)px_ci <- "bootstrap"
      
      nm <- switch(px_ci, bootstrap="Boot", parametric="Norm")
      px_ci <- coef(x)["PX",ci_names(nm,coverage=x$coverage)]
      
      u <- par()$usr
      dx <- (u[2] - u[1])/150
      
      if(px_ci_type == "vertical"){
        abline(v=xsign * px_ci, col=pxlinecol, lty=5)

        if(px_ci_label){
          lab <- paste(label_coverage(x$coverage),nm)
          text(xsign * (px_ci[2]+dx), u[3] + 0.96*(u[4] - u[3]),
               lab, cex=0.5*par()$cex, pos=4)
        }
        
      }
      if(px_ci_type == "horizontal") {
        segments(x0=xsign*px_ci[1], x1=xsign*px_ci[2],
                 y0=1-x$x/100, y1=1-x$x/100)
        points(x=xsign * px, y=1-x$x/100, pch=21, bg="white")
        
        if(px_ci_label){
          lab <- bquote(P[.(x$x)])
          text(xsign*(px_ci[2]+dx), 1-x$x/100,
               lab, cex=0.6*par()$cex, pos=4)
        }
      }
      
    }
    
  }
}


#'@export
#'@param onepanel For plotting of many curve fits, plot all curves in one panel (TRUE) or in separate panels (FALSE)
#'@param legend (for fitconds and fitplcs only) Logical (default TRUE), whether to include a simple legend when plotting multiple fits
#'@param legendwhere (for fitconds and fitplcs only) As in \code{\link{legend}}, specification of where to place legend (e.g. 'bottomleft'; coordinates not accepted)
#'@rdname plot.plcfit
#'@importFrom grDevices rainbow
plot.manyplcfit <- function(x, what=c("relk","embol","PLC"), 
                            onepanel=FALSE, 
                            linecol=NULL, 
                            pointcol=NULL,
                            pch=19, 
                            legend=TRUE, 
                            legendwhere="topright", ...){
  
  what <- match.arg(what)
  if(what == "embol")what <- "PLC"
  np <- length(x)
  
  if(length(pch) < np)pch <- rep(pch,np)
  
  if(!onepanel){
    if(is.null(linecol) | length(linecol) < np)linecol <- rep("black",np)
    if(is.null(pointcol) | length(pointcol) < np)pointcol <- rep("black",np)
    
    for(i in 1:np)plot(x[[i]], linecol=linecol[i], pch=pch[i], what=what, ...)
  } else {
    if(is.null(linecol))linecol <- rainbow(np)
    if(is.null(pointcol))pointcol <- rainbow(np)
    
    plot(x[[1]], pch=pch[1], pointcol=pointcol[1], linecol=linecol[1], what=what,...)
    if(np > 1){
      for(i in 2:np){
        plot(x[[i]], add=TRUE, linecol=linecol[i], pointcol=pointcol[i], pch=pch[i], what=what, ...)
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







