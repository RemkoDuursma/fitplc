#'@rdname fitplc
#'@export


#'@rdname fitplc
#'@export
plot.plcfit <- function(x, xlab=NULL, ylab=NULL, ylim=NULL, pch=19, 
                        plotPx=TRUE, plotci=TRUE, plotdata=TRUE, add=FALSE,
                        selines=c("parametric","bootstrap"),
                        plotrandom=FALSE,linecol="black",
                        linecol2="blue",
                        pxlinecol="red",
                        pxcex=0.7,
                        citype=c("polygon","lines"),
                        what=c("relk","embol"), ...){
  
  
  selines <- match.arg(selines)
  citype <- match.arg(citype)
  
  if(is.null(xlab))xlab <- expression(Water~potential~~(-MPa))
  
  type <- ifelse(plotdata, 'p', 'n')
  what <- match.arg(what)
  
  if(plotrandom && !x$fitran)
    stop("To plot random effects predictions, refit with 'random' argument.")
  
  if(what == "relk"){
    if(is.null(ylab))ylab <- "Relative conductivity (0 - 1)"
    x$data$Y <- x$data$relK
    if(is.null(ylim))ylim <- c(0,1)
  }
  toEmbol <- function(k)100 - 100*k
  if(what == "embol"){
    if(is.null(ylab))ylab <- "% Embolism"
    
    x$data$Y <- toEmbol(x$data$relK)
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
    
  }
  
  if(!add){
    with(x, {
      plot(data$P, data$Y, ylim=ylim, pch=pch,
           xlab=xlab,
           type=type,
           ylab=ylab, ...)
    })
  } else {
    with(x, {
      points(data$P, data$Y, pch=pch, type=type,...)
    })
  }
  if(!plotrandom){
    if(x$bootci && plotci){
      if(citype == "lines"){
        with(x$pred,{
          lines(x, lwr, type='l', lty=5, col=linecol)
          lines(x, upr, type='l', lty=5, col=linecol)
        })
      }
      if(citype == "polygon"){
        with(x$pred, addpoly(x,lwr,upr))
        # replot points
        if(plotdata){
          with(x, {
            points(data$P, data$Y, pch=pch, type=type,...)
          })
        }
      }
      
    }
    with(x$pred,{
      lines(x, pred, type='l', lty=1, col=linecol)
    })
  }
  
  if(plotrandom){
    for(i in 1:length(x$prednlme)){
      with(x$prednlme[[i]], lines(x,y,type='l'))
    }  
    with(x$prednlmefix, lines(x,y,type='l',lwd=2, col=linecol2))
  }
  
  if(plotPx){
    if(!x$fitran){
      px <- coef(x$fit)["PX"]
      
      if(selines == "bootstrap"){
        px_ci <- x$bootpars[2,2:3]
      } else {
        px_ci <- x$ci[2,]
      }
    } else {
      px <- fixef(x$nlmefit)["PX"]
      px_ci <- x$cinlme[2,]
    }
    
    abline(v=px, col=pxlinecol)
    abline(v=px_ci, col=pxlinecol, lty=5)
    mtext(side=3, at=px, text=bquote(P[.(x$x)]), 
          line=0, col=pxlinecol, cex=pxcex)
  }
  
}