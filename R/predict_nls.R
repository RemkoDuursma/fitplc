predict_nls <- function(object, xvarname=NULL, from=NULL, to=NULL, x=NULL,interval = c("none", "confidence"), 
                        level=0.95, 
                        n=101, nboot=999, add=TRUE, data, startList, weights=NULL, ...){
  
  
  interval <- match.arg(interval)
  
  if(is.null(x)){
    
    if(is.null(from) || is.null(to)){
      e <- object$m$getEnv()
      xval <- get(xvarname, envir=e)
      
      if(is.null(from))from <- min(xval)
      if(is.null(to))to <- max(xval)
    } 
    xi <- seq(from,to, length=n)
  } else {
    xi <- x
  }
  
  f <- object$m$formula()[[3]]
  
  dfr <- as.data.frame(xi)
  names(dfr) <- xvarname
  pred <- predict(object, dfr)
  l <- list()
  l$x <- xi
  l$pred <- pred
  class(l) <- "nlspred"
  
  if(interval == "confidence"){
    d <- (1-level)/2
    
    b <- bootfit(object, n=nboot, Data=data, startList=startList, weights=weights)

    parnames <- names(coef(object))
    npars <- length(parnames)
    preds <- list()
    
    for(i in seq_len(nboot)){
      p <- c(as.list(b[i,]),list(xi))
      names(p)[length(p)] <- xvarname
      preds[[i]] <- with(data, eval(f,p))
    }
    preds <- do.call(rbind,preds)
    
    l$lwr <- apply(preds,2,function(x)quantile(x, d))
    l$upr <- apply(preds,2,function(x)quantile(x, level+d))
    
    l$boot <- b
  }
  
  return(l)
   
}


plot.nls <- function(x, add=FALSE, 
                     lwd=c(1,1), lty=c(1,5), col=c("black", "red"), 
         xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,  ...){
  
  p <- predict_nls(x,...)
  
  if(!add){
    with(p, plot(rep(x,3), c(pred,upr,lwr), type='n'))
  }
  
  with(p,{
    lines(x, pred, lty=lty[1], lwd=lwd[1], col=col[1])
    if("upr" %in% names(p)){
      lines(x, upr, lty=lty[2], lwd=lwd[2], col=col[2])
      lines(x, lwr, lty=lty[2], lwd=lwd[2], col=col[2])
    }
  })
   
}

# 
# plot(x,y)
# plot(nls1, add=T)
# 
# plot(x,y,pch=19,cex=1.3,
#      panel.first=plot(nls1,add=T,lwd=2,interval="confidence",col=c("grey","red")))
# 

