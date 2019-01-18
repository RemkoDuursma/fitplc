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
  l <- list()
  l$x <- xi
  l$fit <- predict(object, dfr)
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



