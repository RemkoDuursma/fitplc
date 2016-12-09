#' Extract Px from fitted objects
#' @description Extract esimates of Px from an object returned by \code{\link{fitplc}}. This allows extraction of estimates of P88 or other values when the fit estimated P50 (or other), for example. When the Weibull model is used, it is especially recommended to fit the P50 and estimate other points of the curve with \code{getPx}. 
#' 
#' See examples for use of this function. Note that the confidence interval is based on the bootstrap resampling performed by \code{\link{fitplc}}. This function only works when \code{bootci=TRUE} when the curve was fit.
#' @param object Object returned by \code{\link{fitplc}}
#' @param x The x in Px, that is, if P50 should be returned, x=50.
#' @details Note that this function does not return a standard error, because the bootstrap confidence interval will be rarely symmetrical. If you like, you can calculate it as the mean of the half CI width (and note it as an 'approximate standard error'). Or, perhaps better, just report the CI and not the SE.
#' 
#' Also, frequently only the lower CI will be reported - sometimes the upper CI cannot be calculated (this will be more common when x is large, say for P88). In that case, assume symmetry and construct the CI with the lower confidence limit that will be reported. 
#' 
#' Finally, sometimes when x=88, it will return a missing value also for the predicted parameter. This happens when the predictions of the fitted model within the range of the data don't reach 88% embolism. We should probably avoid extrapolating the fit.
#' 
#' @examples
#' \dontrun{
#' # Make sure bootci=TRUE (the default)
#' somefit <- fitplc(dfr, x=50)
#' 
#' # Extract P12
#' getPx(somefit, x=12)
#' 
#' }
#'@export
getPx <- function(object, x=50){
  
  getpx_fun <- function(object, x){
    X <- 1 - x/100
    
    if(object$model == "Weibull"){
      
      px <- coef(object)["PX","Estimate"]
      sx <- coef(object)["SX","Estimate"]
      v <- (object$x - 100)*log(1 - object$x/100)
      p <- px*(log(1 - x/100)/log(1 - object$x/100))^(v/(px*sx))
    } else {
      p <- approx(x=object$pred$fit, y=object$pred$x, xout=X)$y
    }
    
    haveci <- "lwr" %in% names(object$pred)
    
    if(haveci){
      lwrci <- approx(x=object$pred$lwr, y=object$pred$x, xout=X)$y
      uprci <- approx(x=object$pred$upr, y=object$pred$x, xout=X)$y
      
      vec <- c(p, lwrci, uprci)
      names(vec) <- c(paste0("P",x),"2.5%","97.5%")
      
    } else {
      lwrci <- NA
      uprci <- NA
      vec <- p
      names(vec) <- paste0("P",x)
      
    }
    
  return(vec)
  }

  
  l <- lapply(x, function(val)getpx_fun(x=val, object=object))
  l <- as.data.frame(do.call(rbind,l))
  
  l <- cbind(data.frame(x=x), l)
  names(l)[1:2] <- c("x","Px")
  
return(l)
}
