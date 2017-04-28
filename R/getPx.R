#' Extract Px from fitted objects
#' @description Extract esimates of Px from an object returned by \code{\link{fitplc}}. This allows extraction of estimates of P88 or other values when the fit estimated P50 (or other), for example. When the Weibull model is used, it is especially recommended to fit the P50 and estimate other points of the curve with \code{getPx}. 
#' 
#' See examples for use of this function. Note that the confidence interval is based on the bootstrap resampling performed by \code{\link{fitplc}}. This function only works when \code{bootci=TRUE} when the curve was fit.
#' @param object Object returned by \code{\link{fitplc}}
#' @param x The x in Px, that is, if P50 should be returned, x=50.
#' @param coverage The desired coverage of the confidence interval (0.95 is the default).
#' @details Note that this function does not return a standard error, because the bootstrap confidence interval will be rarely symmetrical. If you like, you can calculate it as the mean of the half CI width (and note it as an 'approximate standard error'). A better approach is to only report the CI and not the SE.
#' 
#' Sometimes the upper CI cannot be calculated and will be reported as \code{NA}. This indicates that the upper confidence bound is outside the range of the data, and can therefore not be reliably reported. It is especially common when \code{x} is large, say for P88. 
#' 
#' @examples
#' # A fit
#' somefit <- fitplc(stemvul, x=50, model="sigmoid")
#' 
#' # Extract P12, P88
#' # Note NA for upper CI for P88; this is quite common
#' # and should be interpreted as falling outside the range of the data.
#' getPx(somefit, x=c(12,88))
#' 
#'@export
getPx <- function(object, x=50, coverage=0.95, sigmoid_rescale_Px = FALSE){
  
  resc_cons <- 1
  
  getpx_fun <- function(object, x){
    X <- 1 - x/100
    
    if(object$model == "Weibull"){
      
      px <- coef(object)["PX","Estimate"]
      sx <- coef(object)["SX","Estimate"]
      v <- (object$x - 100)*log(1 - object$x/100)
      p <- px*(log(1 - x/100)/log(1 - object$x/100))^(v/(px*sx))
    } else {
      
      if(!sigmoid_rescale_Px){
        p <- approx(x=object$pred$fit, y=object$pred$x, xout=X)$y
      } else {
        resc_cons <- object$Kmax / object$K0
        p <- approx(x=object$pred$fit * resc_cons, y=object$pred$x, xout=X)$y
      }
    }
    
    haveci <- "lwr" %in% names(object$pred)
    
    if(haveci){
      lwrci <- approx(x=object$pred$lwr * resc_cons, y=object$pred$x, xout=X)$y
      uprci <- approx(x=object$pred$upr * resc_cons, y=object$pred$x, xout=X)$y
      
      vec <- c(p, lwrci, uprci)
      names(vec) <- c(paste0("P",x),label_lowci(coverage), label_upci(coverage))
      
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
