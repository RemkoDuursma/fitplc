#' Extract Px from fitted objects
#' @description Extract esimates of Px from an object returned by \code{\link{fitplc}}, \code{\link{fitplcs}}, \code{\link{fitcond}} or \code{\link{fitconds}}. This function allows extraction of estimates of P88 or other values when the fit estimated P50 (or other). 
#' 
#' With the Weibull model, it appears to be more robust to set \code{x=50} when fitting the curve, and extracting other points with \code{getPx}.
#' 
#' See examples for use of this function. Note that the confidence interval is based on the bootstrap resampling performed by \code{\link{fitplc}}. If the bootstrap was not performed durinf the fit (i.e. \code{boot=FALSE} in \code{fitplc} or elsewhere), it only returns the fitted values, and not the confidence intervals.
#' @param object Object returned by any of the fitting functions (e.g. \code{\link{fitplc}})
#' @param x The x in Px, that is, if P50 should be returned, x=50. Can be a vector, to return multiple points at once.
#' @param coverage The desired coverage of the confidence interval (0.95 is the default).
#' @param rescale_Px Logical (default FALSE). If TRUE, rescales calculation of Px for the sigmoidal model, by finding water potential relative to K at zero water potential (which for the sigmoidal model, is not equal to Kmax). If you fitted \code{fitcond} with \code{rescale_Px = TRUE}, make sure to set TRUE here as well to be consistent (it is not assumed from the fitted model, yet).
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
#' # Extract P88 from multiple fitted curves
#' fits <- fitplcs(stemvul, "Species", boot=FALSE)
#' getPx(fits, 88)
#' 
#'@export
getPx <- function(object, ...)UseMethod("getPx")

#'@export
getPx.default <- function(object, x=50, coverage=0.95, rescale_Px = FALSE){
  
  resc_cons <- 1
  
  getpx_fun <- function(object, x){
    X <- 1 - x/100
    
    if(object$model == "Weibull"){
      
      px <- coef(object)["PX","Estimate"]
      sx <- coef(object)["SX","Estimate"]
      v <- (object$x - 100)*log(1 - object$x/100)
      p <- px*(log(1 - x/100)/log(1 - object$x/100))^(v/(px*sx))
    } else {
      
      if(!rescale_Px){
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


#'@export
getPx.manyplcfit <- function(object,  ...){
  
  l <- lapply(object, getPx, ...)
  dfr <- cbind(data.frame(Group=names(l)),
               do.call(rbind, l))
  rownames(dfr) <- NULL
return(dfr)
}
