
#' Extract Px from fitted objects
#' @description In cases where fitplc does not converge (as is common when x=88 or other large/small values), this
#' function can be used to extract esimates of Px from a fitted object. See examples. Note that the CI is approximate, and is based on the bootstrap resampling performed by fitplc. This function only works when \code{bootci=TRUE} when the curve was fit, luckily this is the default behaviour (see also examples below).
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
getPx <- function(fit, x=50){
  
  X <- 1 - x/100
  p <- approx(x=fit$p$pred, y=fit$p$x, xout=X)$y
  lwrci <- approx(x=fit$p$lwr, y=fit$p$x, xout=X)$y
  uprci <- approx(x=fit$p$upr, y=fit$p$x, xout=X)$y
  
  vec <- c(p, lwrci, uprci)
  names(vec) <- c(paste0("P",x),"2.5%","97.5%")
return(vec)
}