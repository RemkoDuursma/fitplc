#' @export
#' @rdname fitplc
fitcond <- function(dfr, varnames = c(K="K", WP="MPa"), 
                    Kmax=NULL,
                    WP_Kmax=NULL, 
                    rescale_Px = FALSE,
                    ...){

  
  if(!rescale_Px && (is.null(WP_Kmax) && is.null(Kmax))){
    stop("Provide either the maximum hydraulic conductance (Kmax), or the water potential threshold (WP_Kmax), or set rescale_Px=TRUE (see help file).", .call=FALSE)
  }

  varnames <- unlist(varnames)
  
  # Get variables out of dataframe
  if(!varnames["K"] %in% names(dfr))
    stop("Check variable name for K!")
  if(!varnames["WP"] %in% names(dfr))
    stop("Check variable name for water potential!")
  
  K <- dfr[[varnames["K"]]]
  P <- dfr[[varnames["WP"]]]
  
  # Need absolute values of water potential
  if(mean(P) < 0)P <- -P
  if(!is.null(WP_Kmax) && WP_Kmax < 0)WP_Kmax <- -WP_Kmax
  
  # Calculate Kmax based on WP threshold
  if(is.null(Kmax) & !rescale_Px){
    Kmax <- mean(K[P < WP_Kmax], na.rm=TRUE)
  }
  if(rescale_Px){
    Kmax <- K[which.min(P)]  # Don't need Kmax, but safe to set a reasonable value.
  }
      
  # Now calculate PLC
  dfr$plc <- 100 * (1 - K/Kmax)
  
  # Now do fit
  KmaxVal <- Kmax
  f <- fitplc(dfr, varnames=c(PLC = "plc", WP = varnames[["WP"]]), 
              calledfromfitcond=TRUE, Kmax=KmaxVal,
              ...)
    
  if(rescale_Px & f$model != "loess"){
    stop("'rescale_Px' only implemented (reliably) for loess models, at the moment")
  }
  
  # Store K at WP = min(WP) (in data)
  if(f$model == "sigmoidal"){
    f$K0 <- KmaxVal * sigmoid_untrans(predict(f$fit, data.frame(minP=max(f$data$minP))))
  } else {
    f$K0 <- KmaxVal * f$pred$fit[which.min(f$pred$x)]
  }
  
  # If rescale_Px, recalculate Px relative to K0, not Kmax.
  # For linearish data, K0 << Kmax, thus Px too large.
  if(rescale_Px){
    gp <- getPx(f, x=f$x, rescale_Px = TRUE)
    f$cipars[2,] <- unname(unlist(gp)[2:4])
  }
  f$rescale_Px <- rescale_Px
  
return(f)
}


