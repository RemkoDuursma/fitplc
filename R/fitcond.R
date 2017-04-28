#' @export
#' @rdname fitplc
fitcond <- function(dfr, varnames = c(K="K", WP="MPa"), 
                    Kmax=NULL,
                    WP_Kmax=NULL, 
                    rescale_Px = FALSE,
                    ...){

  
  if(is.null(WP_Kmax) && is.null(Kmax)){
    stop("Provide either the maximum hydraulic conductance (Kmax), or the water potential threshold (WP_Kmax) (see help file).", .call=FALSE)
  }

  if(is.list(varnames))varnames <- unlist(varnames)
  
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
  if(is.null(Kmax)){
    Kmax <- mean(K[P < WP_Kmax], na.rm=TRUE)
  }
      
  # Now calculate PLC
  dfr$plc <- 100 * (1 - K/Kmax)
  
  # Now do fit
  KmaxVal <- Kmax
  f <- fitplc(dfr, varnames=c(PLC = "plc", WP = varnames[["WP"]]), 
              calledfromfitcond=TRUE, Kmax=KmaxVal,
              ...)
  
  # Store K at WP = 0
  f$K0 <- KmaxVal * sigmoid_untrans(predict(f$fit, data.frame(minP=0)))
  
  # If rescale_Px, recalculate Px relative to K0, not Kmax.
  # For linearish data, K0 << Kmax, thus Px too large.
  if(rescale_Px){
    gp <- getPx(f, x=f$x, sigmoid_rescale_Px = TRUE)
    f$cipars[2,] <- unname(unlist(gp)[2:4])
  }
  
return(f)
}


