# This works for this object only
bootfit <- function(fit, n=999, maxnbad=50, Data, startList){
  

  
  nrow <- length(residuals(fit))
  
  refits <- list()
  ndone <- 0
  nbad <- 0
  
  while(ndone < n && nbad < maxnbad){
    bootsample <- sample(nrow, replace=TRUE)
    tryfit <- try(update(fit, subset=bootsample, data=Data, start=startList), silent=TRUE)
    if(!inherits(tryfit, "try-error")){
      ndone <- ndone + 1
      nbad <- 0
      refits[[ndone]] <- tryfit  
    } else {
      nbad <- nbad + 1
    }
  
  if(nbad == maxnbad)warning("Too many consecutive fits have failed.")
  }
  
  
  coefs <- t(sapply(refits,coef))
  return(coefs)
}