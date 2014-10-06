#'@export
coef.manyplcfit <- function(object, ...){
  
  x <- lapply(object,coef)
  vf <- function(m)as.vector(t(m))
  dfr <-  as.data.frame(do.call(rbind,lapply(x,vf)))
  names(dfr) <- c("S","S_SE","S_lower","S_upper",
                  "Px","Px_SE","P50_lower","P50_upper")
  
  return(dfr)
}
