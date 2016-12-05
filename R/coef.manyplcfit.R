#'@export
coef.manyplcfit <- function(object, ...){
  
  x <- do.call(rbind, lapply(object,coef))
  rn <- rownames(x)
  rownames(x) <- NULL
  dfr <- cbind(data.frame(Group=rep(names(object), each=2), Parameter=rn),
               as.data.frame(x))
  
  return(dfr)
}
