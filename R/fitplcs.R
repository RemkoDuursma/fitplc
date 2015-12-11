#' @export
#' @rdname fitplc
#' @param group Character; variable in the dataframe that specifies groups. The curve will be fit for every group level.
#' @param \dots Further parameters passed to \code{fitplc}.
fitplcs <- function(dfr, group, ...){
  
  if(!group %in% names(dfr))
    stop("You must provide a name in the dataframe to fit by.")
  
  dfrs <- split(dfr, dfr[,group])
  
  fits <- lapply(dfrs, function(x)fitplc(x, ...))
  class(fits) <- "manyplcfit"
  
  return(fits)
}


#' @rdname fitplc
#' @export
fitconds <- function(dfr, group, ...){
  
  if(!group %in% names(dfr))
    stop("You must provide a name in the dataframe to fit by.")
  
  dfrs <- split(dfr, dfr[,group])
  
  fits <- lapply(dfrs, function(x)fitcond(x, ...))
  
  return(fits)
  
}

