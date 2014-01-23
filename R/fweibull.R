#' Weibull vulnerability curve
#' @description As re-parameterized by Ogle et al.
#' @export
fweibull <- function(P, SX, PX, X=50){
  
  X <- X[1] # when fitting; vector may be passed but X cannot actually vary.
  V <- (X-100)*log(1-X/100)
  p <- (P/PX)^((PX*SX)/V)
  relk <- (1-X/100)^p
  
  return(relk)
}