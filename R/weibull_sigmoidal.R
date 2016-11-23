#' Weibull vulnerability curve
#' @description The Weibull function, as re-parameterized by Ogle et al. (2009), which describes the relative conductivity as a function of the plant water potential. The relative conductivity (relK) is scaled to be 1 when water potential is zero. The slope of relK versus P at the inflection point can be calculated from the shape parameter (SX) as slope = -SX/100.
#' @param P Water potential (positive-valued MPa)
#' @param SX Shape parameter
#' @param PX Water potential at X loss of conductivity.
#' @param X If 50, PX is the P50.
#' @references Ogle, K., Barber, J.J., Willson, C., Thompson, B., 2009. Hierarchical statistical modeling of xylem vulnerability to cavitation. New Phytologist 182, 541-554. 
#' @examples
#' curve(fweibull(x, SX=30, PX=2), from=0, to=5)
#' @export
fweibull <- function(P, SX, PX, X=50){
  
  X <- X[1] # when fitting; vector may be passed but X cannot actually vary.
  V <- (X-100)*log(1-X/100)
  p <- (P/PX)^((PX*SX)/V)
  relk <- (1-X/100)^p
  
  return(relk)
}

#' Sigmoidal vulnerability curve
#' @description A sigmoidal-exponential function, which describes the relative conductivity as a function of the plant water potential. The relative conductivity is scaled to be 1 when water potential is zero. This function was used by Pammenter and vander Willigen (1998), but note that this implementation gives the relative conductivity, not the PLC (but relK = 1 - PLC). The slope of relK versus P at the inflection point can be calculated from the shape parameter (a) as slope = -a/4.
#' @param P Water potential (positive-valued MPa)
#' @param a Shape parameter, related to the slope at the inflection point (see Description).
#' @param PX Water potential at X loss of conductivity (positive valued).
#' @param X If 50, PX is the P50.
#' @references Pammenter, N.W., Willigen, C.V. der, 1998. A mathematical and statistical analysis of the curves illustrating vulnerability of xylem to cavitation. Tree Physiol 18, 589-593. doi:10.1093/treephys/18.8-9.589
#' @examples
#' curve(fsigmoidal(x, PX=-2, a=5), from=0, to=-5)
#' curve(fsigmoidal(x, PX=-2, a=2), add=T)
#' 
#' # Comparison to Weibull
#' curve(fweibull(x, PX=3, SX=40), from=0, to=6)
#' curve(fsigmoidal(x, PX=3, a=4*(40/100)), add=TRUE, col="red")
fsigmoidal <- function(P, PX, a, X=50){
  
  X <- X[1] # vector might have been passed but X cannot actually vary.
  P <- -P
  PX <- -PX
  b <- PX - (1/a)*(50/X - 1)
  
  1 - 1/(1 + exp(a*(P - b)))
}


