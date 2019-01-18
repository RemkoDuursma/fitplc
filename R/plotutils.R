

#' @importFrom graphics polygon
#' @importFrom grDevices col2rgb
addpoly <- function(x,y1,y2,col="#D3D3D3CC",...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}

