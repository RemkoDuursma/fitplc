
#'@export
print.plcfit <- function(x,...){
  cat("Class of object 'plcfit' as returned by 'fitplc'.\n")
  cat("-------------------------------------------------\n\n")
  
  if(x$condfit){
    cat("Fit to unscaled conductance or conductivity data.\nKmax = ",x$Kmax,"\n\n")
  }
  if(x$fitran){
    cat("Random effects estimated for ",x$ranvar,"\n")
  }
  cat("Parameters and 95% confidence interval:\n\n")
  
  print(coef(x))
  cat("\n")
}
