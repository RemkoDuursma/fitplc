
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
  cat("Parameters and %s%% confidence interval:\n\n", label_coverage(x$coverage))
  
  print(coef(x))
  cat("\n")
}
