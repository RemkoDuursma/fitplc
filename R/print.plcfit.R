
#'@export
print.plcfit <- function(x,...){
  cat("Class of object 'plcfit' as returned by 'fitplc'.\n")
  cat("-------------------------------------------------\n\n")
  if(x$fitran){
    cat("Random effects estimated for ",x$ranvar,"\n")
  }
  cat("Parameters, SE, and 95% confidence interval:\n\n")
  cat()
  
  print(coef(x))
  cat("\n")
}