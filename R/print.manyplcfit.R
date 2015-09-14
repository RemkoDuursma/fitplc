#'@export
print.manyplcfit <- function(x,...){
  
  
  coefs <- lapply(x, coef)
  
  cat("Object of class 'manyplcfit'\n")
  cat("------------------------------\n\n")
  cat("Parameter estimates and 95% confidence intervals:\n\n")
  
  for(i in 1:length(x)){
    cat("Group: ",names(x)[i],"\n")
    print(coefs[[i]])
    cat("\n")
  }
  
}
