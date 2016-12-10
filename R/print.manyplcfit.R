#'@export
print.manyplcfit <- function(x,...){
  
  
  coefs <- lapply(x, coef)
  
  cat("Object of class 'manyplcfit'\n")
  cat("------------------------------\n\n")
  cat(sprintf("Parameter estimates and %s%% confidence intervals:\n\n", label_coverage(x$coverage)))
  
  for(i in 1:length(x)){
    cat("Group: ",names(x)[i],"\n")
    print(coefs[[i]])
    cat("\n")
  }
  
}
