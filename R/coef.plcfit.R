#'@export
coef.plcfit <- function(object, which=c("parametric","bootstrap"), ...){
  
  which <- match.arg(which)
  
  if(object$fitran){
    Estimate <- object$cinlme[,2]
    SE <- summary(object$nlmefit)$tTable[,2]
    Table <- cbind(Estimate, SE, object$cinlme[,c(1,3)])
    colnames(Table) <- c("Estimate","Std. Error","2.5%","97.5%")
  } else {
    if(which == "parametric"){
      Estimate <- summary(object$fit)$coefficients[,1:2]
      Table <- cbind(Estimate, object$ci)
    } else {
      if(!object$bootci)stop("First refit model with bootci=TRUE")
      Table <- object$bootpars
    }
  }
  
  return(Table)
}
