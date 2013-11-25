#' Fit PLC curves
#' @description Description goes here.
#' @param dfr A dataframe that contains water potential and plc data.
#' @param varnames A vector specifying the names of the PLC and water potential data (WP) in the dataframe.
#' @export
#' @examples
#' \dontrun{
#' dfr <- read.csv("test/stemvul-ros.csv")
#' dfr <- subset(dfr, Species =="EuTe")
#' 
#' pfit <- fitplc(dfr)
#' pfit
#' plot(pfit)
#' }
fitplc <- function(dfr, varnames = c(PLC="PLC", WP="MPa"), model="Weibull", bootci=FALSE){

                   
    # Get variables out of dataframe
    Y <- dfr[[varnames["PLC"]]]
    X <- dfr[[varnames["WP"]]]
    
    # check for NA
    if(any(is.na(c(Y,X))))stop("Remove missing values first.")
    
    # Need absolute values of water potential
    if(mean(X) < 0)X <- -X
    
    
    # Weibull vulnerability curve, as re-parameterized by Ogle et al.
    fweibull <- function(P, SX, PX, X=50){
      
      V <- (X-100)*log(1-X/100)
      p <- (P/PX)^((PX*SX)/V)
      relk <- (1-X/100)^p
      
      return(relk)
    }
    
    # Calculate relative conductivity:
    relK <- (100 - Y)/100
    
    # fit
    message("Fitting nls ...", appendLF=FALSE)
    nlsfit <- nls(relK ~ fweibull(X, S, P50),
                               start=list(S=5, P50=2.5))
    message("done.")
    
    # bootstrap
    if(bootci){
      message("Fitting to bootstrap replicates ...", appendLF=FALSE)
      p <- predict_nls(nlsfit, interval="confidence")
      message("done.")
    } else {
      p <- predict_nls(nlsfit, interval="none")
    }
    
    # ci on pars.
    cipars <- suppressMessages(confint(nlsfit))
    
    l <- list()
    l$fit <- nlsfit
    l$pred <- p
    l$ci <- cipars
    l$data <- data.frame(X=X, Y=Y, relK=relK)
    class(l) <- "plcfit"
    
return(l)
}    
    

plot.plcfit <- function(x, ...){
  
    with(x, {
    plot(data$X, data$relK, ylim=c(0,1), pch=19,
         xlab=expression(Water~potential~~(-MPa)),
         ylab="Relative conductivity (0 - 1)")
    })
    if("lwr" %in% names(x$p)){
      with(x$p,{
        lines(x, lwr, type='l', lty=5)
        lines(x, upr, type='l', lty=5)
        lines(x, pred, type='l', lty=1)
      })
    } else {
      with(x$p,{
        lines(x, pred, type='l', lty=1)
      })
      
      
    }
    p50 <- coef(x$fit)["P50"]
    abline(v=p50, col="red")
    
    # method 1 : bootstrap
    #p50_ci <- quantile(x$p$boot[,2], c(0.025, 0.975))
    # method 2: confint (parametric)
    p50_ci <- x$ci[2,]
    
    abline(v=p50_ci, col="red", lty=5)
    mtext(side=3, at=p50, text=expression(P[50]), line=1, col="red")

}


print.plcfit <- function(x,...){
  cat("Class of object 'plcfit' as returned by 'fitplc'.\n")
  cat("-------------------------------------------------\n\n")
  cat("Parameters and 95% confidence interval:\n")
  cat()
  
  Estimate <- coef(x$fit)
  Table <- cbind(Estimate, x$ci)
  print(Table)
}




