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
#' pfit <- fitplc(dfr, bootci=TRUE)
#' pfit
#' plot(pfit)
#' 
#' }
#' 
#' 
#' 
#' 
#' 
fitplc <- function(dfr, varnames = c(PLC="PLC", WP="MPa"), model="Weibull", bootci=FALSE){

                   
    # Get variables out of dataframe
    Y <- dfr[[varnames["PLC"]]]
    X <- dfr[[varnames["WP"]]]
    
    # check for NA
    if(any(is.na(c(Y,X))))stop("Remove missing values first.")
    
    # Need absolute values of water potential
    if(mean(X) < 0)X <- -X
    
    # Calculate relative conductivity:
    relK <- (100 - Y)/100
    
    Data <- data.frame(X=X, relK=relK)
    
    # fit
    message("Fitting nls ...", appendLF=FALSE)
    nlsfit <- nls(relK ~ fweibull(X, S, P50),
                  data=Data,start=list(S=5, P50=2.5))
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
    
#'@export
plot.plcfit <- function(x, xlab=NULL, ylab=NULL, ylim=NULL, pch=19, 
                        plotP50=TRUE, plotci=TRUE, plotdata=TRUE, add=FALSE,
                        linecol="black", ...){
  
    if(is.null(xlab))xlab <- expression(Water~potential~~(-MPa))
    if(is.null(ylab))ylab <- "Relative conductivity (0 - 1)"
    if(is.null(ylim))ylim <- c(0,1)
    type <- ifelse(plotdata, 'p', 'n')
    
    if(!add){
      with(x, {
      plot(data$X, data$relK, ylim=ylim, pch=pch,
           xlab=xlab,
           type=type,
           ylab=ylab, ...)
      })
    } else {
      with(x, {
        points(data$X, data$relK, pch=pch, type=type,...)
      })
    }
    if("lwr" %in% names(x$p) && plotci){
      with(x$p,{
        lines(x, lwr, type='l', lty=5, col=linecol)
        lines(x, upr, type='l', lty=5, col=linecol)
      })      
    }
    with(x$p,{
      lines(x, pred, type='l', lty=1, col=linecol)
    })
    
    if(plotP50){
      p50 <- coef(x$fit)["P50"]
      abline(v=p50, col="red")
      
      # method 1 : bootstrap
      #p50_ci <- quantile(x$p$boot[,2], c(0.025, 0.975))
      # method 2: confint (parametric)
      p50_ci <- x$ci[2,]
      
      abline(v=p50_ci, col="red", lty=5)
      mtext(side=3, at=p50, text=expression(P[50]), line=1, col="red")
    }
    
}

#'@export
print.plcfit <- function(x,...){
  cat("Class of object 'plcfit' as returned by 'fitplc'.\n")
  cat("-------------------------------------------------\n\n")
  cat("Parameters and 95% confidence interval:\n\n")
  cat()
  
  print(coef(x))
  cat("\n")
}

#'@export
coef.plcfit <- function(object,...){
  
  Estimate <- coef(object$fit)
  Table <- cbind(Estimate, object$ci)
  
return(Table)
}





