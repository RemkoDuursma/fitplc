#' Fit PLC curves
#' @description This function fits 'percent loss conductivity' (PLC) curves. At the moment, it 
#' only fits the Weibull curve, as reparameterized by Ogle et al. (2009), and only returns
#' P50 (although soon it will return P88 or whatever values).
#' The function \code{fitplcs} can be used for batch fitting. See examples below for usage of \code{fitplc}
#' or \code{fitplcs}.
#' @param dfr A dataframe that contains water potential and plc data.
#' @param varnames A vector specifying the names of the PLC, water potential data (WP), and optional Weights in the dataframe.
#' @param model At the moment, only 'Weibull' is allowed.
#' @param startvalues A list of starting values. If set to NULL, \code{fitplc} will attempt to guess starting values.
#' @param bootci If TRUE, also computes the bootstrap confidence interval.
#' @param x A fitted curve returned by \code{fitplc}
#' @param plotPx Logical (default TRUE), whether to plot a vertical line for the P50.
#' @param plotci Logical (default TRUE), whether to plot the confidence interval (if computed with bootci=TRUE).
#' @param plotdata Logical (default TRUE), whether to add the data to the plot.
#' @param add Logical (default FALSE), whether to add the plot to a current device. This is useful to overlay two plots or curves, for example.
#' @param linecol the color of the line
#' @param what Either 'relk' or 'embol'; it will plot either relative conductivity or percent embolism.
#' @details If a variable with the name Weights is present in the dataframe, this variable will be used as the \code{weights} argument in \code{\link{nls}} to perform weighted non-linear regression. See the final example on how to use this.
#' 
#' A plot method is available for the fitted object, see examples on how to use it.
#' @export
#' @rdname fitplc
#' @examples
#' \dontrun{
#' 
#' # First read a dataframe (in this example, from the folder 'test')
#' dfr <- read.csv("test/stemvul-ros.csv")
#' 
#' # 1. Fit one species (or fit all, see next example)
#' dfr_eute <- subset(dfr, Species =="EuTe")
#' 
#' # Make fit. Store results in object 'pfit'
#' # 'varnames' specifies the names of the 'PLC' variable in the dataframe,
#' # and water potential (WP). 
#' pfit <- fitplc(dfr_eute, varnames=c(PLC="PLC", WP="MPa"))
#' 
#' # Look at fit
#' pfit
#' 
#' # Make a standard plot. The default plot is 'relative conductivity',
#' # (which is 1.0 where PLC = 0).
#' plot(pfit)
#' 
#' # Or plot the percent embolism
#' plot(pfit, what="embol")
#' 
#' # Get the coefficients of the fit.
#' coef(pfit)
#' 
#' # 2. Fit all species in the dataset.
#' allfit <- fitplcs(dfr, "Species", varnames=c(PLC="PLC", WP="MPa"))
#' 
#' # Make three plots
#' # windows(10,8) # optional : open up a window and split.
#' # par(mfrow=c(3,1), mar=c(4,4,2,2))
#' for(i in 1:3)plot(allfit[[i]], xlim=c(0,7), main=names(allfit)[i])
#' 
#' # Coefficients show the estimates and 95% CI (given by 'lower' and 'upper')
#' # Based on the CI's, species differences can be decided.
#' coef(allfit)
#' 
#' # 3. Specify Weights. The default variable name is Weights, if present in the dataset
#' # it will be used for weighted non-linear regression
#' dfr_eute$Weights <- abs(50-dfr_eute$PLC)^1.2
#' pfit <- fitplc(dfr_eute, varnames=c(PLC="PLC", WP="MPa"))
#' coef(pfit)
#' 
#' # NOTE: to turn weights off, remove from the dataframe
#' dfr_eute$Weights <- NULL
#' pfit <- fitplc(dfr_eute, varnames=c(PLC="PLC", WP="MPa"))
#' coef(pfit)
#' }
#' 
fitplc <- function(dfr, varnames = c(PLC="PLC", WP="MPa"),
                   weights=NULL,
                   model="Weibull", 
                   startvalues=list(Px=3, S=20), x=50,
                   Weights=NULL,
                   bootci=TRUE){

                   
    # Get variables out of dataframe
    if(!varnames["PLC"] %in% names(dfr))
      stop("Check variable name for PLC!")
    if(!varnames["WP"] %in% names(dfr))
      stop("Check variable name for water potential!")
    
    Y <- dfr[[varnames["PLC"]]]
    P <- dfr[[varnames["WP"]]]
    
    W <- eval(substitute(weights), dfr)
    
    # check for NA
    if(any(is.na(c(Y,P))))stop("Remove missing values first.")
    
    # Need absolute values of water potential
    if(mean(P) < 0)P <- -P
    
    # Calculate relative conductivity:
    relK <- (100 - Y)/100
    
    Data <- data.frame(P=P, relK=relK)
    
    # guess starting values
    if(is.null(startvalues)){
      pxstart <- (1-x/100)*(max(P) - min(P))
      Sh <- 15
    } else {
      pxstart <- startvalues$Px
      Sh <- startvalues$S
    }
    
    # fit
    Data$X <- x
    message("Fitting nls ...", appendLF=FALSE)

    if(!is.null(W)){
      nlsfit <- nls(relK ~ fweibull(P,SX,PX,X),
                    data=Data, start=list(SX=Sh, PX=pxstart),
                    weights=W)
    } else {
      nlsfit <- nls(relK ~ fweibull(P,SX,PX,X),
                    data=Data, start=list(SX=Sh, PX=pxstart))
    }
    message("done.")
    
    # bootstrap
    if(bootci){
      message("Fitting to bootstrap replicates ...", appendLF=FALSE)
      p <- predict_nls(nlsfit, xvarname="P", interval="confidence", data=Data, 
                       startList=list(SX=Sh, PX=pxstart), weights=W)
      message("done.")
    } else {
      p <- predict_nls(nlsfit, xvarname="P", interval="none", data=Data, 
                       startList=list(SX=Sh, PX=pxstart), weights=W)
    }
    
    # ci on pars.
    cipars <- suppressMessages(confint(nlsfit))
    
    l <- list()
    l$fit <- nlsfit
    l$pred <- p
    l$ci <- cipars
    l$data <- data.frame(P=P, Y=Y, relK=relK)
    l$x <- x
    l$bootci <- bootci
    
    class(l) <- "plcfit"
    
return(l)
}    



#'@rdname fitplc
#'@export
plot.plcfit <- function(x, xlab=NULL, ylab=NULL, ylim=NULL, pch=19, 
                        plotPx=TRUE, plotci=TRUE, plotdata=TRUE, add=FALSE,
                        linecol="black", what=c("relk","embol"), ...){
  
    if(is.null(xlab))xlab <- expression(Water~potential~~(-MPa))
      
    type <- ifelse(plotdata, 'p', 'n')
    what <- match.arg(what)
    
    if(what == "relk"){
      if(is.null(ylab))ylab <- "Relative conductivity (0 - 1)"
      x$data$Y <- x$data$relK
      if(is.null(ylim))ylim <- c(0,1)
    }
    if(what == "embol"){
      if(is.null(ylab))ylab <- "% Embolism"
      
      x$data$Y <- 100 - 100*x$data$relK
      if(x$bootci){
        x$p$lwr <- 100 - 100*x$p$lwr
        x$p$upr <- 100 - 100*x$p$upr
      }
      x$p$pred <- 100 - 100*x$p$pred
      if(is.null(ylim))ylim <- c(0,100)
    }
    
    if(!add){
      with(x, {
      plot(data$P, data$Y, ylim=ylim, pch=pch,
           xlab=xlab,
           type=type,
           ylab=ylab, ...)
      })
    } else {
      with(x, {
        points(data$P, data$Y, pch=pch, type=type,...)
      })
    }
    if(x$bootci && plotci){
      with(x$p,{
        lines(x, lwr, type='l', lty=5, col=linecol)
        lines(x, upr, type='l', lty=5, col=linecol)
      })      
    }
    with(x$p,{
      lines(x, pred, type='l', lty=1, col=linecol)
    })
    
    if(plotPx){
      px <- coef(x$fit)["PX"]
      abline(v=px, col="red")
      
      # method 1 : bootstrap
      #p50_ci <- quantile(x$p$boot[,2], c(0.025, 0.975))
      # method 2: confint (parametric)
      px_ci <- x$ci[2,]
      
      abline(v=px_ci, col="red", lty=5)
      mtext(side=3, at=px, text=expression(Px), line=0, col="red", cex=0.7)
    }
    
}

#'@export
print.plcfit <- function(x,...){
  cat("Class of object 'plcfit' as returned by 'fitplc'.\n")
  cat("-------------------------------------------------\n\n")
  cat("Parameters, SE, and 95% confidence interval:\n\n")
  cat()
  
  print(coef(x))
  cat("\n")
}

#'@export
coef.plcfit <- function(object,...){
  
  Estimate <- summary(object$fit)$coefficients[,1:2]
  Table <- cbind(Estimate, object$ci)
  
return(Table)
}

#' @export
#' @rdname fitplc
#' @param group Character; variable in the dataframe that specifies groups. The curve will be fit for every group level.
#' @param \dots Further parameters passed to \code{fitplc}.
fitplcs <- function(dfr, group, ...){
  
  if(!group %in% names(dfr))
    stop("You must provide a name in the dataframe to fit by.")
  
  dfrs <- split(dfr, dfr[,group])
  
  fits <- lapply(dfrs, function(x)fitplc(x, ...))
  class(fits) <- "manyplcfit"

return(fits)
}

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


#'@export
coef.manyplcfit <- function(object, ...){
  
  x <- lapply(object,coef)
  vf <- function(m)as.vector(t(m))
  dfr <-  as.data.frame(do.call(rbind,lapply(x,vf)))
  names(dfr) <- c("S","S_SE","S_lower","S_upper",
                  "Px","Px_SE","P50_lower","P50_upper")
  
  return(dfr)
}






