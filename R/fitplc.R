#' Fit PLC curves
#' @description This function fits 'percent loss conductivity' (PLC) curves. At the moment, it 
#' only fits the Weibull curve, as reparameterized by Ogle et al. (2009), and only returns
#' P50 (although soon it will return P88 or whatever values).
#' The function \code{fitplcs} can be used for batch fitting. See examples below for usage of \code{fitplc}
#' or \code{fitplcs}.
#' @param dfr A dataframe that contains water potential and plc data.
#' @param varnames A vector specifying the names of the PLC and water potential data (WP).
#' @param weights A variable used as weights in weighted non-linear regression that must be present in the dataframe (unquoted, see examples).
#' @param random Variable that specified random effects (unquoted; must be present in dfr).
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
#' @details If a variable with the name Weights is present in the dataframe, 
#' this variable will be used as the \code{weights} argument in \code{\link{nls}} to perform 
#' weighted non-linear regression. See the final example on how to use this.
#' 
#' If the \code{random} argument specifies a factor variable present in the dataframe, random effects will 
#' be estimated both for SX and PX. This affects \code{coef} as well as the confidence intervals for the fixed effects.
#'
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
#' pfit <- fitplc(dfr_eute, varnames=c(PLC="PLC", WP="MPa"), weights=Weights)
#' coef(pfit)
#' 
#' }
#' 
fitplc <- function(dfr, varnames = c(PLC="PLC", WP="MPa"),
                   weights=NULL,
                   random=NULL,
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
    
    if(!is.null(substitute(random))){
      G <- eval(substitute(random), dfr)
      fitran <- TRUE
      if(bootci){
        bootci <- FALSE
        message("Not performing bootstrap when random effects present.")
      }
    } else {
      fitran <- FALSE
    }
    
    W <- eval(substitute(weights), dfr)
    
    # check for NA
    if(any(is.na(c(Y,P))))stop("Remove missing values first.")
    
    # Need absolute values of water potential
    if(mean(P) < 0)P <- -P
    
    # Calculate relative conductivity:
    relK <- (100 - Y)/100
    
    if(!fitran){
      Data <- data.frame(P=P, relK=relK)
    } else {
      Data <- data.frame(P=P, relK=relK, G=G)
    }
    
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

    # Weighted NLS
    if(!is.null(W)){
        nlsfit <- nls(relK ~ fweibull(P,SX,PX,X),
                    data=Data, start=list(SX=Sh, PX=pxstart),
                    weights=W)
        if(fitran){
          nlmefit <- nlme(relK ~ fweibull(P, SX, PX, X),
                     fixed=list(SX ~ 1, PX ~ 1),
                     random= SX + PX ~ 1|G,
                     start=list(fixed=c(SX=coef(nlsfit)["SX"], 
                                        PX=coef(nlsfit)["PX"])),
                     weights=W,
                     data=Data)
        }
    } else {
      
    # Ordinary NLS
      nlsfit <- nls(relK ~ fweibull(P,SX,PX,X),
                    data=Data, start=list(SX=Sh, PX=pxstart))
      if(fitran){
        nlmefit <- nlme(relK ~ fweibull(P, SX, PX, X),
                        fixed=list(SX ~ 1, PX ~ 1),
                        random= SX + PX ~ 1|G,
                        start=list(fixed=c(SX=coef(nlsfit)["SX"], 
                                           PX=coef(nlsfit)["PX"])),
                        data=Data)
      }
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
    
    # Predictions at innermost random effect
    if(fitran){
      
      d <- split(Data, Data$G)
      pm <- list()
      for(i in 1:length(d)){
        ps <- seq(min(d[[i]]$P),max(d[[i]]$P),length=101)
        newdat <- data.frame(P=ps, 
                             G=unique(d[[i]]$G), X=x)
        y <- predict(nlmefit, newdat)
        pm[[i]] <- data.frame(x=ps, y=y) 
      }
      ps <- seq(min(P),max(P),length=101)
      newdat <- data.frame(P=ps, X=x)
      pmf <- data.frame(x=ps, y=predict(nlmefit, newdat, level=0))
      
    } else {
      pm <- NA
      pmf <- NA
    }
    
    # ci on pars.
    cipars <- suppressMessages(confint(nlsfit))
    
    if(bootci){
      cisx <- quantile(p$boot[,"SX"], c(0.025,0.975))
      cipx <- quantile(p$boot[,"PX"], c(0.025,0.975))

      bootpars <- matrix(c(coef(nlsfit),cisx[1],cipx[1],cisx[2],cipx[2]), nrow=2,
                         dimnames=list(c("SX","PX"),c("Estimate","2.5%","97.5%")))
    } else {
      bootpars <- NA
    }               
    
    l <- list()
    l$fit <- nlsfit
    l$pred <- p
    l$prednlme <- pm
    l$prednlmefix <- pmf
    l$ci <- cipars
    l$bootpars <- bootpars
    l$data <- data.frame(P=P, Y=Y, relK=relK)
    l$x <- x
    l$fitran <- fitran
    if(fitran){
      l$nlmefit <- nlmefit
      l$cinlme <- intervals(nlmefit,which="fixed")$fixed
      l$ranvar <- substitute(random)
    } else {
      l$nlmefit <- NA
      l$cinlme <- NA
      l$ranvar <- NA
    }
      
    l$bootci <- bootci
    
    class(l) <- "plcfit"
    
return(l)
}    



#'@rdname fitplc
#'@export
plot.plcfit <- function(x, xlab=NULL, ylab=NULL, ylim=NULL, pch=19, 
                        plotPx=TRUE, plotci=TRUE, plotdata=TRUE, add=FALSE,
                        selines=c("parametric","bootstrap"),
                        plotrandom=FALSE,
                        linecol="black", what=c("relk","embol"), ...){
  
  
    selines <- match.arg(selines)
  
    if(is.null(xlab))xlab <- expression(Water~potential~~(-MPa))
      
    type <- ifelse(plotdata, 'p', 'n')
    what <- match.arg(what)
    
    if(plotrandom && !x$fitran)
      stop("To plot random effects predictions, refit with 'random' argument.")
    
    if(what == "relk"){
      if(is.null(ylab))ylab <- "Relative conductivity (0 - 1)"
      x$data$Y <- x$data$relK
      if(is.null(ylim))ylim <- c(0,1)
    }
    toEmbol <- function(k)100 - 100*k
    if(what == "embol"){
      if(is.null(ylab))ylab <- "% Embolism"
      
      x$data$Y <- toEmbol(x$data$relK)
      if(x$bootci){
        x$pred$lwr <- toEmbol(x$pred$lwr)
        x$pred$upr <- toEmbol(x$pred$upr)
      }
      x$p$pred <- toEmbol(x$pred$pred)
      if(is.null(ylim))ylim <- c(0,100)
      
      if(x$fitran && plotrandom){
        
        ng <- length(x$prednlme)
        for(i in 1:ng){
          x$prednlme[[i]]$y <- toEmbol(x$prednlme[[i]]$y)
        }
        x$prednlmefix$y <- toEmbol(x$prednlmefix$y)
      }
      
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
    if(!plotrandom){
      if(x$bootci && plotci){
        with(x$pred,{
          lines(x, lwr, type='l', lty=5, col=linecol)
          lines(x, upr, type='l', lty=5, col=linecol)
        })      
      }
      with(x$pred,{
        lines(x, pred, type='l', lty=1, col=linecol)
      })
    }
    
    if(plotrandom){
      for(i in 1:length(x$prednlme)){
        with(x$prednlme[[i]], lines(x,y,type='l'))
      }  
      with(x$prednlmefix, lines(x,y,type='l',lwd=2, col="blue"))
    }
    
    if(plotPx){
      if(!x$fitran){
        px <- coef(x$fit)["PX"]
        
        if(selines == "bootstrap"){
          px_ci <- x$bootpars[2,2:3]
        } else {
          px_ci <- x$ci[2,]
        }
      } else {
        px <- fixef(x$nlmefit)["PX"]
        px_ci <- x$cinlme[2,]
      }
      
      abline(v=px, col="red")
      abline(v=px_ci, col="red", lty=5)
      mtext(side=3, at=px, text=expression(Px), line=0, col="red", cex=0.7)
    }
    
}

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






