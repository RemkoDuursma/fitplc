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
#' @param citype Either 'polygon' (default), or 'lines', specifying formatting of the confidence interval in the plot.
#' @param linecol The color of the fitted curve (or color of the random effects curves if plotrandom=TRUE).
#' @param linecol2 The color of the fixed effects curve (if plotrandom=TRUE; otherwise ignored).
#' @param pxlinecol The color of the lines indicating Px and its confidence interval 
#' @param pxcex Character size for the Px label above the Y-axis.
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
#' # Here we also set the starting values (which is sometimes needed).
#' allfit <- fitplcs(dfr, "Species", varnames=c(PLC="PLC", WP="MPa"),
#' startvalues=list(Px=4, Sx=10))
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
#' @importFrom nlme fixef
#' @importFrom nlme nlme
#' @importFrom nlme intervals
fitplc <- function(dfr, varnames = c(PLC="PLC", WP="MPa"),
                   weights=NULL,
                   random=NULL,
                   model="Weibull", 
                   startvalues=list(Px=3, S=20), x=50,
                   bootci=TRUE, ...){

                   
  
  # Find out if called from fitcond.
  mc <- names(as.list(match.call()))
  
  condfit <- "calledfromfitcond" %in% mc
  
  # Get Kmax value
  if(!"Kmax" %in% mc)Kmax <- 1
  
  
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
    cipars <- try(suppressMessages(confint(nlsfit)), silent=TRUE)
    if(inherits(cipars, "try-error")){
      cipars <- matrix(rep(NA,4),ncol=2)
      dimnames(cipars) <- list(c("SX","PX"), c("2.5%","97.5%")) 
    }
    
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
    l$condfit <- condfit
    l$Kmax <- Kmax
    
    class(l) <- "plcfit"
    
return(l)
}    


