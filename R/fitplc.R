#' Fit a PLC curve
#' @description This function fits the Weibull curve to measurements of stem or leaf conductivity 
#' measurements at various water potentials. If measurements are organized as 'percent loss conductivity' (PLC), use the \code{fitplc}
#' function. If they are organized as the actual conductance or conductivity (as is common for leaf hydraulic conductance data, for example),
#' use the \code{fitcond} function. See Details and Examples for more information on how to use these functions. 
#' 
#' It is also possible to fit multiple curves at once, for example one for each species or site, 
#' with the \code{fitplcs} and \code{fitconds} functions.
#' @param dfr A dataframe that contains water potential and plc or conductivity/conductance data.
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
#' @param plotrandom If TRUE, and the model was fit with a random effect, plots the random effects predictions.
#' @param pxlinecol The color of the lines indicating Px and its confidence interval 
#' @param pxcex Character size for the Px label above the Y-axis.
#' @param what Either 'relk' or 'embol'; it will plot either relative conductivity or percent embolism.
#' @param quiet Logical (default FALSE), if TRUE, don't print any messages.
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
#' @importFrom nlme fixef
#' @importFrom nlme nlme
#' @importFrom nlme intervals

#' @rdname fitplc
#' @examples
#'
#' # We use the built-in example dataset 'stemvul' in the examples below. See ?stemvul.
#'   
#' # 1. Fit one species (or fit all, see next example)
#' dfr1 <- subset(stemvul, Species =="dpap")
#' 
#' # Make fit. Store results in object 'pfit'
#' # 'varnames' specifies the names of the 'PLC' variable in the dataframe,
#' # and water potential (WP). 
#' pfit <- fitplc(dfr1, varnames=c(PLC="PLC", WP="MPa"))
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
#' allfit <- fitplcs(stemvul, "Species", varnames=c(PLC="PLC", WP="MPa"))
#' 
#' 3. Plot the fits.
#' plot(allfit, onepanel=TRUE, plotci=FALSE, selines="none", pxlinecol="dimgrey")
#'
#' # Coefficients show the estimates and 95% CI (given by 'lower' and 'upper')
#' # Based on the CI's, species differences can be decided.
#' coef(allfit)
#' 
#' # 3. Specify Weights. The default variable name is Weights, if present in the dataset
#' # it will be used for weighted non-linear regression
#' dfr1$Weights <- abs(50-dfr1$PLC)^1.2
#' pfit <- fitplc(dfr1, varnames=c(PLC="PLC", WP="MPa"), weights=Weights)
#' coef(pfit)
#' 
#' # 4. Fit the Weibull curve directly to the raw conductance data. Use this option when you don't want to transform your data to PLC. You have two options: specify the 'maximum' conductance yourself (and provide Kmax), or set the threshold water potential (Kmax_WP), which is then used to calculate Kmax (from the average of the conductance values where WP > Kmax_WP). 
#' 
#' # Option 1 : maximum conductivity (i.e. at full hydration) is known, and used as input.
#' kfit1 <- fitcond(dfr1, varnames=c(K="Cond", WP="MPa"), Kmax=7.2)
#'
#' # Option 2 : calculate maximum cond. from data where water potential : -0.3 MPa.
#' kfit2 <- fitcond(dfr1, varnames=c(K="Cond", WP="MPa"), WP_Kmax = -0.3)
#' # Use plot(kfit1) as for fitplc, as well as coef() etc.
#' 
#' # 5. Random effects.
#' # This example takes into account the fact that the individual data points for a species are not 
#' # independent, but rather clustered by branch. 
#' fitr <- fitplc(dfr1, random=Branch)
#' 
#' # Visualize the random effects.
#' plot(fitr, plotrandom=TRUE)
fitplc <- function(dfr, varnames = c(PLC="PLC", WP="MPa"),
                   weights=NULL,
                   random=NULL,
                   model="Weibull", 
                   startvalues=list(Px=3, S=20), x=50,
                   bootci=TRUE,
                   quiet=FALSE,
                   ...){

                   
  
    # Find out if called from fitcond.
    mc <- names(as.list(match.call()))
    
    condfit <- "calledfromfitcond" %in% mc
    
    # Get Kmax value
    if(!"Kmax" %in% mc){
      Kmax <- 1
    } else {
      Kmax <- list(...)$Kmax
    }
    
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
        if(!quiet)message("Not performing bootstrap when random effects present.")
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
    if(!quiet)message("Fitting nls ...", appendLF=FALSE)

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
    if(!quiet)message("done.")
    
    # bootstrap
    if(bootci){
      if(!quiet)message("Fitting to bootstrap replicates ...", appendLF=FALSE)
      p <- predict_nls(nlsfit, xvarname="P", interval="confidence", data=Data, 
                       startList=list(SX=Sh, PX=pxstart), weights=W)
      if(!quiet)message("done.")
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


