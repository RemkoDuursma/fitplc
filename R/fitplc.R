#' Fit a PLC curve
#' @description This function fits the Weibull curve to measurements of stem or leaf conductivity 
#' measurements at various water potentials. If measurements are organized as 'percent loss conductivity' (PLC), use the \code{fitplc} function. If they are organized as the actual conductance or conductivity (as is common for leaf hydraulic  conductance data, for example), use the \code{fitcond} function. See Details and Examples for more information on how to use these functions. 
#' 
#' It is also possible to fit multiple curves at once, for example one for each species or site, 
#' with the \code{fitplcs} and \code{fitconds} functions.
#' 
#' See \code{\link{plot.plcfit}} for documentation on plotting methods for the fitted objects, and the examples below.
#' @param dfr A dataframe that contains water potential and plc or conductivity/conductance data.
#' @param varnames A vector specifying the names of the PLC and water potential data (WP).
#' @param weights A variable used as weights in weighted non-linear regression that must be present in the dataframe (unquoted, see examples).
#' @param random Variable that specified random effects (unquoted; must be present in dfr).
#' @param x If the P50 is to be returned, x = 50. Set this value if other points of the PLC curve should be estimated (although probably more robustly done via \code{\link{getPx}}).
#' @param model At the moment, only 'Weibull' is allowed.
#' @param startvalues A list of starting values. If set to NULL, \code{fitplc} will attempt to guess starting values.
#' @param bootci If TRUE, also computes the bootstrap confidence interval.
#' @param nboot The number of bootstrap replicates (only relevant when \code{bootci=TRUE}).
#' @param quiet Logical (default FALSE), if TRUE, don't print any messages.
#' @param Kmax Maximum conduct(ance)(ivity), optional (and only when using \code{fitcond}). See Examples.
#' @param WP_Kmax Water potential above which Kmax will be calculated from the data. Optional (and only when using \code{fitcond}). See Examples.
#' @details If a variable with the name Weights is present in the dataframe, 
#' this variable will be used as the \code{weights} argument in \code{\link{nls}} to perform 
#' weighted non-linear regression. See the final example on how to use this.
#' 
#' If the \code{random} argument specifies a factor variable present in the dataframe, random effects will 
#' be estimated both for SX and PX. This affects \code{coef} as well as the confidence intervals for the fixed effects.
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
#' # In this example, we use only 50 bootstrap replicates but recommend you set this
#' # to 1000 or so.
#' pfit <- fitplc(dfr1, varnames=c(PLC="PLC", WP="MPa"), nboot=50)
#' 
#' # Look at fit
#' pfit
#' 
#' # Make a standard plot. The default plot is 'relative conductivity',
#' # (which is 1.0 where PLC = 0). For plotting options, see ?plot.plcfit
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
#' # In this example, we use only 50 bootstrap replicates but recommend you set this
#' # to 1000 or so. 
#' allfit <- fitplcs(stemvul, "Species", varnames=c(PLC="PLC", WP="MPa"), nboot=50)
#' 
#' # 3. Plot the fits.
#' plot(allfit, onepanel=TRUE, plotci=FALSE, selines="none", pxlinecol="dimgrey")
#'
#' # Coefficients show the estimates and 95% CI (given by 'lower' and 'upper')
#' # Based on the CI's, species differences can be decided.
#' coef(allfit)
#' 
#' # 3. Specify Weights. The default variable name is Weights, if present in the dataset
#' # it will be used for weighted non-linear regression
#' # In this example, we use only 50 bootstrap replicates but recommend you set this
#' # to 1000 or so. 
#' dfr1$Weights <- abs(50-dfr1$PLC)^1.2
#' pfit <- fitplc(dfr1, varnames=c(PLC="PLC", WP="MPa"), weights=Weights, nboot=50)
#' coef(pfit)
#' 
#' # 4. Fit the Weibull curve directly to the raw conductance data. 
#' # Use this option when you don't want to transform your data to PLC. 
#' # You have two options: specify the 'maximum' conductance yourself (and provide Kmax), 
#' # or set the threshold water potential (Kmax_WP), which is then used to calculate Kmax
#' # (from the average of the conductance values where WP > Kmax_WP).
#' 
#' # Option 1 : maximum conductivity (i.e. at full hydration) is known, and used as input.
#' kfit1 <- fitcond(dfr1, varnames=c(K="Cond", WP="MPa"), Kmax=7.2, nboot=50)
#'
#' # Option 2 : calculate maximum cond. from data where water potential : -0.3 MPa.
#' # In this example, we use only 50 bootstrap replicates but recommend you set this
#' # to 1000 or so. 
#' kfit2 <- fitcond(dfr1, varnames=c(K="Cond", WP="MPa"), WP_Kmax = -0.3, nboot=50)
#' # Use plot(kfit1) as for fitplc, as well as coef() etc.
#' 
#' # Fit multiple conductivity curves at once (bootstrap omitted for speed).
#' kfits3 <- fitconds(stemvul, "Species", varnames=list(K="Cond", WP="MPa"), WP_Kmax=-0.3, boot=FALSE)
#' plot(kfits3, onepanel=TRUE, ylim=c(0,12), selines="none")
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
                   model=c("Weibull","sigmoidal"), 
                   #startvalues=list(Px=3, S=20), 
                   x=50,
                   coverage=0.95,
                   bootci=TRUE,
                   nboot=999,
                   quiet=FALSE,
                   ...){

    # sigmoid_lm
    # sigmoid_lme
    # weibull_nls
    # weibull_nlme
    
    # return list with
    # fit
    # cipars (estimates, CI for Sx, Px with norm and/or boot)
    # predictions fixed effect, confidence interval (lwr,upr)
    # predictions innermost random effect
  
    model <- match.arg(model)
  
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
    

    if(!is.null(substitute(random))){
      G <- eval(substitute(random), dfr)
      fitran <- TRUE
      if(bootci){
        bootci <- FALSE
        if(!quiet)message("Not performing bootstrap when random effects present.")
      }
    } else {
      G <- NA
      fitran <- FALSE
    }
    

    # Extract data
    plc <- dfr[[varnames["PLC"]]]
    P <- dfr[[varnames["WP"]]]
    if(any(is.na(c(plc,P))))stop("Remove missing values first.")
    relK <- plc_to_relk(plc)
    
    # Need absolute values of water potential
    if(mean(P) < 0)P <- -P
    
    # weights, if provided
    W <- eval(substitute(weights), dfr)
    
    # Dataset tidied
    Data <- data.frame(P=P, PLC=plc, relK=relK, G=G)
    Data$minP <- -Data$P  # negative valued water potential
    
    
    
    if(model == "sigmoidal"){
      
      if(!fitran){
        f <- do_sigmoid_fit(Data, boot=TRUE, nboot=nboot)
        
        cf <- sigfit_coefs(f$boot[,1], f$boot[,2],x=x)
        boot_Sx <- cf$Sx  
        boot_Px <- cf$Px
        
        p <- coef(f$fit)
        mf <- sigfit_coefs(p[1],p[2],x=x)
        ml_Sx <- mf$Sx
        ml_Px <- mf$Px
          
        # Coefficients matrix
        cipars <- rbind(c(ml_Sx, boot_ci(boot_Sx, coverage)),
                        c(ml_Px, boot_ci(boot_Px, coverage)))
        
        dimnames(cipars) <- list(c("SX","PX"), 
                                 c("Estimate", "Boot - 2.5%","Boot - 97.5%"))
        
        # f must be component with 'fit' and 'boot'
        pred <- get_boot_pred_sigmoid(f, Data, coverage)
        
        
      
      } else {
        
        
        #--> to subfunction
        # This is necessary - might have to revisit this method.
        #Data$PLCf <- pmax(0.1, pmin(99.9, Data$PLC))
        fit <- do_sigmoid_lme_fit(Data, W)
        
        Px_ci <- deltaMethod(fit, "b0/b1", parameterNames=c("b0","b1"))
        Sx_ci <- deltaMethod(fit, "100*b1/4", parameterNames=c("b0","b1"))
        cipars <- rbind(Sx_ci, Px_ci)
        cipars$SE <- NULL
        names(cipars)[2:3] <- c("Norm - 2.5%","Norm - 97.5%")
        rownames(cipars) <- c("SX","PX")
        
        #--> to here
        
        predran <- lapply(split(Data, Data$G), function(x){
          
          ps <- seq_within(x$minP, n=101)
          newdat <- data.frame(minP=ps, G=unique(x$G))
          
          list(x=-ps, fit=sigmoid_untrans(unname(predict(fit, newdat)))) 
        })
        
        ps <- seq_within(Data$minP, n=101)
        newdat <- data.frame(minP=ps, X=x)
        pred <- list(x=-ps, fit=predict(fit, newdat, level=0), ran=predran)
        pred$fit <- sigmoid_untrans(pred$fit)
        
      }
      
      l <- list(fit=fit, pred=pred, cipars=cipars, data=Data, x=x, Kmax=1)
      
      # l$lmefit <- fit
      # l$b <- f$boot
      # l$model <- model
      # l$data <- Data[,c("P","PLC","relK")]
      # l$cipars <- cipars
      # 
      # # boot CI - used in plotting
      # l$pred <- pred
      # 
      # l$condfit <- condfit
      # l$fitran <- fitran
      # l$bootci <- TRUE
      
      
    }

    
    # move to subfunction
    if(model == "Weibull"){
    
      # guess starting values from sigmoidal
      f <- do_sigmoid_fit(Data, boot=FALSE)
      p <- coef(f$fit)
      sp <- sigfit_coefs(p[1],p[2],x=x)
      
      # fit
      Data$X <- x
      if(!quiet)message("Fitting nls ...", appendLF=FALSE)
  
      # Weighted NLS
      if(!is.null(W)){
          nlsfit <- nls(relK ~ fweibull(P,SX,PX,X),
                      data=Data, start=list(SX=sp$Sx, PX=sp$Px),
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
                      data=Data, start=list(SX=sp$Sx, PX=sp$Px))
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
                         startList=list(SX=sp$Sx, PX=sp$Px), weights=W, nboot=nboot)
        if(!quiet)message("done.")
      } else {
        p <- predict_nls(nlsfit, xvarname="P", interval="none", data=Data, 
                         startList=list(SX=sp$Sx, PX=sp$Px), weights=W, nboot=nboot)
      }
      
      # Predictions at innermost random effect
      if(fitran){
        
        # pm: innermost random effects (list)
        # pmf : fixed effects predictions
        
      } else {
        pm <- NA
        pmf <- NA
      }
      
      # ci on pars.
      cipars <- try(suppressMessages(confint(nlsfit)), silent=TRUE)
      if(inherits(cipars, "try-error")){
        cipars <- matrix(rep(NA,4),ncol=2)
      }
      cipars <- cbind(coef(nlsfit), cipars)
      dimnames(cipars) <- list(c("SX","PX"), c("Estimate", "Norm - 2.5%","Norm - 97.5%"))
      
      if(bootci){
        cisx <- quantile(p$boot[,"SX"], c(0.025,0.975))
        cipx <- quantile(p$boot[,"PX"], c(0.025,0.975))
  
        bootpars <- matrix(c(cisx[1],cipx[1],cisx[2],cipx[2]), nrow=2,
                           dimnames=list(c("SX","PX"),c("Boot - 2.5%","Boot - 97.5%")))
        cipars <- cbind(cipars, bootpars)
      } else {
        bootpars <- NA
      }               
    
      l <- list()
      l$fit <- nlsfit
      l$pred <- p
      l$prednlme <- pm
      l$prednlmefix <- pmf
      l$cipars <- cipars
      l$data <- data.frame(P=P, PLC=plc, relK=relK)
      l$x <- x
      l$model <- model
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
      l$nboot <- nboot
      

    }
    
    class(l) <- "plcfit"
    l$condfit <- condfit
    l$fitran <- fitran
    l$bootci <- bootci
    return(l)
}    



do_sigmoid_fit <- function(data, W=NULL, boot=FALSE, nboot){
  
  # This is necessary - might have to revisit this method.
  #data$PLCf <- pmax(0.1, pmin(99.9, data$PLC))
  
  data <- subset(data, PLC < 100 & PLC > 0)

  # Transformation as per P&vW
  data$logPLC <- log(100/data$PLC - 1)
  
  if(!is.null(W)){
    lmfit <- lm(logPLC ~ minP, data=data, weights=W)
    br <- if(boot) suppressWarnings(bootfit(lmfit, n=nboot, Data=data, startList=NULL, weights=W)) else NA
  } else {
    lmfit <- lm(logPLC ~ minP, data=data)
    br <- if(boot) suppressWarnings(bootfit(lmfit, n=nboot, Data=data, startList=NULL)) else NA
  }
  
  return(list(fit=lmfit, boot=br))
}

do_sigmoid_lme_fit <- function(data, W=NULL){
  
  data <- subset(data, PLC > 0 & PLC < 100)
  
  data$logPLC <- log(100/data$PLC - 1)
  
  fit <- lme(logPLC ~ minP,
             random= ~1|G,
             weights=W,
             data=data)
  
return(fit)
}



# Calculate Sx, Px, given log-linear fit of sigmoidal model
sigfit_coefs <- function(c1,c2,x){
  b <- c1 / c2
  Sx <- 100 * c2/4
  Px <- ab_to_px(c2, b, x)
  
  list(Px=unname(Px), Sx=unname(Sx))
}

get_boot_pred_sigmoid <- function(f, data, coverage){
  
  preddfr <- data.frame(minP=seq(min(data$minP), max(data$minP), length=101))
  
  normpred <- sigmoid_untrans(predict(f$fit, preddfr, interval="none"))
  
  bootm <- apply(f$boot,1, function(x)x[1] + x[2]*preddfr$minP)
  bootpred <- as.data.frame(t(apply(bootm, 1, boot_ci, coverage=coverage)))
  names(bootpred) <- c("lwr","upr")
  bootpred <- lapply(bootpred, sigmoid_untrans)
  bootpred$x <- -preddfr$minP
  bootpred$fit <- normpred
  
  return(bootpred)
}

