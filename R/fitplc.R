#' Fit a PLC curve
#' @description Fit a curve to measurements of stem or leaf conductivity at various water potentials. If measurements are organized as 'percent loss conductivity' (PLC), use the \code{fitplc} function. If they are organized as the actual conductance or conductivity (as is common for leaf hydraulic  conductance data, for example), use the \code{fitcond} function. You can choose to either fit the Weibull function (the default), or the sigmoidal-exponential model. See Details and Examples for more information on how to use these functions. 
#' 
#' It is also possible to fit multiple curves at once, for example one for each species or site, with the \code{fitplcs} and \code{fitconds} functions. This is useful when you have data for multiple curves organized in one file.
#' 
#' Random effects may be incorporated via the \code{random} argument (see Examples), in which case \code{nlme} will be used (in case of the Weibull), or \code{lme} (in case of the sigmoidal model).
#'
#' See \code{\link{plot.plcfit}} for documentation on plotting methods for the fitted objects, and the examples below.
#'
#' @param dfr A dataframe that contains water potential and plc or conductivity/conductance data.
#' @param varnames A vector specifying the names of the PLC and water potential data (see Examples).
#' @param weights A variable used as weights that must be present in the dataframe (unquoted, see examples).
#' @param random Variable that specifies random effects (unquoted; must be present in dfr).
#' @param x If the P50 is to be returned, x = 50. Set this value if other points of the PLC curve should be estimated (although probably more robustly done via \code{\link{getPx}}).
#' @param coverage The coverage of the confidence interval for the parameters (0.95 is the default).
#' @param model Either 'Weibull' or 'sigmoidal'. See Details.
#' @param startvalues Obsolete - starting values for Weibull now estimated from sigmoidal model fit.
#' @param bootci If TRUE, also computes the bootstrap confidence interval.
#' @param nboot The number of bootstrap replicates used for calculating confidence intervals.
#' @param quiet Logical (default FALSE), if TRUE, don't print any messages.
#' @param Kmax Maximum conduct(ance)(ivity), optional (and only when using \code{fitcond}). See Examples.
#' @param WP_Kmax Water potential above which Kmax will be calculated from the data. Optional (and only when using \code{fitcond}). See Examples.
#' @param rescale_Px Logical (default FALSE). If TRUE, rescales calculation of Px relative to the fitted value of conductance/PLC at the maximum (least negative) water potential in the dataset. Use this argument only when you know exactly what that means. Identical to \code{rescale_Px} argument in \code{\link{getPx}}.
#' @param shift_zero_min Logical (default FALSE). If TRUE, shifts the water potential data so that the highest (least negative) value measured is set to zero. This has consequences for estimation of Kmax, and is only used for \code{fitcond}. 
#' @param loess_span Only used when \code{model="loess"}, the span parameter setting the desired degree of smoothness (see \code{\link{loess}}).
#' @param msMaxIter Maximum iterations for \code{\link{nlminb}}. Only change when needed.
#'
#' @details 
#' \strong{Models} - 
#' The Weibull model is fit as reparameterized by Ogle et al. (2009), using non-linear regression (\code{\link{nls}}) or a non-linear mixed-effects model if a random effect is present (\code{\link{nlme}}). The sigmoidal-exponential model follows the 
#' specification by Pammenter and van Willigen (1998) : PLC is log-transformed so a linear fit can be obtained with \code{\link{lm}} or \code{\link{lme}} in the presence of a random effect. 
#' Parameters estimated are PX (water potential at which X% conductivity is lost) and SX 
#' (slope of PLC vs. water potential at P50, MPa per percent). For the sigmoidal model, SX is a parameter combination (and so is PX when x is not 50), so only bootstrap estimates of the confidence intervals are given. 
#'
#' \strong{Bootstrap} - 
#' We recommend, where possible, to use the bootstrapped confidence intervals for inference (use at least ca 1000 resamples). For the Weibull model, this is only possible when a sufficiently large sample size is available for a single curve (otherwise too many nls fits will fail). For the sigmoidal model, however, bootstrap is always possible and will always be employed (it cannot be turned off).
#' 
#' \strong{Confidence intervals} - 
#' Calculation of confidence intervals (CI) depends on the method chosen. For the Weibull model, the CI based on profiling ('Normal approximation') is always performed, and a non-parametric bootstrap when \code{bootci=TRUE}. Both are output in \code{coef}, and the bootstrap CI is used in plotting unless otherwise specified (see \code{\link{plot.plcfit}}). When a random effect is specified (for the Weibull model), the CI is calculated with \code{\link{intervals.lme}}. For the sigmoidal model, PX and SX are functions of parameters of a linearized fit, and we thus always use the bootstrap when no random effect is present (it cannot be switched off). When a random effect is included in the sigmoidal model, we use \code{\link{deltaMethod}} from the \code{car} package.
#' 
#' \strong{Weights} - 
#' If a variable with the name Weights is present in the dataframe, this variable will be used as the \code{weights} argument to perform weighted (non-linear) regression. See Examples on how to use this.
#' 
#' \strong{Random effects} - 
#' If the \code{random} argument specifies a factor variable present in the dataframe, random effects will be estimated both for SX and PX. This affects \code{coef} as well as the confidence intervals for the fixed effects. For both the Weibull model and the sigmoidal model, only the random intercept terms are estimated (i.e. \code{random=~1|group}).
#'
#' A plot method is available for the fitted object, see Examples below.
#' @export
#' @importFrom nlme fixef
#' @importFrom nlme nlme
#' @importFrom nlme intervals
#' @importFrom nlme nlmeControl
#' @importFrom car deltaMethod
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics text
#' @importFrom stats lm
#' @importFrom nlme lme
#' @importFrom stats loess
#' @importFrom stats resid
#' @rdname fitplc
#' 
#' @examples
#'
#' # We use the built-in example dataset 'stemvul' in the examples below. See ?stemvul.
#' # Most examples will fit the Weibull model (the default); try running some of the examples
#' # with 'model="sigmoidal"' and compare the results.
#'   
#' # 1. Fit one species (or fit all, see next example)
#' dfr1 <- subset(stemvul, Species =="dpap")
#' 
#' # Fit Weibull model. Store results in object 'pfit'
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
#' # Repeat for the sigmoidal model
#' # Note that varnames specification above is the same as the default, so it 
#' # can be omitted.
#' pfit2 <- fitplc(dfr1, model="sigmoid")
#' plot(pfit2)
#' coef(pfit2)
#' 
#' # 2. Fit all species in the dataset.
#' # Here we also set the starting values (which is sometimes needed).
#' # In this example, we use only 50 bootstrap replicates but recommend you set this
#' # to 1000 or so. 
#' allfit <- fitplcs(stemvul, "Species", varnames=c(PLC="PLC", WP="MPa"), nboot=50)
#' 
#' # 3. Plot the fits.
#' plot(allfit, onepanel=TRUE, plotci=FALSE, px_ci="none", pxlinecol="dimgrey")
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
#' plot(kfits3, onepanel=TRUE, ylim=c(0,12), px_ci="none")
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
                   model=c("Weibull","sigmoidal","loess","nls_sigmoidal"), 
                   x=50,
                   coverage=0.95,
                   bootci=TRUE,
                   nboot=999,
                   quiet=TRUE,
                   startvalues=NULL,
                   shift_zero_min = FALSE,
                   loess_span = 0.7, 
                   msMaxIter = 1000,
                   ...){
    
    
    if(!is.null(startvalues)){
      if(!quiet)warning("startvalues ignored - starting values now estimated from sigmoidal fit.")
    }
  
    model <- match.arg(model)
    
    if(!bootci && model == "sigmoidal"){
      warning("Cannot switch off bootstrap with sigmoidal model - ignored.")
      bootci <- TRUE
    }
    
    if(model == "nls_sigmoidal" && x != 50){
      warning("For nls_sigmoidal, can only set x = 50 (for now anyway).")
      x <- 50
    }
  
    # Find out if called from fitcond.
    mc <- names(as.list(match.call()))
    condfit <- "calledfromfitcond" %in% mc
    
    # Get Kmax value (set in fitcond)
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
      
      if(model == "loess"){
        stop("Cannot estimate random effects with the loess model.")
      }
      if(model == "nls_sigmoidal"){
        stop("Random effects for sigmoidal_nls not yet implemented.")
      }
      
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
    if(any(is.na(c(plc,P))))stop("Missing values found in PLC or WP - remove first!")
    relK <- plc_to_relk(plc)
    
    # Need absolute values of water potential
    if(mean(P) < 0)P <- -P

    # Set least negative to zero, if requested.
    # Useful for very linear data.
    if(shift_zero_min){
      shift_val <- min(P)
      P <- P - shift_val
    } else {
      shift_val <- 0
    }
    
    # weights, if provided
    W <- eval(substitute(weights), dfr)
    if(is.null(W))W <- rep(1, nrow(dfr))
    
    # Dataset tidied
    Data <- data.frame(P=P, PLC=plc, relK=relK, G=G)
    Data$minP <- -Data$P  # negative valued water potential
    
    
    model2 <- paste0(model, ifelse(fitran, "_random", "_fixed"))
    
    fitting_function <- get(model2)
    
    #- fitting_function is one of:
    # weibull_fixed
    # weibull_random
    # loess_fixed
    # sigmoidal_fixed
    # sigmoidal_random
    # nls_sigmoidal_fixed
    #- all fitting functions return a list like:
    #list(fit = fit, pred = pred, cipars = cipars)
    
    structure(c(fitting_function(Data, W, x, coverage, condfit, Kmax,
                                 bootci, nboot, quiet, 
                                 loess_span, msMaxIter),
                
                list(data = Data,
                     x = x,
                     condfit = condfit,
                     Kmax = Kmax,
                     fitran = fitran,
                     bootci = bootci,
                     nboot = nboot,
                     model = model, 
                     coverage = coverage,
                     shiftval = shift_val)),
              class = "plcfit")
    
}    




Weibull_fixed <- function(Data, W, x, coverage, condfit, Kmax,
                          bootci, nboot, quiet, 
                          loess_span, msMaxIter){
  
  # guess starting values from sigmoidal
  f <- do_sigmoid_fit(Data, boot=FALSE, W=W)
  p <- coef(f$fit)
  sp <- sigfit_coefs(p[1], p[2], x=x)
  
  # fit
  Data$X <- x  # Necessary for bootstrap - I think.
  if(!quiet)message("Fitting nls ...", appendLF=FALSE)
  
  fit <- nls(relK ~ fweibull(P,SX,PX,X),
             data=Data, start=list(SX=sp$Sx, PX=sp$Px),
             weights=W)
  
  if(!quiet)message("done.")
  
  inter_val <- ifelse(bootci, "confidence", "none")
  
  pred <- predict_nls(fit, xvarname="P", interval=inter_val, data=Data, 
                      startList=list(SX=sp$Sx, PX=sp$Px), weights=W, 
                      level=coverage,
                      nboot=nboot)
  
  cipars <- try(suppressMessages(confint(fit, level=coverage)), silent=TRUE)
  if(inherits(cipars, "try-error")){
    cipars <- matrix(rep(NA,4),ncol=2)
  }
  cipars <- cbind(coef(fit), cipars)
  dimnames(cipars) <- list(c("SX","PX"), c("Estimate", 
                                           sprintf("Norm - %s",label_lowci(coverage)),
                                           sprintf("Norm - %s",label_upci(coverage))))
  
  if(bootci){
    cisx <- quantile(pred$boot[,"SX"], c((1-coverage)/2, 1 - (1-coverage)/2))
    cipx <- quantile(pred$boot[,"PX"], c((1-coverage)/2, 1 - (1-coverage)/2))
    
    bootpars <- matrix(c(cisx[1],cipx[1],cisx[2],cipx[2]), nrow=2,
                       dimnames=list(c("SX","PX"),
                                     c(sprintf("Boot - %s",label_lowci(coverage)),
                                       sprintf("Boot - %s",label_upci(coverage)))))
    cipars <- cbind(cipars, bootpars)
  }    
  
list(fit = fit, pred = pred, cipars = cipars)
}


Weibull_random <- function(Data, W, x, coverage, condfit, Kmax,
                           bootci, nboot, quiet, 
                           loess_span, msMaxIter){
  
  # guess starting values from sigmoidal
  f <- do_sigmoid_fit(Data, boot=FALSE, W=W)
  p <- coef(f$fit)
  sp <- sigfit_coefs(p[1],p[2],x=x)
  
  Data$X <- x  # Necessary for bootstrap - I think.
  
  fit <- nlme(relK ~ fweibull(P, SX, PX, X),
              fixed=list(SX ~ 1, PX ~ 1),
              random= SX + PX ~ 1|G,
              start=list(fixed=c(SX=sp$Sx, 
                                 PX=sp$Px)),
              control=nlmeControl(msMaxIter = msMaxIter, eval.max=1e06),
              data=Data)
  
  predran <- lapply(split(Data, Data$G), function(group){
    
    ps <- seq_within(group$P, n=101)
    newdat <- data.frame(P=ps, G=unique(group$G), X=x)
    
    list(x=ps, fit=unname(predict(fit, newdat))) 
  })
  
  ps <- seq_within(Data$P, n=101)
  newdat <- data.frame(P=ps, X=x)
  pred <- list(x=ps, fit=predict(fit, newdat, level=0), ran=predran)
  
  cipars <- intervals(fit,which="fixed")$fixed[,c(2,1,3)]
  colnames(cipars) <- c("Estimate", ci_names("Norm", coverage))
  
  attributes(cipars)$label <- NULL
  
list(fit = fit, pred = pred, cipars = cipars)
}



loess_fixed <- function(Data, W, x, coverage, condfit, Kmax,
                        bootci, nboot, quiet, 
                        loess_span, msMaxIter){
  
  Data$W <- W
  fit <- loess(relK ~ P, data=Data, span=loess_span, weights=Data$W, degree=1)
  
  pred <- get_loess_pred(fit, coverage)
  
  ml_Px <- get_px_loessfit(pred, x, rescale = condfit)
  boot_Px <- boot_px_loess(fit, Data, B=999, loess_span, x, rescale = condfit)
  
  cipars <- rbind(c(NA, NA, NA),
                  c(ml_Px, boot_ci(boot_Px, coverage)))
  
  dimnames(cipars) <- list(c("SX","PX"), 
                           c("Estimate", ci_names("Boot",coverage)))
  
list(fit = fit, pred = pred, cipars = cipars)
}



sigmoidal_fixed <- function(Data, W, x, coverage, condfit, Kmax,
                            bootci, nboot, quiet, 
                            loess_span, msMaxIter){
  
  f <- do_sigmoid_fit(Data, boot=TRUE, nboot=nboot, W=W)
  
  # Bootstrap
  cf <- sigfit_coefs(f$boot[,1], f$boot[,2],x=x)
  boot_Sx <- cf$Sx  
  boot_Px <- cf$Px
  
  # Maximum likelihood
  p <- coef(f$fit)
  
  mf <- sigfit_coefs(p[1],p[2],x=x)
  ml_Sx <- mf$Sx
  ml_Px <- mf$Px
  
  # Coefficients matrix
  cipars <- rbind(c(ml_Sx, boot_ci(boot_Sx, coverage)),
                  c(ml_Px, boot_ci(boot_Px, coverage)))
  
  dimnames(cipars) <- list(c("SX","PX"), 
                           c("Estimate", ci_names("Boot",coverage)))
  
  # f must be component with 'fit' and 'boot'
  pred <- get_boot_pred_sigmoid(f, Data, coverage)
  
list(fit = f$fit, pred = pred, cipars = cipars)
}


nls_sigmoidal_fixed <- function(Data, W, x, coverage, condfit, Kmax,
                                bootci, nboot, quiet, 
                                loess_span, msMaxIter){
  
  # guess starting values from sigmoidal
  f <- do_sigmoid_fit(Data, boot=FALSE, W=W)
  p <- coef(f$fit)
  
  # As in sigfit_coefs, but only partial conversion
  a <- p[2]
  b <- p[1]/p[2]
  
  # fit
  Data$X <- x  # Necessary for bootstrap - I think.
  if(!quiet)message("Fitting nls ...", appendLF=FALSE)
  
  # Weighted NLS
  if(!is.null(W)){
    
    fit <- nls(relK ~ 1 - 1/(1 + exp(a*(P - b))),
               data=Data, start=list(a=a, b=b),
               weights=W)
    
  } else {
    
    # Ordinary NLS
    fit <- nls(relK ~ 1 - 1/(1 + exp(a*(P - b))),
               data=Data, start=list(a=a, b=b))
  }
  
  nls_sig_convert_coef <- function(coefs, x){
    
    a <- coefs[1]
    b <- coefs[2]
    Px <- ab_to_px(a, b, x)
    
    # Derivative of sigmoid
    sig2d <- function(Px, a,b)-(exp(a * (Px - b)) * a/(1 + exp(a * (Px - b)))^2)
    Sx <- 100 * sig2d(Px,a,b)
    
    list(Px=Px, Sx=Sx)
  }
  
  fit_ab <- coef(fit)
  sp <- nls_sig_convert_coef(fit_ab, x=x)
  
  if(bootci){
    pred <- predict_nls(fit, xvarname="P", interval="confidence", data=Data, 
                        startList=list(a=fit_ab[1], b=fit_ab[2]), weights=W, 
                        level=coverage,
                        nboot=nboot)
    
    z <- apply(pred$boot, 1, nls_sig_convert_coef, x=x)
    pred$boot <- as.data.frame(do.call(rbind, lapply(z, unlist)))
    names(pred$boot) <- c("Px","Sx")
    
    cisx <- quantile(pred$boot[,"Sx"], c((1-coverage)/2, 1 - (1-coverage)/2))
    cipx <- quantile(pred$boot[,"Px"], c((1-coverage)/2, 1 - (1-coverage)/2))
    
    bootpars <- matrix(c(cisx[1],cipx[1],cisx[2],cipx[2]), nrow=2,
                       dimnames=list(c("SX","PX"),
                                     c(sprintf("Boot - %s",label_lowci(coverage)),
                                       sprintf("Boot - %s",label_upci(coverage)))))
    
    cipars <- matrix(c(sp$Sx, sp$Px), ncol=1, nrow=2)
    dimnames(cipars) <- list(c("SX","PX"), "Estimate")
    cipars <- cbind(cipars, bootpars)
  }
  
list(fit = fit, pred = pred, cipars = cipars)
}


sigmoidal_random <- function(Data, weights, x, coverage, condfit, Kmax,
                             bootci, nboot, quiet, 
                             loess_span, msMaxIter){
  
  # With random effect
  fit <- do_sigmoid_lme_fit(Data)
  
  Px_ci <- deltaMethod(fit, "b0/b1", parameterNames=c("b0","b1"))
  
  # deltaMethod not needed here but convenient and equivalent
  Sx_ci <- deltaMethod(fit, "100*b1/4", parameterNames=c("b0","b1"))
  cipars <- rbind(Sx_ci, Px_ci)
  cipars$SE <- NULL
  dimnames(cipars) <- list(c("SX","PX"),
                           c("Estimate", ci_names("Norm",coverage)))
  
  predran <- lapply(split(Data, Data$G), function(x){
    
    ps <- seq_within(x$minP, n=101)
    newdat <- data.frame(minP=ps, G=unique(x$G))
    
    list(x=-ps, fit=sigmoid_untrans(unname(predict(fit, newdat)))) 
  })
  ps <- seq_within(Data$minP, n=101)
  newdat <- data.frame(minP=ps, X=x)
  pred <- list(x=-ps, fit=predict(fit, newdat, level=0), ran=predran)
  pred$fit <- sigmoid_untrans(pred$fit)
  

list(fit = fit, pred = pred, cipars = cipars)  
}



#----- Subsidiary functions

do_sigmoid_fit <- function(data, W=NULL, boot=FALSE, nboot){
  
  keep <- which(data$PLC < 100 & data$PLC > 0)
  data <- data[keep,]
  W <- W[keep]

  # Transformation as per P&vW
  data$logPLC <- log(100/data$PLC - 1)
  
  lmfit <- lm(logPLC ~ minP, data=data, weights=W)
  if(boot){
    br <- suppressWarnings(bootfit(lmfit, n=nboot, 
                             Data=data, startList=NULL, weights=W))
  } else {
    br <- NA
  }

list(fit=lmfit, boot=br)
}

do_sigmoid_lme_fit <- function(data){
  
  data <- data[data$PLC > 0 & data$PLC < 100,]
  
  data$logPLC <- log(100/data$PLC - 1)
  
  fit <- lme(logPLC ~ minP,
             random= ~1|G,
             data=data)
  
return(fit)
}

# Calculate Sx, Px, given log-linear fit of sigmoidal model
sigfit_coefs <- function(c1,c2,x){
  a <- c2
  b <- c1 / c2
  Px <- ab_to_px(a, b, x)
  
  # Derivative of sigmoid
  sig2d <- function(Px, a,b)-(exp(a * (Px - b)) * a/(1 + exp(a * (Px - b)))^2)
  Sx <- -100 * sig2d(Px,a,b)
  #Sx <- 100 * c2/4  # slope at P50
  
list(Px=unname(Px), Sx=unname(Sx))
}

# Bootstrap predictions for the sigmoidal model
get_boot_pred_sigmoid <- function(f, data, coverage){
  
  preddfr <- data.frame(minP=seq_within(data$minP))
  
  normpred <- sigmoid_untrans(predict(f$fit, preddfr, interval="none"))
  
  bootm <- apply(f$boot,1, function(x)x[1] + x[2]*preddfr$minP)
  bootpred <- as.data.frame(t(apply(bootm, 1, boot_ci, coverage=coverage)))
  names(bootpred) <- c("lwr","upr")
  bootpred <- lapply(bootpred, sigmoid_untrans)
  bootpred$x <- -preddfr$minP
  bootpred$fit <- normpred
  
  return(bootpred)
}




get_loess_pred <- function(fit, coverage){
 
  preddf <- data.frame(P=seq_within(fit$x))
  
  alpha <- 1 - coverage
  qv <- coverage + alpha/2
  
  normpred <- predict(fit, preddf, se=TRUE)
  
  normpred$lwr <- with(normpred, fit - qt(qv, df)*se.fit)
  normpred$upr <- with(normpred, fit + qt(qv, df)*se.fit)
  
return(data.frame(x = preddf$P, 
                  fit = normpred$fit, 
                  lwr = normpred$lwr, 
                  upr = normpred$upr))  
}


# get Px from a fitted and predicted loess model object
get_px_loessfit <- function(pred, x, rescale = FALSE){
  
  if(rescale){
    K0 <- pred$fit[which.min(pred$x)]
    target <- (x/100) * K0
    
    px <- approx(x=pred$fit, y=pred$x, xout=target)$y
    
  } else {
    X <- 1 - x/100
    px <- approx(x=pred$fit, y=pred$x, xout=X)$y
  }
  
return(px)
}

boot_px_loess <- function(fit, Data, B=999, span, x, rescale = FALSE){
  
  rws <- seq(length(resid(fit)))
  
  px <- c()
  for(i in 1:B){
    u <- update(fit, data=Data[sample(rws, replace=TRUE),], span=span)
    p <- get_loess_pred(u, coverage=0.95)
    px[i] <- get_px_loessfit(p, x, rescale = rescale)
  }
  
  # Missing value when Px not actually reached.
  ii <- is.na(px)
  if(any(ii)){
    px <- px[!ii]
    if(sum(ii) > 0.1*B){
      warning("More than 10% of bootstrap replicates produced missing value - CI on Px probably not reliable.")
    }
  } 
  
return(px)
}







