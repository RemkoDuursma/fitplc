\name{fitplc}
\alias{fitplc}
\alias{fitplcs}
\title{Fit PLC curves}
\usage{
  fitplc(dfr, varnames = c(PLC = "PLC", WP = "MPa"),
    model = "Weibull", bootci = TRUE)

  fitplcs(dfr, group, ...)
}
\arguments{
  \item{dfr}{A dataframe that contains water potential and
  plc data.}

  \item{varnames}{A vector specifying the names of the PLC
  and water potential data (WP) in the dataframe.}

  \item{model}{At the moment, only 'Weibull' is allowed.}

  \item{bootci}{If TRUE, also computes the bootstrap
  confidence interval.}

  \item{group}{Character; variable in the dataframe that
  specifies groups. The curve will be fit for every group
  level.}

  \item{\dots}{Further parameters passed to \code{fitplc}.}
}
\description{
  This function fits 'percent loss conductivity' (PLC)
  curves. At the moment, it only fits the Weibull curve, as
  reparameterized by Ogle et al. (2009), and only returns
  P50 (although soon it will return P88 or whatever
  values). The function \code{fitplcs} can be used for
  batch fitting. See examples below for usage of
  \code{fitplc} or \code{fitplcs}.
}
\examples{
\dontrun{

# First read a dataframe (in this example, from the folder 'test')
dfr <- read.csv("test/stemvul-ros.csv")

# 1. Fit one species (or fit all, see next example)
dfr <- subset(dfr, Species =="EuTe")

# Make fit. Store results in object 'pfit'
pfit <- fitplc(dfr)

# Look at fit
pfit

# Make a standard plot. At the moment, can only plot 'relative conductivity',
# (which is 1.0 where PLC = 0).
plot(pfit)

# Get the coefficients of the fit.
coef(pfit)

# 2. Fit all species in the dataset.
allfit <- fitplcs(dfr, )
}
}

