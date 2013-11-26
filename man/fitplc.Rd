\name{fitplc}
\alias{fitplc}
\title{Fit PLC curves}
\usage{
  fitplc(dfr, varnames = c(PLC = "PLC", WP = "MPa"),
    model = "Weibull", bootci = FALSE)
}
\arguments{
  \item{dfr}{A dataframe that contains water potential and
  plc data.}

  \item{varnames}{A vector specifying the names of the PLC
  and water potential data (WP) in the dataframe.}
}
\description{
  Description goes here.
}
\examples{
\dontrun{
dfr <- read.csv("test/stemvul-ros.csv")
dfr <- subset(dfr, Species =="EuTe")

pfit <- fitplc(dfr, bootci=TRUE)
pfit
plot(pfit)

}
}

