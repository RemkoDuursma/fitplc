Fit vulnerability curves in R
--------------------------------------

[![](http://www.r-pkg.org/badges/version/fitplc)](http://cran.rstudio.com/web/packages/fitplc/index.html) ![](http://cranlogs.r-pkg.org/badges/fitplc) ![](https://img.shields.io/bitbucket/issues/remkoduursma/fitplc.svg)
![](https://travis-ci.org/RemkoDuursma/fitplc.svg?branch=master) ![](https://img.shields.io/codecov/c/github/RemkoDuursma/fitplc/master.svg)


This page describes the `fitplc` package, which can be used to fit curves to measurements of plant stem, leaf or root conductivity (or conductance) at varying water potentials (so called 'PLC curves'). 
The package either fits the Weibull curve as reparameterized by Ogle et al. (2009), or a sigmoidal model proposed by Pammenter and van Willigen (1998). 

The package calculate confidence intervals for the parameters (e.g., P50, slope at P50) (using the bootstrap or normal approximations), make standard plots of 'relative conductivity' or the 'percentage embolized' vs. water potential (including a bootstrap 95% CI of the fit and the location of the P50).

## Citation

[Duursma R.A., Choat B. 2017. fitplc - an R package to fit hydraulic vulnerability curves. Journal of Plant Hydraulics. doi:10.20870/jph.2017.e002](http://jplanthydro.org/article/view/1541) (Open Access)

## Help

Please read `?fitplc` for more instructions, and the examples on that page.

To report bugs or suggest new features, please [open a new issue by following this link](https://bitbucket.org/remkoduursma/fitplc/issues/new).


## Installation instructions

The `fitplc` package is now on CRAN, so just do

```
install.packages("fitplc")
library(fitplc)
```

To install the development version, use this command:
```
library(devtools)
install_bitbucket("remkoduursma/fitplc")
library(fitplc)
```

Windows users must have [Rtools](http://cran.r-project.org/bin/windows/Rtools/) installed for this to work.

