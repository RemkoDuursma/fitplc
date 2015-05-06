Fit vulnerability curves in R
--------------------------------------

This page describes the `fitplc` package, which can be used to fit PLC curves.
At the moment, simply fits the Weibull curve as reparameterized by Ogle et al. (2009),
computes confidence intervals for P50 (using the bootstrap), and makes a standard plot of 
'relative conductivity' or the 'percentage embolized' vs. water potential (including a bootstrap 95% CI of the fit and the location of the P50).

Please read `?fitplc` for more instructions, and the examples on that page.

More features will be added soon. To let me know which features should be implemented first, please submit an 'Issue', look in the menu above.

## Installation instructions

To install `fitplc`, use this command in R:
```
library(devtools)
install_bitbucket("remkoduursma/fitplc")
library(fitplc)
```

Windows users must have Rtools installed for this to work (http://cran.r-project.org/bin/windows/Rtools/) (Mac users should be OK).


## Example

Example for a leaf vulnerability curve.

![leafvulexample](https://bitbucket.org/remkoduursma/fitplc/raw/master/test/plcfitexample.jpg)