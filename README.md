Fit vulnerability curves in R
--------------------------------------

This page describes the `fitplc` package, which can be used to fit PLC curves.
At the moment, simply fits the Weibull curve as reparameterized by Ogle et al. (2009),
computes confidence intervals for P50 (using the bootstrap), and makes a standard plot of 
'relative conductivity' vs. water potential (including a bootstrap 95% CI).

More features will be added soon. To let me know which features should be implemented first, please submit an 'Issue', look in the menu above.

## Installation instructions

To install `fitplc`, use this command in R:
```
library(devtools)
install_bitbucket("fitplc","remkoduursma")
library(fitplc)
```

Windows users must have Rtools installed for this to work (http://cran.r-project.org/bin/windows/Rtools/) (Mac users should be OK).


## Example

Example for a leaf vulnerability curve.

![leafvulexample](plcfitexample.jpg)


Below I have simply pasted the example from the help file `?fitplc`. Please read that help file for more instructions. 

```
 # First read a dataframe (in this example, from the folder 'test')
 dfr <- read.csv("test/stemvul-ros.csv")
 
 # 1. Fit one species (or fit all, see next example)
 dfr <- subset(dfr, Species =="EuTe")
 
 # Make fit. Store results in object 'pfit'
 # 'varnames' specifies the names of the 'PLC' variable in the dataframe,
 # and water potential (WP). 
 pfit <- fitplc(dfr, varnames=c(PLC="PLC", WP="MPa"))
 
 # Look at fit
 pfit
 
 # Make a standard plot. At the moment, can only plot 'relative conductivity',
 # (which is 1.0 where PLC = 0).
 plot(pfit)
 
 # Get the coefficients of the fit.
 coef(pfit)
 
 # 2. Fit all species in the dataset.
 allfit <- fitplcs(dfr, "Species", varnames=c(PLC="PLC", WP="MPa"))
 
 # Make three plots
 # windows(10,8) # optional : open up a window and split.
 # par(mfrow=c(3,1), mar=c(4,4,2,2))
 for(i in 1:3)plot(f[[i]], xlim=c(0,7), main=names(f)[i])
 
 # Coefficients show the estimates and 95% CI (given by 'lower' and 'upper')
 # Based on the CI's, species differences can be decided.
 coef(allfit)
```
