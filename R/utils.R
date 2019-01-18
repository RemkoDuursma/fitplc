
# Relative conductance/conductivity to PLC.
relk_to_plc <- function(relk)100 - 100*relk

# PLC to relative conductance/conductivity.
plc_to_relk <- function(plc)(100 - plc)/100

# Convert a and b parameters in P&vW curve to Px.
ab_to_px <- function(a,b,x)(log(1/(1 - x/100) - 1)/a) + b

sigmoid_untrans <- function(x)(100 - 100/(exp(x) + 1))/100

# Derivative of sigmoid
sig2d <- function(Px, a,b)-(exp(a * (Px - b)) * a/(1 + exp(a * (Px - b)))^2)

# Simple parametric bootstrap confidence intervals (quantiles).
boot_ci <- function(b, coverage){
  
  a <- (1 - coverage)/2
  quantile(b, probs=c(a, 1-a))
  
}


seq_within <- function(v, n=101){
  
  Min <- min(v, na.rm=TRUE)
  Max <- max(v, na.rm=TRUE)
  seq(Min,Max,length=n)
  
}


label_lowci <- function(coverage) sprintf("%s%%", 100*(1 - coverage)/2)
label_upci <- function(coverage) sprintf("%s%%", 100*(coverage + (1 - coverage)/2))
label_coverage <- function(coverage)paste0(100*coverage,"%")


ci_names <- function(prefix="", coverage=0.95, sep = " - "){
  c(sprintf("%s%s%s",prefix,
            ifelse(prefix == "", "", sep),
            label_lowci(coverage)),
    sprintf("%s%s%s",prefix,
            ifelse(prefix == "", "", sep),
            label_upci(coverage)))
}

