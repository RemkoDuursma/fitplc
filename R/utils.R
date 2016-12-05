
relk_to_plc <- function(relk)100 - 100*relk

plc_to_relk <- function(plc)(100 - plc)/100

ab_to_px <- function(a,b,x)(log(1/(1 - x/100) - 1)/a) + b

sigmoid_untrans <- function(x)(100 - 100/(exp(x) + 1))/100

boot_ci <- function(b, coverage){
  
  a <- (1 - coverage)/2
  quantile(b, probs=c(a, 1-a))
  
}


seq_within <- function(v, n=101){
  
  Min <- min(v, na.rm=TRUE)
  Max <- max(v, na.rm=TRUE)
  seq(Min,Max,length=n)
  
}
