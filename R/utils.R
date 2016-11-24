
relk_to_plc <- function(relk)100 - 100*k

plc_to_relk <- function(plc)(100 - plc)/100


ab_to_px <- function(a,b,x)(1/a)*(50/x - 1) + b


boot_ci <- function(b, coverage){
  
  a <- (1 - coverage)/2
  quantile(b, probs=c(a, 1-a))
  
}
