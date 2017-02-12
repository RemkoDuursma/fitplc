


X <- runif(25,1,20)
Y <- 1.8*X^1.4 + rnorm(25,mean=0.3)
nlsfit <- nls(Y ~ a*X^b, start=list(a=2, b=3))


nlsWrapper <- function(X,Y){
  
  fit <- nls(Y ~ a*X^b, start=list(a=2, b=3))
  
return(fit)
}

nlsWrapper2 <- function(X,Y, aStart=2, bStart=4){
  
  func <- function(x, a, b)a*x^b
  
  fit <- nls(Y ~ func(X,a,b), start=list(a=aStart, b=bStart))

  return(fit)
}

