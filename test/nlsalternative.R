

# test optim nls alternative (more flexible)



X <- runif(25,1,20)
Y <- 1.8*X^1.4 + rnorm(25,mean=0.3)

O <- function(par=c(1,1)){
  
  pred <- par[1]*X^par[2]  
  sum((Y-pred)^2)
  
}
opt <- optim(c(1,1),O,hessian=TRUE,method="Nelder")

# magical equation to estimate standard errors:
sqrt(diag(2*opt$value/(length(Y) - 2) * solve(opt$hessian)))

# compare to nls:
summary(nls(Y ~ a*X^b, start=list(a=1,b=1)))$coefficients[,2]




dfr <- read.csv("test/stemvul-ros.csv")
dfr <- subset(dfr, Species =="EuTe")

pfit <- fitplc(dfr, varnames=c(PLC="PLC", WP="MPa"))


Y <- dfr$PLC
relK <- (100 - Y)/100
X <- -dfr$MPa

O <- function(pars=c(22,4)){
  
  pred <- fweibull(X, pars[1], pars[2])
  sum((relK - pred)^2)
  
}
opt <- optim(c(20,5),O,hessian=TRUE) #,method="Nelder")
sqrt(diag(2*opt$value/(length(Y) - 2) * solve(opt$hessian)))





nlmfit <- nlm(O, c(20,5), hessian=TRUE)
sqrt(diag(2*nlmfit$minimum/(length(Y) - 2) * solve(nlmfit$hessian)))  





