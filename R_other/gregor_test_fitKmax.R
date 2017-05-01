
greg <- read.csv("R_other/gregor_13_species_set1.csv")

library(ggplot2)

ggplot(greg, aes(MPa,K)) + geom_point() + facet_wrap(~Species)


dat <- subset(greg, Species == "DI")

library(fitplc)

dat$MPa_re <- with(dat, MPa - max(MPa))
vn <- c(K="K", WP="MPa_re")

f <- fitcond(dat, Kmax=34.8, varnames=vn, model="sigm")
plot(f)

O <- function(kmax, ...){
  
  f <- fitcond(dat, Kmax=kmax, boot=TRUE, model="sigmoidal", varnames=vn)
summary(f$fit)$sigma

}
Ov <- Vectorize(O)


  
f <- fitcond(dat, Kmax=60, varnames=vn, model="sigm", rescale_Px = TRUE)

f$pred$fit[f$pred$x == 0]
plot(f)



curve(Ov(x), from=30, to=60, n=50)

g <- fitcond(dat, Kmax=60, boot=FALSE, varnames=vn)
plot(g)


library(nlshelper)
l <- loess(K ~ MPa, data=dat, span=0.9)
plot_loess(l, xlim=c(-6,0))




f <- fitcond(dat, WP_Kmax=1)
f2 <- fitcond(dat, WP_Kmax=1, shift_zero_min=TRUE)
plot(f2)

f3 <- fitcond(dat, WP_Kmax=1, shift_zero_min=TRUE, recale_Px=TRUE)
plot(f3)


O <- function(kmax, ...){
  
  f <- fitcond(dat, Kmax=kmax, model="sigmoidal", shift_zero_min=TRUE)
  summary(f$fit)$sigma
  
}
Ov <- Vectorize(O)


system.time(opt <- optimize(Ov, c(1,10^3)))



curve(Ov(x), from=25, to=50, n=50)


f <- fitcond(dat, Kmax=100, model="sigmoidal", shift_zero_min=TRUE, rescale_Px = TRUE)
plot(f)


vals <- seq(30, 200, by=10)
l <- list()
for(i in seq_along(vals)){
  l[[i]] <- fitcond(dat, Kmax=vals[i], model="sigmoidal", shift_zero_min=TRUE, rescale_Px = TRUE)
}

k0 <- sapply(l, "[", "K0")
plot(vals, k0)

px <- sapply(l, function(x)coef(x)["PX","Estimate"])

sig <- sapply(l, function(x)summary(x$fit)$sigma)




g <- fitcond(dat, shift_zero_min=TRUE, Kmax=30, bootci=FALSE)

vals <- seq(20, 60, by=1)
l <- list()
for(i in seq_along(vals)){
  l[[i]] <- fitcond(dat, Kmax=vals[i], bootci=FALSE, shift_zero_min=TRUE)
}
px <- sapply(l, function(x)coef(x)["PX","Estimate"])

sig <- sapply(l, function(x)summary(x$fit)$sigma)

plot(vals, sig)


