
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
