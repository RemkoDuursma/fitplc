
library(fitplc)



dpap <- subset(stemvul, Species == "dpap")

Xval <- 50

f <- fitplc(dpap, x=Xval)
g <- fitplc(dpap, bootci=FALSE, x=Xval)
h <- fitplc(dpap, model="sigm", x=50)

k <- fitplc(stemvul, random=Species, x=Xval)
m <- fitplc(stemvul, random=Species, model="sigm", x=50)
m0 <- fitplc(stemvul, model="sigm", x=Xval)


plot(f)
plot(f, what="PLC")
plot(g)
plot(g, what="PLC")
plot(h)
plot(h, what="PLC")



plot(k, plotrandom=TRUE)
plot(k, plotrandom=FALSE)
coef(k)

plot(m, plotrandom=F)
plot(m, plotrandom=T)
coef(m)

plot(m0)
coef(m0)



fc1 <- fitcond(dpap, WP_Kmax=0.5, varnames=c(K="Cond", WP="MPa"))
plot(fc1)
coef(fc1)
summary(fc1)

fc2 <- fitcond(dpap, WP_Kmax=0.5, varnames=c(K="Cond", WP="MPa"), model="sigmoid")
plot(fc2)
coef(fc2)
summary(fc2)



h1 <- fitplcs(stemvul, group="Species")
plot(h1, onepanel=TRUE, selines="none", plotci=FALSE)
abline(h=0.5)
coef(h1)

h2 <- fitplcs(stemvul, group="Species", model="sigmoid")
plot(h2, onepanel=TRUE, selines="none", plotci=FALSE)
abline(h=0.5)
coef(h2)


pira <- read.csv("c:/repos/fitplcpaper/data/Pinus_radiata_centrifuge.csv")
f4 <- fitplc(pira, varnames=c(PLC="PLC",WP="Psi"),
             random=Rep)


