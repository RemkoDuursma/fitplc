
library(fitplc)



dpap <- subset(stemvul, Species == "dpap")

Xval <- 50

f <- fitplc(dpap, x=Xval)
g <- fitplc(dpap, bootci=FALSE, x=Xval)
h <- fitplc(dpap, model="sigm", x=Xval)

k <- fitplc(stemvul, random=Species, x=Xval)
m <- fitplc(stemvul, random=Species, model="sigm", x=Xval)
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



