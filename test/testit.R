

dfr <- read.csv("test/stemvul-ros.csv")
dfr <- subset(dfr, Species =="EuTe")

pfit <- fitplc(dfr, bootci=FALSE)
pfit2 <- fitplc(dfr, bootci=TRUE)
pfit
plot(pfit)
coef(pfit)
 
pfit2
plot(pfit2)
coef(pfit2)




dfr <- read.csv("test/stemvul-ros.csv")

dfrs <- split(dfr, dfr$Species)


# no bootstrap.
pfits <- lapply(dfrs, fitplc)


plot(pfits[[1]])
plot(pfits[[2]])
plot(pfits[[3]])

coef(pfits[[1]])
coef(pfits[[2]])
coef(pfits[[3]])


windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(pfits[[i]], xlim=c(0,7), main=names(pfits)[i])


# only curves:
plot(pfits[[1]], plotdata=F, plotP50=F, linecol="blue")
plot(pfits[[2]], plotdata=F, plotP50=F, add=T, linecol="red")
plot(pfits[[3]], plotdata=F, plotP50=F, add=T, linecol="forestgreen")


pfits2 <- lapply(dfrs, function(x)fitplc(x, bootci=TRUE))
plot(pfits2[[1]], plotdata=F, plotP50=F, linecol="blue")
plot(pfits2[[2]], plotdata=F, plotP50=F, add=T, linecol="red")
plot(pfits2[[3]], plotdata=F, plotP50=F, add=T, linecol="forestgreen")


# something wrong!!!
windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(pfits2[[i]], xlim=c(0,7), main=names(pfits2)[i])

