




stemros <- read.csv("test/stemvul-ros.csv")

stemros$Cond <- 5.4 * (1 - stemros$PLC/100)

f <- fitcond(stemros, varnames=c(K="Cond", WP="MPa"), WP_Kmax=-0.3)
g <- fitplc(stemros, varnames=c(PLC="PLC", WP="MPa"))

h <- fitcond(stemros, varnames=c(K="Cond", WP="MPa"), WP_Kmax=-0.3, random = Branch)

windows()
par(mfrow=c(2,2))
plot(f)
plot(g, what="embol")
plot(g, what="relk")
plot(h, plotrandom=TRUE)



# Stem at ROS
stemros <- read.csv("test/stemvul-ros.csv")

# fit one, change starting values
dfr <- subset(stemros, Species == "CaCu")
cacufit <- fitplc(dfr, bootci=FALSE, startvalues=list(Px=1.5, S=10))
plot(cacufit, what="embol")


# 
dfr$Weights <- abs(50-dfr$PLC)^1
eutefit1 <- fitplc(dfr,  start=list(Px=4.36,S=22.78 ),
                   weights=Weights)

eutefit1 <- fitplc(dfr,  start=list(Px=4.36,S=22.78 ),
                   weights=NULL)




fit2 <- fitplc(dfr, x=12)
fit3 <- fitplc(dfr, x=50)
getPx(fit3, 12)

# fit all
# no bootstrap
stemrosfits <- fitplcs(stemros, "Species", bootci=FALSE, startvalues=list(Px=1.5, S=10))
stemrosfits
coef(stemrosfits)
summary(stemrosfits)
summary(stemrosfits$CaCu)

# with bootstrap
stemrosfits2 <- fitplcs(stemros, "Species", bootci=TRUE, startvalues=list(Px=1.5, S=10))
stemrosfits2
coef(stemrosfits2)
summary(stemrosfits2)
summary(stemrosfits2$CaCu)

# plot
windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(stemrosfits2[[i]], xlim=c(0,7), main=names(stemrosfits2)[i])

windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(stemrosfits[[i]], xlim=c(0,7), main=names(stemrosfits)[i])

# plot 1
plot(stemrosfits, onepanel=TRUE, plotPx=FALSE, plotdata=TRUE, plotci=FALSE, pch=c(15,19,1))

# plot 2
plot(stemrosfits, onepanel=TRUE, plotPx=FALSE, plotdata=FALSE, 
     legendwhere="topleft",
     plotci=TRUE, what="embol")


# Stemvul
stemvul <- read.csv("test/stemvul-points.csv", sep=";")
f <- fitplcs(stemvul, "Species", bootci=TRUE)

# quiet
g <- fitplcs(stemvul, "Species", bootci=FALSE, quiet=TRUE)

windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
plot(f, xlim=c(0,7), linecol="black") #, main=names(f)[i])


# leafvul
leafvul <- read.csv("test/leafvul-points.csv", sep=";")
leafvulfits <- fitplcs(leafvul, "Species", bootci=TRUE)
windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
plot(leafvulfits, xlim=c(0,7))

# usleafvul
usleafvul <- read.csv("test/usleafvul-points.csv", sep=";")
f <- fitplcs(usleafvul, "Species", bootci=TRUE, 
             varnames=c(PLC="CumulPercent", WP="MPa"))
windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
plot(f, xlim=c(0,7), main=names(f)[i])


plot(f, onepanel=TRUE, xlim=c(0,3), selines="none", linecol=c("#0000FFFF","#13F2FFFF","#FF3300FF"))


# usstemvul
usstemvul <- read.csv("test/usstemvul-points.csv", sep=";")
f <- fitplcs(usstemvul, "Species", bootci=TRUE, 
             varnames=c(PLC="CumulPercent", WP="MPa"))
windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
plot(f, xlim=c(0,7))



# make figure
jpeg("test/plcfitexample.jpg", width=500, height=250, quality=100, pointsize=8)
par(xaxs="i", yaxs="i", cex.lab=2, mar=c(6,6,2,1), cex.axis=1.6)
plot(leafvulfits[[2]], xlim=c(1,7), ylim=c(0,1.1))
dev.off()


jpeg("test/plcfitexample2.jpg", width=368, height=256, quality=100, pointsize=10,
     antialias="cleartype")
par(xaxs="i", yaxs="i", cex.lab=1.2, mar=c(5,5,2,1), cex.axis=1.2)
plot(leafvulfits[[2]], xlim=c(1,7), ylim=c(0,1.1))
dev.off()

windows(5,4)
par(xaxs="i", yaxs="i", cex.lab=1.2, mar=c(5,5,2,1), cex.axis=1.2)
plot(leafvulfits[[2]], xlim=c(1,7), ylim=c(0,1.1))
dev.copy2pdf(file="test/plcfitexample.pdf")


# overlap <- function(vec1, vec2){
#   
#   min1 <- min(vec1)
#   max1 <- max(vec1)
#   vec2 <- c(min(vec2),max(vec2))
#   
#   res <- findInterval(min1, vec2) == 1|| findInterval(max1, vec2) == 1
# return(res)
# }
# 
# x <- c(1.5,3,4)
# y <- c(1,2)
# overlap(x,y)





