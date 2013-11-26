

# Stem at ROS
stemros <- read.csv("test/stemvul-ros.csv")
stemrosfits <- fitplcs(stemros, "Species", bootci=TRUE)

windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(stemrosfits[[i]], xlim=c(0,7), main=names(stemrosfits)[i])
dev.copy2pdf(file="output/stemrosfits.pdf")


# only curves:
plot(stemrosfits[[1]], plotdata=F, plotP50=F, linecol="blue")
plot(stemrosfits[[2]], plotdata=F, plotP50=F, add=T, linecol="red")
plot(stemrosfits[[3]], plotdata=F, plotP50=F, add=T, linecol="forestgreen")
legend("topright", names(stemrosfits), fill=c("blue","red","forestgreen"))
dev.copy2pdf(file="output/stemrosfits_curvesonly.pdf")


# Stemvul
stemvul <- read.csv("test/stemvul-points.csv", sep=";")
f <- fitplcs(stemvul, "Species", bootci=TRUE)

windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(f[[i]], xlim=c(0,7), main=names(f)[i])
dev.copy2pdf(file="output/leafvulfits.pdf")


# leafvul
leafvul <- read.csv("test/leafvul-points.csv", sep=";")
f <- fitplcs(leafvul, "Species", bootci=TRUE)
windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(f[[i]], xlim=c(0,7), main=names(f)[i])
dev.copy2pdf(file="output/leafvulfits.pdf")





