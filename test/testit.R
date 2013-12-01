

# Stem at ROS
stemros <- read.csv("test/stemvul-ros.csv")
stemrosfits <- fitplcs(stemros, "Species", bootci=TRUE)
stemrosfits

windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(stemrosfits[[i]], xlim=c(0,7), main=names(stemrosfits)[i])
dev.copy2pdf(file="output/stemrosfits.pdf")


# only curves:
windows()
par(mfrow=c(1,1))
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
dev.copy2pdf(file="output/stemvulfits.pdf")


# leafvul
leafvul <- read.csv("test/leafvul-points.csv", sep=";")
leafvulfits <- fitplcs(leafvul, "Species", bootci=TRUE)
windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(leafvulfits[[i]], xlim=c(0,7), main=names(f)[i])
dev.copy2pdf(file="output/leafvulfits.pdf")

# usleafvul
usleafvul <- read.csv("test/usleafvul-points.csv", sep=";")
f <- fitplcs(usleafvul, "Species", bootci=TRUE, 
             varnames=c(PLC="CumulPercent", WP="MPa"))
windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(f[[i]], xlim=c(0,7), main=names(f)[i])
dev.copy2pdf(file="output/usleafvulfits.pdf")


# usstemvul
usstemvul <- read.csv("test/usstemvul-points.csv", sep=";")
f <- fitplcs(usstemvul, "Species", bootci=TRUE, 
             varnames=c(PLC="CumulPercent", WP="MPa"))
windows(10,8)
par(mfrow=c(3,1), mar=c(4,4,2,2))
for(i in 1:3)plot(f[[i]], xlim=c(0,7), main=names(f)[i])
dev.copy2pdf(file="output/usstemvulfits.pdf")



# make figure
jpeg("test/plcfitexample.jpg", width=500, height=250, quality=100, pointsize=8)
par(xaxs="i", yaxs="i", cex.lab=2, mar=c(6,6,2,1), cex.axis=1.6)
plot(leafvulfits[[2]], xlim=c(1,7), ylim=c(0,1.1))
dev.off()




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





