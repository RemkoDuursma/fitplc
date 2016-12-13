

dfr1 <- subset(stemvul, Species =="dpap")

windows(8,6)
with(dfr1, plot(-MPa, PLC, pch=19, col="darkgrey"))

nboot <- 100
pause <- 0.8
last <- NULL

for(i in seq_len(nboot)){
  
  ii <- sample(seq_len(nrow(dfr1)), replace=TRUE)
  dfr <- dfr1[ii,]
  
  with(dfr1, points(-MPa, PLC, pch=19, col="darkgrey"))
  with(dfr, points(-MPa, PLC, pch=19, col="red"))
  
  f <- fitplc(dfr, bootci=FALSE)
  plot(f, add=T, plotci=FALSE, what="PLC", px_ci="none", plotPx=F, plotdata=F,
       linecol="red2")

  Sys.sleep(pause)
  plot(f, add=T, plotci=FALSE, what="PLC", px_ci="none", plotPx=F, plotdata=F,
       linecol="black")
}
plot(fitplc(dfr1), what="PLC")
