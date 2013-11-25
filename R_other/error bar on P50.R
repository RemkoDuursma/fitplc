eleo

setwd("G:/Work/People/Brendan Choat")


# Ogle
f <- function(P, SX, PX, X=50){
  
  V <- (X-100)*log(1-X/100)
  p <- (P/PX)^((PX*SX)/V)
  relk <- (1-X/100)^p
  
  return(relk)
}
curve(f(x, SX=18, PX=6), from=0.1, to=14,
      xlab=expression(-Psi),
      ylab="Relative K")
abline(v=6, lty=3)
abline(h=0.5, lty=3)




# Choat data
plc <- read.csv("PLCthreespecies.csv")

# plc by species
palette(c("blue","red","forestgreen"))
plot(PLC ~ PSI, col=species, data=plc, pch=15)

# relative conductivity
plc$relK <- (100 - plc$PLC)/100


dysox <- subset(plc, species=="Dysoxylum papauanum")
fit_dysox <- nls(relK ~ f(PSI, S, P50),
        start=list(S=5, P50=2.5),
        data=dysox)

eleo <- subset(plc, species=="Eleocarpus grandis")
fit_eleo <- nls(relK ~ f(PSI, S, P50),
                 start=list(S=5, P50=2.5),
                 data=eleo)

syzy <- subset(plc, species=="Syzygium sayeri")
fit_syzy <- nls(relK ~ f(PSI, S, P50),
                start=list(S=35, P50=2.25),
                data=syzy)


# estimate and SE for P50
v <- function(fit){
  c(coef(fit)[[2]], summary(fit)$coefficients[2,2])
}

df <- as.data.frame(rbind(v(fit_dysox), v(fit_eleo), v(fit_syzy)))
rownames(df) <- levels(plc$species)
names(df) <- c("P50","SE")



