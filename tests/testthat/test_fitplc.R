library(fitplc)
context("Fit PLC curves")

dpap <- subset(stemvul, Species == "dpap")

# Don't change - tests won't make sense
Xval <- 50
cover <- 0.95
n_boot <- 101

f <- fitplc(dpap, x=Xval, coverage=cover, nboot=n_boot)
f2 <- fitplc(dpap, x=Xval, coverage=cover, nboot=n_boot, shift_zero_min=0.1)

g <- fitplc(dpap, bootci=FALSE, x=Xval, coverage=cover, nboot=n_boot)
h <- fitplc(dpap, model="sigm", x=Xval, coverage=cover, nboot=n_boot)
h4 <- fitplc(dpap, model="loess", x=Xval, coverage=cover, nboot=n_boot)

k <- fitplc(stemvul, random=Species, x=Xval, coverage=cover, nboot=n_boot)
m <- fitplc(stemvul, random=Species, model="sigm", x=Xval, coverage=cover, nboot=n_boot)
m0 <- fitplc(stemvul, model="sigm", x=Xval, coverage=cover, nboot=n_boot)

h1 <- fitplcs(stemvul, group="Species", x=Xval, coverage=cover, nboot=n_boot)
h2 <- fitplcs(stemvul, group="Species", x=Xval, model="sigmoid", coverage=cover, nboot=n_boot)

fc1 <- fitcond(dpap, WP_Kmax=0.5, varnames=c(K="Cond", WP="MPa"), 
               x=Xval, coverage=cover, nboot=n_boot)
fc2 <- fitcond(dpap, WP_Kmax=0.5, varnames=c(K="Cond", WP="MPa"), 
               model="sigmoid", x=Xval, coverage=cover, nboot=n_boot)
fc3 <- fitcond(dpap, WP_Kmax=-0.5, varnames=c(K="Cond", WP="MPa"), 
               model="loess", x=Xval, coverage=cover, nboot=n_boot,
               rescale_Px = TRUE)
h3 <- fitconds(stemvul, group="Species", 
               WP_Kmax=0.5, varnames=c(K="Cond", WP="MPa"), 
               x=Xval, model="sigmoid", coverage=cover, nboot=n_boot)

s1 <- fitplc(dpap, x=Xval, coverage=cover, nboot=n_boot, 
             model = "nls_sigmoidal")

# More test data
test_cav <- read.table(text="Species    WP       PLC
      B3 -0.50  0.000000
      B3 -0.80  1.019108
      B3 -1.10  5.987261
      B3 -1.50  6.242038
      B3 -2.00  8.789809
      B3 -2.30 60.636943
      B3 -2.49 90.522293
      B3 -2.79 98.598726", header=TRUE)
s2 <- fitplc(test_cav, model= "nls_sigmoidal", varnames=c(WP="WP", PLC="PLC"), boot=FALSE)
s3 <- fitplc(test_cav, model= "nls_sigmoidal", varnames=c(WP="WP", PLC="PLC"), boot=TRUE)


dpap$Weights <- abs(50-dpap$PLC)^1.2
w1 <- fitplc(dpap, model = "Weibull", weights=Weights, nboot=n_boot)
w2 <- fitplc(dpap, model = "sigmoidal", weights=Weights, nboot=n_boot)

# getPx
fun <- function(object)print(getPx(object, x=c(12,50,88)))
gpx <- lapply(list(f,g,h,k,m,h1,h3),fun)


test_that("Estimated coefficients", {
  expect_gt(coef(f)["SX","Estimate"], 25)
  expect_lt(coef(f)["SX","Estimate"], 30)
  expect_gt(coef(f)["PX","Estimate"], 2.5)
  expect_lt(coef(f)["PX","Estimate"], 2.7)
  
  expect_equal(coef(f)[,c("Estimate","Norm - 2.5%","Norm - 97.5%")],
               coef(g)[,c("Estimate","Norm - 2.5%","Norm - 97.5%")])
  
  expect_equal(nrow(f$pred$boot), n_boot)
  
  expect_gt(coef(h)["SX","Estimate"], 24)
  expect_lt(coef(h)["SX","Estimate"], 29)
  expect_gt(coef(h)["PX","Estimate"], 2.5)
  expect_lt(coef(h)["PX","Estimate"], 3.0)
  
  expect_equal(colnames(coef(g)), colnames(coef(k)))
  
  expect_gt(coef(k)["SX","Estimate"], 30)
  expect_lt(coef(k)["SX","Estimate"], 40)
  expect_gt(coef(k)["PX","Estimate"], 2.5)
  expect_lt(coef(k)["PX","Estimate"], 3.0)
  
  expect_equal(colnames(coef(k)), colnames(coef(m)))
  
  expect_gt(mean(coef(h2)[seq(1,5,by=2),"Estimate"]), 25)
  expect_lt(mean(coef(h2)[seq(1,5,by=2),"Estimate"]), 30)
  expect_gt(mean(coef(h2)[seq(2,6,by=2),"Estimate"]), 2.5)
  expect_lt(mean(coef(h2)[seq(2,6,by=2),"Estimate"]), 3.0)

  expect(all(colnames(coef(fc2)) %in% colnames(coef(fc1))), 
         "Coefficient names not consistent")

  expect_gt(coef(fc1)["SX","Estimate"], 15)
  expect_lt(coef(fc1)["SX","Estimate"], 40)
  expect_gt(coef(fc1)["PX","Estimate"], 2.5)
  expect_lt(coef(fc1)["PX","Estimate"], 3.5)
  
  expect_gt(coef(fc2)["SX","Estimate"], 15)
  expect_lt(coef(fc2)["SX","Estimate"], 40)
  expect_gt(coef(fc2)["PX","Estimate"], 2.5)
  expect_lt(coef(fc2)["PX","Estimate"], 3.5)
  
})


test_that("Breaking things", {
  expect_error(fitplcs(stemvul, group="NOGROUP", x=Xval, coverage=cover, nboot=n_boot))
  expect_error(fitconds(stemvul, group="NOGROUP", x=Xval, coverage=cover, nboot=n_boot))
  expect_warning(fitplc(dpap, x=Xval, coverage=cover, nboot=n_boot, quiet=FALSE,
                        start=list(a=1, b=2)))
  expect_warning(fitplc(dpap, model="sigm", x=Xval, coverage=cover, nboot=n_boot, boot=FALSE))
  expect_error(fitplc(dpap, model="sigm", x=Xval, coverage=cover, nboot=n_boot, 
                      varnames = c(PLC = "NOT", WP = "MPa")))
  expect_error(fitplc(dpap, model="sigm", x=Xval, coverage=cover, nboot=n_boot, 
                      varnames = c(PLC = "PLC", WP = "NOT")))
  expect_error(fitplc(stemvul, random=Species, model="loess",
                      x=Xval, coverage=cover, nboot=n_boot))
  expect_error(fitplc(stemvul, random=Species, model="nls_sigmoidal",
                      x=Xval, coverage=cover, nboot=n_boot))
  expect_message(fitplc(stemvul, random=Species, quiet=FALSE, bootci=TRUE,
                      x=Xval, coverage=cover, nboot=n_boot))
})

context("Print methods")
print(f)
print(g)
print(h)
print(h4)
print(h1)
print(h2)
print(h3)
print(fc1)
print(fc2)
print(fc3)
summary(f) # equals (print(f))
print(k)

context("Plot curves")
# Not sure how to test except to run a few bits.

curve(fsigmoidal(x, 2,10), from=0, to=5)

plot(f)
plot(g, add=TRUE,px_ci_type="horizontal")

plot(f, what = "embol")

plot(w1, xaxis="negative")
plot(f, citype = "lines")
plot(g, px_ci = "bootstrap")

plot(f, xlab="Hello", ylab="Hey", ylim=c(0,2), pch=15, plotPx=FALSE,
     plotci=FALSE, plotdata=FALSE, plotfit=FALSE, multiplier=1.01,
     px_ci="parametric", px_ci_type="horizontal", px_ci_label=FALSE,
     linecol="grey",linetype=3, linecol2="red",pxlinecol="blue",
     pxcex=0.9, citype="lines", what="PLC",xaxis="negative",
     main="intentionally blank figure")

plot(f, px_ci="parametric")
plot(h, px_ci="parametric")

plot(fc1)

plot(k, plotrandom=TRUE)
plot(m, plotrandom=TRUE)

# Jen bug, 2018-05-29
plot(k, what="PLC", plotrandom=TRUE)
plot(h1, what="embol")
plot(h1, onepanel=TRUE)

