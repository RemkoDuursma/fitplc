library(fitplc)
context("Fit PLC curves")

dpap <- subset(stemvul, Species == "dpap")

# Don't change - tests won't make sense
Xval <- 50
cover <- 0.95
n_boot <- 101

f <- fitplc(dpap, x=Xval, coverage=cover, nboot=n_boot)
g <- fitplc(dpap, bootci=FALSE, x=Xval, coverage=cover, nboot=n_boot)
h <- fitplc(dpap, model="sigm", x=Xval, coverage=cover, nboot=n_boot)

k <- fitplc(stemvul, random=Species, x=Xval, coverage=cover, nboot=n_boot)
m <- fitplc(stemvul, random=Species, model="sigm", x=Xval, coverage=cover, nboot=n_boot)
m0 <- fitplc(stemvul, model="sigm", x=Xval, coverage=cover, nboot=n_boot)

h1 <- fitplcs(stemvul, group="Species", x=Xval, coverage=cover, nboot=n_boot)
h2 <- fitplcs(stemvul, group="Species", x=Xval, model="sigmoid", coverage=cover, nboot=n_boot)

fc1 <- fitcond(dpap, WP_Kmax=0.5, varnames=c(K="Cond", WP="MPa"), 
               x=Xval, coverage=cover, nboot=n_boot)
fc2 <- fitcond(dpap, WP_Kmax=0.5, varnames=c(K="Cond", WP="MPa"), 
               model="sigmoid", x=Xval, coverage=cover, nboot=n_boot)


# getPx
fun <- function(object)print(getPx(object, x=c(12,50,88)))
gpx <- lapply(list(f,g,h,k,m),fun)


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



context("Plot curves")


# Not sure how to test except to run a few bits.
plot(f)
plot(g, add=TRUE,px_ci_type="horizontal")

plot(f, xlab="Hello", ylab="Hey", ylim=c(0,2), pch=15, plotPx=FALSE,
     plotci=FALSE, plotdata=FALSE, plotfit=FALSE, multiplier=1.01,
     px_ci="parametric", px_ci_type="horizontal", px_ci_label=FALSE,
     linecol="grey",linetype=3, linecol2="red",pxlinecol="blue",
     pxcex=0.9, citype="lines", what="PLC",xaxis="negative")

plot(f, px_ci="parametric")
plot(h, px_ci="parametric")







