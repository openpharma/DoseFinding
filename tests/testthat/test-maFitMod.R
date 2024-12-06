library(DoseFinding)
source("maFitMod.R")

data(biom)
anMod <- lm(resp~factor(dose)-1, data=biom)
drFit <- coef(anMod)
S <- vcov(anMod)
doses <- c(0,0.05,0.2,0.6,1)

## Test input parameters
## models
maFitMod(doses, drFit, S, models = c("linear", "emax", "XYZ"))  # should fail
maFitMod(doses, drFit, S, models = c("ABC", "emax", "XYZ"))  # should fail
## bnds
bnds <- defBnds(1)
bnds$emax <- NULL
maFitMod(doses, drFit, S, models = c("emax"), nSim = 10,
         bnds = bnds)  # should show a message

bnds$emax <- c(10,20)
tmp <- maFitMod(doses, drFit, S, models = c("emax"), nSim = 10,
         bnds = bnds)
sapply(tmp$fits, function(x) coef(x)[3]) # all should be in [10,20]
bnds$betaMod <- rbind(c(1,2), c(1,2))
tmp <- maFitMod(doses, drFit, S, models = c("betaMod"), nSim = 10,
         bnds = bnds)
sapply(tmp$fits, function(x) coef(x)[3:4]) # all should be in [1,2]
## addArgs
tmp <- maFitMod(doses, drFit, S, models = c("betaMod"), nSim = 10)
sapply(tmp$fits, function(x) attr(x, "scal")) # should be 1.2
tmp <- maFitMod(doses, drFit, S, models = c("betaMod"), nSim = 10,
                addArgs = list(scal = 200)) 
sapply(tmp$fits, function(x) attr(x, "scal")) # should be 200
tmp <- maFitMod(doses, drFit, S, models = c("linlog"), nSim = 10,
                addArgs = list(off = 123)) 
sapply(tmp$fits, function(x) attr(x, "off")) # should be 123

## test model fitting
fits1 <- maFitMod(doses, drFit, S, models = c("linear"), nSim = 10)
fits2 <- maFitMod(doses, drFit, S, models = c("linear", "sigEmax"), nSim = 10)
fits3 <- maFitMod(doses, drFit, S, models = c("linear", "emax", "betaMod"), nSim = 10)
fits4 <- maFitMod(doses, drFit, S, models = c("linear", "emax", "betaMod"), nSim = 10)
builtin <- c("linlog", "linear", "quadratic", "linInt", "emax",
             "exponential", "logistic", "betaMod", "sigEmax")
fits5 <- maFitMod(doses, drFit, S, models = builtin, nSim = 10)

## test print method
print(fits1)
print(fits2, digits=10)
print(fits3)
print(fits4)
print(fits5)

## test prediction
predict(fits1) # should fail
dsq <- seq(0,1,length=101)
p1 <- predict(fits1, doseSeq = dsq)
p2 <- predict(fits2, doseSeq = dsq)
p3 <- predict(fits3, doseSeq = dsq)
p4 <- predict(fits4, doseSeq = dsq)
p5 <- predict(fits5, doseSeq = dsq)

plot(fits1)
plot(fits2, xlab = "ABC")
plot(fits3, ylab = "XYZ")
plot(fits4)
plot(fits5)

plot(fits1, plotData = "none")
plot(fits2, plotData = "none")
plot(fits3, plotData = "none")
plot(fits4, plotData = "none")
plot(fits5, plotData = "none")

plot(fits1, plotData = "meansCI")
plot(fits2, plotData = "meansCI")
plot(fits3, plotData = "meansCI")
plot(fits4, plotData = "meansCI")
plot(fits5, plotData = "meansCI")

## ED calculation based on median dose-response
ED.maFit(fits5, p = 0.9) # should fail (direction not specified)
ED.maFit(fits5, p = 0.9, direction = "increasing")
ED.maFit(fits5, p = 0.5, direction = "increasing")

## check decreasing direction
drFit2 <- 1-drFit
fits6 <- maFitMod(doses, drFit2, S, models = builtin, nSim = 10)
ED.maFit(fits6, p = 0.9, direction = "increasing")
ED.maFit(fits6, p = 0.9, direction = "decreasing")
ED.maFit(fits6, p = 0.5, direction = "decreasing")
