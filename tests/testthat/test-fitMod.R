context("Model Fitting")

source("generate_test_datasets.R")

# Generate data sets and compare results of fitDRModel to the result of nls and
# lm for AIC function (if these are consistent parameter estimates, residual
# sum of square and degrees of freedom are consistent) and the vcov function
# (if these are consistent parameter estimates, RSS, df and gradient are
# consistent)

# TODO:
# * Against what do we compare the following things from testsFitting.R?
#   - predict(fit0, predType="effect-curve", se.fit=TRUE)
#   - predict(fit0, predType="full-model", se.fit=TRUE)
#   - TD(fit0, Delta = 1)
# * Using `unname` to make all.equal shut up about unequal dimnames is a bit ugly
# * What is an acceptable numerical tolerance for expect_equal()?
#   About 1e-4 seems to be the right magnitude (?)
# * sigemax: set.seed(25) as an example where nls and bndnls find different optimum
# * re-set seeds for every model type?


# beta model -------------------------------------------------------------------
set.seed(2000)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)
bnds <- matrix(c(0.05, 0.05, 6, 6), nrow=2)

test_that("the beta model can be fitted (without covariates)", {
  fit0 <- fitMod(x, y, datset, model = "betaMod", addCovars = ~1,
                 addArgs=list(scal=1.2*max(datset$x)), bnds=bnds, start=c(0.6, 0.6))
  fitnls <- nls(y~betaMod(x, e0, eMax, delta1, delta2, 1.2*max(datset$x)),
                start=c(e0=15, eMax=14, delta1=0.8, delta2=0.5), data=datset)
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  expect_equal(coef(fit0), coef(fitnls))
  expect_equal(vcov(fit0), vcov(fitnls))
})

test_that("the beta model can be fitted (with covariates)", {
  fit0 <- fitMod(x, y, datset, model="betaMod", addCovars = ~age+center,
                 addArgs=list(scal=1.2*max(datset$x)), bnds=bnds)
  XX <- model.matrix(~center+age, data=datset)
  scl <- 1.2*max(datset$x)
  fitnls <- nls(y~cbind(XX, betaMod(x, 0, 1, delta1, delta2, scl)),
                data=datset, start=c(delta1=1, delta2=0.2),
                algorithm = "plinear")
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  ord <- c(3, 9, 1, 2, 8, 4, 5, 6, 7)
  expect_equal(unname(coef(fit0)), unname(coef(fitnls))[ord])
  expect_equal(unname(vcov(fit0)), unname(vcov(fitnls))[ord, ord])
})

# emax model -------------------------------------------------------------------
set.seed(15)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)
bnds <- c(1e-5, max(datset$x))

test_that("the emax model can be fitted (without covariates)", {
  fit0 <- fitMod(x,y, datset, model="emax", addCovars = ~1, bnds=bnds)
  fitnls <- nls(y~emax(x, e0, eMax, ed50), start=c(e0=-1, eMax=1.3, ed50=0.1), data=datset)
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  expect_equal(coef(fit0), coef(fitnls))
  expect_equal(vcov(fit0), vcov(fitnls))
})

test_that("the emax model can be fitted (with covariates)", {
  fit0 <- fitMod(x,y, datset, model="emax", addCovars = ~age+center, bnds=bnds)
  XX <- model.matrix(~center+age, data=datset)
  fitnls <- nls(y~cbind(XX, emax(x, 0, 1, ed50)),
                data=datset, start=list(ed50=1), algorithm = "plinear")
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  ord <- c(2, 8, 1, 7, 3, 4, 5, 6)
  expect_equal(unname(coef(fit0)), unname(coef(fitnls))[ord])
  expect_equal(unname(vcov(fit0)), unname(vcov(fitnls))[ord, ord])
})

# sigEmax model ----------------------------------------------------------------
set.seed(13)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)
bnds <- matrix(c(1e-5, 1e-5, max(datset$x), 30), nrow=2)

test_that("the sigEmax model can be fitted (without covariates)", {
  fit0 <- fitMod(x,y, datset, model = "sigEmax", addCovars = ~1, bnds=bnds)
  fitnls <- nls(y~sigEmax(x, e0, eMax, ed50, h),
                start=c(e0=6, eMax=17, ed50=240, h=2), data=datset)
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  expect_equal(coef(fit0), coef(fitnls))
  expect_equal(vcov(fit0), vcov(fitnls))
})

test_that("the sigEmax model can be fitted (with covariates)", {
  fit0 <- fitMod(x,y, datset, model="sigEmax", addCovars = ~age+center, bnds=bnds)
  XX <- model.matrix(~center+age, data=datset)
  fitnls <- nls(y~cbind(XX, sigEmax(x, 0, 1, ed50, h)),
                data=datset, start=list(ed50=368, h=2), algorithm = "plinear")
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  ord <- c(3, 9, 1, 2, 8, 4, 5, 6, 7)
  expect_equal(unname(coef(fit0)), unname(coef(fitnls))[ord])
  expect_equal(unname(vcov(fit0)), unname(vcov(fitnls))[ord, ord])
})

# logistic model ---------------------------------------------------------------
set.seed(200)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)
bnds <- matrix(c(1e-5, 1e-5, max(datset$x), max(datset$x)/2), nrow=2)

test_that("the logistic model can be fitted (without covariates)", {
  fit0 <- fitMod(x,y, datset, model="logistic", addCovars = ~1, bnds=bnds)
  fitnls <- nls(y~logistic(x, e0, eMax, ed50, delta),
                start=c(e0=0, eMax=16, ed50=250, delta=90), data=datset)
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  expect_equal(coef(fit0), coef(fitnls))
  expect_equal(vcov(fit0), vcov(fitnls))
})

test_that("the logistic model can be fitted (with covariates)", {
  fit0 <- fitMod(x,y, datset, model="logistic", addCovars = ~age+center, bnds=bnds)
  XX <- model.matrix(~center+age, data=datset)
  fitnls <- nls(y~cbind(XX, logistic(x, 0, 1, ed50, delta)),
                data=datset, start=list(ed50=220, delta=48), algorithm = "plinear")
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  ord <- c(3, 9, 1, 2, 8, 4, 5, 6, 7)
  expect_equal(unname(coef(fit0)), unname(coef(fitnls))[ord])
  expect_equal(unname(vcov(fit0)), unname(vcov(fitnls))[ord, ord])
})

# exponential model ------------------------------------------------------------
set.seed(4)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)
bnds <- c(0.1, 2)*max(datset$x)

test_that("the exponential model can be fitted (without covariates)", {
  fit0 <- fitMod(x,y, datset, model = "exponential", addCovars = ~1, bnds=bnds)
  fitnls <- nls(y~exponential(x, e0, e1, delta), start=coef(fit0), data=datset)
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  expect_equal(coef(fit0), coef(fitnls))
  expect_equal(vcov(fit0), vcov(fitnls))
})

test_that("the exponential model can be fitted (with covariates)", {
  fit0 <- fitMod(x,y, datset, model = "exponential", addCovars = ~age+center,
                 bnds=bnds)
  XX <- model.matrix(~center+age, data=datset)
  fitnls <- nls(y~cbind(XX, exponential(x, 0, 1, delta)),
                data=datset, start=c(delta=450), algorithm = "plinear")
  expect_equal(AIC(fit0), AIC(fitnls))
  expect_equal(fit0$df, summary(fitnls)$df[2])
  ord <- c(2, 8, 1, 7, 3, 4, 5, 6)
  expect_equal(unname(coef(fit0)), unname(coef(fitnls))[ord])
  expect_equal(unname(vcov(fit0)), unname(vcov(fitnls))[ord, ord])
})

# linear model -----------------------------------------------------------------
# FIXME: seed?
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

test_that("the linear model can be fitted (without covariates)", {
  fit0 <- fitMod(x,y, datset, model = "linear", addCovars = ~1)
  fitlm <- lm(y~x, data=datset)
  expect_equal(AIC(fit0), AIC(fitlm))
  expect_equal(fit0$df, summary(fitlm)$df[2])
  expect_equal(unname(coef(fit0)), unname(coef(fitlm)))
  expect_equal(unname(vcov(fit0)), unname(vcov(fitlm)))
})

test_that("the linear model can be fitted (with covariates)", {
  fit0 <- fitMod(x,y, datset, model = "linear", addCovars = ~age+center)
  fitlm <- lm(y~x+age+center, data=datset)
  expect_equal(AIC(fit0), AIC(fitlm))
  expect_equal(fit0$df, summary(fitlm)$df[2])
  expect_equal(unname(coef(fit0)), unname(coef(fitlm)))
  expect_equal(unname(vcov(fit0)), unname(vcov(fitlm)))
})

# linlog model -----------------------------------------------------------------
# FIXME: seed?
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)
off <- 0.05*max(datset$x)

test_that("the linlog model can be fitted (without covariates)", {
  fit0 <- fitMod(x,y, datset, model = "linlog", addCovars = ~1,addArgs=list(off=off))
  fitlm <- lm(y~log(x+off), data=datset)
  expect_equal(AIC(fit0), AIC(fitlm))
  expect_equal(fit0$df, summary(fitlm)$df[2])
  expect_equal(unname(coef(fit0)), unname(coef(fitlm)))
  expect_equal(unname(vcov(fit0)), unname(vcov(fitlm)))
})

test_that("the linlog model can be fitted (with covariates)", {
  fit0 <- fitMod(x,y, datset, model = "linlog", addCovars = ~age+center,
                 addArgs=list(off=off))
  fitlm <- lm(y~log(x+off)+age+center, data=datset)
  expect_equal(AIC(fit0), AIC(fitlm))
  expect_equal(fit0$df, summary(fitlm)$df[2])
  expect_equal(unname(coef(fit0)), unname(coef(fitlm)))
  expect_equal(unname(vcov(fit0)), unname(vcov(fitlm)))
})

# quadratic model --------------------------------------------------------------
# FIXME: seed?
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

test_that("the quadratic model can be fitted (without covariates)", {
  fit0 <- fitMod(x,y, datset, model = "quadratic", addCovars = ~1)
  fitlm <- lm(y~x+I(x^2), data=datset)
  expect_equal(AIC(fit0), AIC(fitlm))
  expect_equal(fit0$df, summary(fitlm)$df[2])
  expect_equal(unname(coef(fit0)), unname(coef(fitlm)))
  expect_equal(unname(vcov(fit0)), unname(vcov(fitlm)))
})

test_that("the quadratic model can be fitted (with covariates)", {
  fit0 <- fitMod(x,y, datset, model = "quadratic", addCovars = ~age+center)
  fitlm <- lm(y~x+I(x^2)+age+center, data=datset)
  expect_equal(AIC(fit0), AIC(fitlm))
  expect_equal(fit0$df, summary(fitlm)$df[2])
  expect_equal(unname(coef(fit0)), unname(coef(fitlm)))
  expect_equal(unname(vcov(fit0)), unname(vcov(fitlm)))
})


# ------------------------------------------------------------------------------
# ensure that predict with no argument uses the original data not the sorted
# data that were used for fitting

test_that("predict with no argument uses the original data", {
  data(IBScovars)
  ff <- fitMod(dose, resp, data=IBScovars, model="quadratic", addCovars = ~gender)
  expect_equal(predict(ff, predType = "ls-means"),
               predict(ff, predType = "ls-means", doseSeq = IBScovars[,3]))
  expect_equal(predict(ff, predType = "full-model"),
               predict(ff, predType = "full-model", newdata = IBScovars[,-2]))
  expect_equal(predict(ff, predType = "effect-curve"),
               predict(ff, predType = "effect-curve", doseSeq = IBScovars[,3]))
  ff2 <- fitMod(dose, resp, data=IBScovars, model="quadratic")
  expect_equal(predict(ff2, predType = "ls-means"),
               predict(ff2, predType = "ls-means", doseSeq = IBScovars[,3]))
  expect_equal(predict(ff2, predType = "full-model"),
               predict(ff2, predType = "full-model", newdata = IBScovars[,-2]))
  expect_equal(predict(ff2, predType = "effect-curve"),
               predict(ff2, predType = "effect-curve", doseSeq = IBScovars[,3]))
  dose <- unique(IBScovars$dose)
  ord <- c(2,4,1,3,5)
  mns <- as.numeric(tapply(IBScovars$resp, IBScovars$dose, mean)[ord])
  ff3 <- fitMod(dose, mns, S=diag(5), model="quadratic", type = "general")
  expect_equal(predict(ff3, predType = "ls-means"),
               predict(ff3, predType = "ls-means", doseSeq = dose))
  expect_equal(predict(ff3, predType = "effect-curve"),
               predict(ff3, predType = "effect-curve", doseSeq = dose))
})

# ------------------------------------------------------------------------------
# ensure that S is also sorted when the dose is not entered sorted

test_that("S is also sorted when the dose is not entered sorted", {
  data(IBScovars)
  dose <- sort(unique(IBScovars$dose))
  mns <- as.numeric(tapply(IBScovars$resp, IBScovars$dose, mean))
  S <- c(1000,1,1,1,1)*diag(5)
  ff1 <- fitMod(dose, mns, S = S, model="linear", type="general")
  dose <- unique(IBScovars$dose)
  ord <- c(2,4,1,3,5)
  mns <- as.numeric(tapply(IBScovars$resp, IBScovars$dose, mean)[ord])
  ff2 <- fitMod(dose, mns, S = S, model="linear", type="general")
  ff3 <- fitMod(dose, mns, S = S[ord,ord], model="linear", type="general")
  expect_equal(coef(ff1), coef(ff3))
})

test_that("fitMod complains if `resp` is a row-vector", {
  doses <- seq(0, 100, length.out=5)
  resp_col <- emax(doses, 2, 8, 50)
  resp_row <- t(resp_col)
  cov_mat <- diag(0.5, 5)
  fit <- fitMod(doses, resp_col, model = "emax", S = cov_mat,
                type = "general", bnds = defBnds(max(doses))$emax)
  coefs <- unname(coef(fit))
  expect_equal(coefs, c(2, 8, 50), tolerance = 1e-5)
  expect_warning(fitMod(doses, resp_row, model = "emax", S = cov_mat,
                        type = "general", bnds = defBnds(max(doses))$emax),
                 "resp_row is not a numeric but a matrix, converting with as.numeric()")
})
