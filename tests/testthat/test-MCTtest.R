context("multiple contrast test")

# TODO:
# * think about suitable tolerances
# * maybe define common candidate models outside of test_that() calls?
# * how do we check for equal p-values (calculated with MC algorighm)?
# * the non-psd covariance matrix is actually psd (!?)
# * pull shared code out of test_that() calls

source("generate_test_datasets.R")

require_multcomp <- function() {
  if (!require("multcomp")) {
    skip("multcomp package not available")
  }
}

# helper functions to increase readability of expect_equal() calls
tstat <- function(obj) {
  UseMethod("tstat")
}
tstat.MCTtest <- function(obj) {
  # drop the pVal attribute of obj$tStat
  as.numeric(obj$tStat)
}
tstat.glht <- function(obj) {
  unname(summary(obj)$test$tstat)
}
pval <- function(obj) {
  UseMethod("pval")
}
pval.MCTtest <- function(obj) {
  attr(obj$tStat, "pVal")
}
pval.glht <- function(obj) {
  as.numeric(summary(obj)$test$pvalues)
}

test_that("MCTtest gives the same output as multcomp::glht (beta and sigEmax models)", {
  require_multcomp()
  set.seed(10)
  dd <- getDFdataSet_testsMCT()
  dd_x_factor <- dd
  dd_x_factor$x <- as.factor(dd$x)
  bet <- guesst(0.9*max(dd$x), p=0.8, "betaMod", scal = 1.2*max(dd$x),
                dMax = 0.7*max(dd$x), Maxd = max(dd$x))
  sE <- guesst(c(0.5*max(dd$x), 0.7*max(dd$x)) , p=c(0.5, 0.9), "sigEmax")
  models <- Mods(linear = NULL, betaMod = bet, sigEmax = sE,
                 doses = sort(unique(dd$x)),
                 addArgs=list(scal = 1.2*max(dd$x)))
  # model with covariates
  obj <- MCTtest(x,y, dd, models=models, addCovars = ~cov1+cov2, pVal = TRUE)
  fit <- lm(y~x+cov1+cov2, data=dd_x_factor)
  mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
  expect_equal(tstat(obj), tstat(mcp))
  expect_equal(pval(obj), pval(mcp))
  # model without covariates
  obj <- MCTtest(x,y, dd, models=models, addCovars = ~1, pVal = TRUE)
  fit <- lm(y~x, data=dd_x_factor)
  mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
  expect_equal(tstat(obj), tstat(mcp))
  expect_equal(pval(obj), pval(mcp))
})

test_that("MCTtest gives the same output as multcomp::glht (logistic, exponential, quadratic models)", {
  require_multcomp()
  set.seed(10)
  dd <- getDFdataSet_testsMCT()
  dd_x_factor <- dd
  dd_x_factor$x <- as.factor(dd$x)
  mD <- max(dd$x)
  lg1 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.9), "logistic")
  lg2 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.5), "logistic")
  expo <- guesst(c(0.9*mD), c(0.7), "exponential", Maxd=mD)
  quad <- guesst(c(0.6*mD), c(1), "quadratic")
  models <- Mods(linlog = NULL, logistic = rbind(lg1, lg2),
                 exponential = expo, quadratic = quad,
                 doses = sort(unique(dd$x)), addArgs=list(off = 0.2*max(dd$x)))
  # model with covariates
  obj <- MCTtest(x, y, dd, models=models, addCovars = ~cov1+cov2, pVal = TRUE)
  fit <- lm(y~x+cov1+cov2, data=dd_x_factor)
  mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
  expect_equal(tstat(obj), tstat(mcp))
  expect_equal(pval(obj), pval(mcp))
  # model without covariates
  obj <- MCTtest(x,y, dd, models=models, addCovars = ~1, pVal = TRUE)
  fit <- lm(y~x, data=dd_x_factor)
  mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
  expect_equal(tstat(obj), tstat(mcp))
  expect_equal(pval(obj), pval(mcp))
})

test_that("MCTtest works with contrast matrix handed over", {
  require_multcomp()
  set.seed(23)
  dd <- getDFdataSet_testsMCT()
  mD <- max(dd$x)
  lg1 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.9), "logistic")
  lg2 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.5), "logistic")
  expo <- guesst(c(0.9*mD), c(0.7), "exponential", Maxd=mD)
  quad <- guesst(c(0.6*mD), c(1), "quadratic")
  models <- Mods(linlog = NULL, logistic = rbind(lg1, lg2),
                 exponential = expo, quadratic = quad,
                 doses = dd$x, addArgs=list(off = 0.2*max(dd$x)))
  contMat <- MCTtest(x,y, dd, models=models, addCovars = ~cov1+cov2, pVal = TRUE)$contMat
  obj <- MCTtest(x,y, dd, models=models, addCovars = ~1, pVal = TRUE, contMat = contMat)
  dd$x <- as.factor(dd$x)
  fit <- lm(y~x, data=dd)
  mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
  expect_equal(tstat(obj), tstat(mcp))
  expect_equal(pval(obj), pval(mcp))
})

test_that("MCTtest works with binary data (1)", {
  require_multcomp()
  set.seed(1909)
  dd <- getDFdataSet.bin()
  bet <- guesst(0.9*max(dd$x), p=0.8, "betaMod", scal = 1.2*max(dd$x), dMax = 0.7*max(dd$x),
                Maxd = max(dd$x))
  sE <- guesst(c(0.5*max(dd$x), 0.7*max(dd$x)) , p=c(0.5, 0.9), "sigEmax")
  models <- Mods(linear = NULL, betaMod = bet, sigEmax = sE,
                 doses = sort(unique(dd$x)), addArgs=list(scal = 1.2*max(dd$x)))
  logReg <- glm(y~as.factor(x)-1, family=binomial, data=dd, weights = n)
  dePar <- coef(logReg)
  vCov <- vcov(logReg)
  dose <- sort(unique(dd$x))
  obj <- MCTtest(dose, dePar, S=vCov, models=models, type="general",
                 df=Inf, pVal = TRUE)
  dd$x <- as.factor(dd$x)
  fit <- glm(y~x-1, family = binomial, data=dd, weights = n)
  mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
  expect_equal(tstat(obj), tstat(mcp))
  expect_equal(pval(obj), pval(mcp))
})

test_that("MCTtest works with binary data (2)", {
  require_multcomp()
  set.seed(1997)
  dd <- getDFdataSet.bin()
  bet <- guesst(0.9*max(dd$x), p=0.8, "betaMod", scal = 1.2*max(dd$x),
                dMax = 0.7*max(dd$x), Maxd = max(dd$x))
  sE <- guesst(c(0.5*max(dd$x), 0.7*max(dd$x)) , p=c(0.5, 0.9), "sigEmax")
  models <- Mods(linear = NULL, betaMod = bet, sigEmax = sE,direction = "decreasing",
                 addArgs=list(scal = 1.2*max(dd$x)), doses = sort(unique(dd$x)))
  logReg <- glm(y~as.factor(x)-1, family=binomial, data=dd, weights = n)
  dePar <- coef(logReg)
  vCov <- vcov(logReg)
  dose <- sort(unique(dd$x))
  obj <- MCTtest(dose, dePar, S=vCov, models=models, type = "general",
                 pVal = TRUE, df=Inf)
  dd$x <- as.factor(dd$x)
  fit <- glm(y~x-1, family = binomial, data=dd, weights = n)
  mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
  expect_equal(tstat(obj), tstat(mcp))
  expect_equal(pval(obj), pval(mcp))
})

test_that("MCTtest works with binary data (3)", {
  require_multcomp()
  set.seed(1)
  dd <- getDFdataSet.bin()
  bet <- guesst(0.9*max(dd$x), p=0.8, "betaMod", scal = 1.2*max(dd$x),
                dMax = 0.7*max(dd$x), Maxd = max(dd$x))
  sE <- guesst(c(0.5*max(dd$x), 0.7*max(dd$x)) , p=c(0.5, 0.9), "sigEmax")
  models <- Mods(linear = NULL, betaMod = bet, sigEmax = sE,
                 doses = sort(unique(dd$x)), addArgs=list(scal = 1.2*max(dd$x)))
  logReg <- glm(y~as.factor(x)-1, family=binomial, data=dd, weights = n)
  dePar <- coef(logReg)
  vCov <- vcov(logReg)
  dose <- sort(unique(dd$x))
  obj <- MCTtest(dose, dePar, S=vCov, models=models, type = "general",
                 pVal = TRUE, df=Inf)
  dd$x <- as.factor(dd$x)
  fit <- glm(y~x-1, family = binomial, data=dd, weights = n)
  mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
  expect_equal(tstat(obj), tstat(mcp))
  expect_equal(pval(obj), pval(mcp))
})

test_that("a one-dimensional test works", {
  require_multcomp()
  set.seed(1)
  dd <- getDFdataSet.bin()
  model <- Mods(linear = NULL, doses=sort(unique(dd$x)), addArgs=list(scal = 1.2*max(dd$x)))
  logReg <- glm(y~as.factor(x)-1, family=binomial, data=dd, weights = n)
  dePar <- coef(logReg)
  vCov <- vcov(logReg)
  dose <- sort(unique(dd$x))
  expect_warning(obj <- MCTtest(dose, dePar, S=vCov, models=model, type = "general", pVal = TRUE, df=Inf),
                 "univariate: using pnorm")
  dd$x <- as.factor(dd$x)
  fit <- glm(y~x-1, family = binomial, data=dd, weights = n)
  mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
  expect_equal(tstat(obj), tstat(mcp))
  expect_equal(pval(obj), pval(mcp))
})

test_that("unordered values in MCTtest work (placebo adjusted scale)", {
  require_multcomp()
  data(IBScovars)
  modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                  linInt = c(0, 1, 1, 1), doses = c(0, 1, 2, 3, 4))
  ancMod <- lm(resp~factor(dose)+gender, data=IBScovars)
  drEst <- coef(ancMod)[2:5]
  vc <- vcov(ancMod)[2:5, 2:5]
  doses <- 1:4
  fit_orig <- fitMod(doses, drEst, S=vc, model = "sigEmax", placAdj=TRUE, type = "general")
  test_orig <- MCTtest(doses, drEst, S = vc, models = modlist, placAdj = TRUE,
                       type = "general", df = Inf)
  ord <- c(3,4,1,2)
  drEst2 <- drEst[ord]
  vc2 <- vc[ord,ord]
  doses2 <- doses[ord]
  fit_perm <- fitMod(doses2, drEst2, S=vc2, model = "sigEmax", placAdj=TRUE, type = "general")
  test_perm <- MCTtest(doses2, drEst2, S = vc2, models = modlist, placAdj = TRUE,
                       type = "general", df = Inf)
  # we don't compare stuff we want to be different
  attr(fit_orig, "data") <- attr(fit_perm, "data") <- NULL
  attr(fit_orig, "doseRespNam") <- attr(fit_perm, "doseRespNam") <- NULL
  expect_equal(fit_orig, fit_perm)
  expect_equal(test_orig, test_perm)
})

test_that("unordered values in MCTtest work (unadjusted scale)", {
  require_multcomp()
  data(IBScovars)
  modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                  linInt = c(0, 1, 1, 1), doses = c(0, 1, 2, 3, 4))
  ancMod <- lm(resp~factor(dose)-1, data=IBScovars)
  drEst <- coef(ancMod)
  vc <- vcov(ancMod)
  doses <- 0:4
  fit_orig <- fitMod(doses, drEst, S=vc, model = "sigEmax", type = "general")
  test_orig <- MCTtest(doses, drEst, S = vc, models = modlist, type = "general", df = Inf)
  ord <- c(3,4,1,2,5)
  drEst2 <- drEst[ord]
  vc2 <- vc[ord,ord]
  doses2 <- doses[ord]
  fit_perm <- fitMod(doses2, drEst2, S=vc2, model = "sigEmax", type = "general")
  test_perm <- MCTtest(doses2, drEst2, S = vc2, models = modlist, type = "general", df = Inf)
  # we don't compare stuff we want to be different
  attr(fit_orig, "data") <- attr(fit_perm, "data") <- NULL
  attr(fit_orig, "doseRespNam") <- attr(fit_perm, "doseRespNam") <- NULL
  expect_equal(fit_orig, fit_perm)
  expect_equal(test_orig, test_perm)
})

test_that("a non-psd covariance matrix produces NA p-values", {
  doses <- c(0,10,20,40)
  models <- Mods(doses=doses, emax=0.15, sigEmax=c(0.05, 5), linear=NULL,
               exponential=0.2,quadratic=-0.6,betaMod=c(0.05, 4))
  data.sim <- structure(list(X = structure(1:4, .Label = c("0", "10", "20", "40"), class = "factor"),
                             dose = c(0L, 10L, 20L, 40L),
                             Estimate = c(0.266942236, 3.792703657, 14.69084734, 17.71179102),
                             Cov1 = c(3.685607913, 0.595285049, 0.651289991, 0.742901538),
                             Cov2 = c(0.595285049, 3.31255546, 0.47843908, 0.545737127),
                             Cov3 = c(0.651289991, 0.47843908, 3.398708786, 0.597080557),
                             Cov4 = c(0.742901538, 0.545737127, 0.597080557, 3.556324648)),
                        class = "data.frame", row.names = c(NA, -4L))
  mu <- data.sim[,3]
  S <- data.matrix(data.sim[,4:7],rownames.force = NA)
  tst <- MCTtest(dose=doses,resp=mu,models = models,S=S,type="general")
  # p-value of linear model should be NA
  # FIXME: eigen(S) == 5.327424 2.875258 2.875258 2.875258 => S is psd!
  expect_equal(attr(tst$tStat, "pVal")[3], NA)
})
