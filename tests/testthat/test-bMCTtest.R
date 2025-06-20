# TODO:
# * maybe define common candidate models outside of test_that() calls?
# * how do we check for equal p-values (calculated with MC algorighm)?
# * pull shared code out of test_that() calls

source("generate_test_datasets.R")

require_rbest <- function() {
  if (!require("RBesT")) {
    skip("RBesT package not available")
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
tstat.bMCTtest <- function(obj) {
  # drop the pVal attribute of obj$tStat
  as.numeric(obj$tStat)
}
tstat.glht <- function(obj) {
  unname(summary(obj)$test$tstat)
}
critVal2 <- function(obj) {
  UseMethod("critVal")
}
critVal.MCTtest <- function(obj) {
  as.numeric(obj$critVal)
}
critVal.bMCTtest <- function(obj) {
  as.numeric(obj$critVal)
}
pVal.bMCTtest <- function(obj) {
  as.numeric(attr(obj$tStat, "pVal"))
}

twoarm_rbest <- function(dat, prior1, prior2){
  
  mod <- lm(y ~ as.factor(x) - 1, data = dat)
  mu1 <- coef(mod)[1]
  mu2 <- coef(mod)[2]
  S <- vcov(mod)
  post1 <- postmix(prior1, m = mu1, se = sqrt(S[1,1]))
  post2 <- postmix(prior2, m = mu2, se = sqrt(S[2,2]))
  pmixdiff(post1, post2, 0)
}


test_that("bMCTtest with uninformative prior produces same results as frequentist MCP-Mod", {
  require_rbest()
  set.seed(23)
  dd <- getDFdataSet_testsMCT()
  mD <- max(dd$x)
  nD <- length(unique(dd$x))
  lg1 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.9), "logistic")
  lg2 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.5), "logistic")
  expo <- guesst(c(0.9*mD), c(0.7), "exponential", Maxd=mD)
  quad <- guesst(c(0.6*mD), c(1), "quadratic")
  noninf_prior <- mixnorm(c(1, 0, 10000))
  prior <- vector("list", nD)
  for(i in 1:nD)
    prior[[i]] <- noninf_prior
  
  models <- Mods(linlog = NULL, logistic = rbind(lg1, lg2),
                 exponential = expo, quadratic = quad,
                 doses = dd$x, addArgs=list(off = 0.2*max(dd$x)))
  mcp_freq <- MCTtest(x,y , dd, models = models, df = Inf, critV = TRUE)
  mcp_bayes <- bMCTtest(x,y, dd, models=models, prior = prior)
  expect_equal(tstat(mcp_freq), tstat(mcp_bayes), tolerance = 0.001)
  expect_equal(1-pnorm(critVal2(mcp_freq)), critVal2(mcp_bayes), tolerance = 0.001)
})

test_that("bMCTtest works with contrast matrix handed over and produces same results", {
  require_rbest()
  set.seed(23)
  dd <- getDFdataSet_testsMCT()
  mD <- max(dd$x)
  nD <- length(unique(dd$x))
  lg1 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.9), "logistic")
  lg2 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.5), "logistic")
  expo <- guesst(c(0.9*mD), c(0.7), "exponential", Maxd=mD)
  quad <- guesst(c(0.6*mD), c(1), "quadratic")
  prior <- vector("list", nD)
  for(i in 1:nD)
    prior[[i]] <- mixnorm(c(1, 0, 10000))
  models <- Mods(linlog = NULL, logistic = rbind(lg1, lg2),
                 exponential = expo, quadratic = quad,
                 doses = dd$x, addArgs=list(off = 0.2*max(dd$x)))
  mcp_freq <- MCTtest(x,y , dd, models = models, df = Inf, critV = TRUE)
  mcp_bayes <- bMCTtest(x,y, dd, models=models, prior = prior, contMat = mcp_freq$contMat)
  expect_equal(tstat(mcp_freq), tstat(mcp_bayes), tolerance = 0.001)
  expect_equal(1-pnorm(critVal2(mcp_freq)), critVal2(mcp_bayes), tolerance = 0.001)
})

test_that("bMCTtest works with binary data (1)", {
  require_rbest()
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
  prior <- vector("list", length(dose))
  for(i in 1:length(unique(dd$x)))
    prior[[i]] <- mixnorm(c(1, 0, 10000))
  mcp_freq <- MCTtest(dose, dePar, S=vCov, models=models, type = "general", df = Inf, critV = TRUE)
  mcp_bayes <- bMCTtest(dose, dePar, S=vCov, models=models, prior = prior, type = "general")
  expect_equal(tstat(mcp_freq), tstat(mcp_bayes), tolerance = 0.001)
  expect_equal(1-pnorm(critVal2(mcp_freq)), critVal2(mcp_bayes), tolerance = 0.001)
})

test_that("MCTtest works with binary data (2)", {
  require_rbest()
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
  prior <- vector("list", length(dose))
  for(i in 1:length(dose))
    prior[[i]] <- mixnorm(c(1, 0, 10000))
  mcp_freq <- MCTtest(dose, dePar, S=vCov, models=models, type = "general", df = Inf, critV = TRUE)
  mcp_bayes <- bMCTtest(dose, dePar, S=vCov, models=models, prior = prior, type = "general")
  expect_equal(tstat(mcp_freq), tstat(mcp_bayes), tolerance = 0.001)
  expect_equal(1-pnorm(critVal2(mcp_freq)), critVal2(mcp_bayes), tolerance = 0.001)
})

test_that("MCTtest works with binary data (3)", {
  require_rbest()
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
  prior <- vector("list", length(dose))
  for(i in 1:length(dose))
    prior[[i]] <- mixnorm(c(1, 0, 10000))
  mcp_freq <- MCTtest(dose, dePar, S=vCov, models=models, type = "general", df = Inf, critV = TRUE)
  mcp_bayes <- bMCTtest(dose, dePar, S=vCov, models=models, prior = prior, type = "general")
  expect_equal(tstat(mcp_freq), tstat(mcp_bayes), tolerance = 0.001)
  expect_equal(1-pnorm(critVal2(mcp_freq)), critVal2(mcp_bayes), tolerance = 0.001)
})

test_that("a one-dimensional test works", {
  require_rbest()
  set.seed(1)
  dd <- getDFdataSet.bin()
  model <- Mods(linear = NULL, doses=sort(unique(dd$x)), addArgs=list(scal = 1.2*max(dd$x)))
  logReg <- glm(y~as.factor(x)-1, family=binomial, data=dd, weights = n)
  dePar <- coef(logReg)
  vCov <- vcov(logReg)
  dose <- sort(unique(dd$x))
  prior <- vector("list", length(dose))
  for(i in 1:length(dose))
    prior[[i]] <- mixnorm(c(1, 0, 10000))
  mcp_freq <- expect_warning(MCTtest(dose, dePar, S=vCov, models=model, type = "general", critV = TRUE, df=Inf),
                 "univariate: using pnorm")
  mcp_bayes <- bMCTtest(dose, dePar, S=vCov, models=model, type = "general", prior = prior)

  expect_equal(tstat(mcp_freq), tstat(mcp_bayes), tolerance = 0.001)
  expect_equal(1-pnorm(critVal2(mcp_freq)), critVal2(mcp_bayes), tolerance = 0.001)
})

test_that("unordered values in MCTtest work (unadjusted scale)", {
  require_rbest()
  data(IBScovars)
  modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                  linInt = c(0, 1, 1, 1), doses = c(0, 1, 2, 3, 4))
  ancMod <- lm(resp~factor(dose)-1, data=IBScovars)
  drEst <- coef(ancMod)
  vc <- vcov(ancMod)
  doses <- 0:4
  noninf_prior <- mixnorm(c(1, 0, 10000))
  prior <- vector("list", length(doses))
  for(i in 1:length(doses))
    prior[[i]] <- mixnorm(c(1, 0, 10000))
  bnds <- defBnds(max(doses))$sigEmax
  test_orig <- bMCTtest(doses, drEst, S = vc, models = modlist, type = "general", prior = prior)
  ord <- c(3,4,1,2,5)
  drEst2 <- drEst[ord]
  vc2 <- vc[ord,ord]
  doses2 <- doses[ord]
  test_perm <- bMCTtest(doses2, drEst2, S = vc2, models = modlist, type = "general", prior = prior)
  expect_equal(tstat(test_orig), tstat(test_perm))
  expect_equal(critVal2(test_orig), critVal2(test_perm), tolerance = 0.001)
})

test_that("bMCTtest gives same results as RBesT two-sample analysis with non-informative prior", {
  require_rbest()
  set.seed(23)
  dd <- getDFdataSet_testsMCT()
  ## only keep the highest and lowest dose
  dd <- dd[dd$x %in% range(dd$x), ]
  mD <- max(dd$x)
  model <- Mods(linear = NULL, doses=sort(unique(dd$x)))
  prior <- list(mixnorm(c(1, 0, 1000)), mixnorm(c(1, 0, 1000)))
  twoarm <- twoarm_rbest(dd, prior[[1]], prior[[2]])
  mcp_bayes <- bMCTtest(x,y, dd, models=model, prior = prior)
  expect_equal(twoarm, pVal.bMCTtest(mcp_bayes))
})


test_that("bMCTtest gives same results as RBesT two-sample analysis with informative prior for control", {
  require_rbest()
  set.seed(23)
  dd <- getDFdataSet_testsMCT()
  ## only keep the highest and lowest dose
  dd <- dd[dd$x %in% range(dd$x), ]
  mD <- max(dd$x)
  model <- Mods(linear = NULL, doses=sort(unique(dd$x)))
  noninf_prior <- mixnorm(c(1, 0, 1000))
  inf_prior <- mixnorm(c(1, 0, 1))
  prior <- list(inf_prior, noninf_prior)
  twoarm <- twoarm_rbest(dd, prior[[1]], prior[[2]])
  mcp_bayes <- bMCTtest(x,y, dd, models=model, prior = prior)
  expect_equal(twoarm, pVal.bMCTtest(mcp_bayes))
})

test_that("bMCTtest gives same results as RBesT two-sample analysis with informative prior for both arms", {
  require_rbest()
  set.seed(24)
  dd <- getDFdataSet_testsMCT()
  ## only keep the highest and lowest dose
  dd <- dd[dd$x %in% range(dd$x), ]
  mD <- max(dd$x)
  model <- Mods(linear = NULL, doses=sort(unique(dd$x)))
  inf_prior_cont <- mixnorm(c(0.8, 0, 1), c(0.1, 1, 2), c(0.1, -1, 2))
  inf_prior_trt <- mixnorm(c(0.5, 1, 1), c(0.3, 0.8, 2), c(0.2, 1.5, 2))
  prior <- list(inf_prior_cont, inf_prior_trt)
  twoarm <- twoarm_rbest(dd, prior[[1]], prior[[2]])
  mcp_bayes <- bMCTtest(x,y, dd, models=model, prior = prior)
  expect_equal(twoarm, pVal.bMCTtest(mcp_bayes))
})

test_that("Error message for incorrect prior arguments", {
  data(biom)
  ## define shapes for which to calculate optimal contrasts
  doses <- c(0, 0.05, 0.2, 0.6, 1)
  modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                  linInt = c(0, 1, 1, 1), doses = doses)
  ## specify an informative prior for placebo, weakly informative for other arms
  plc_prior <- mixnorm(inf = c(0.8, 0.4, 0.1), rob = c(0.2, 0.4, 10))
  vague_prior <- mixnorm(c(1, 0, 10))
  ## one component of the list corresponds to each dose
  prior1 <- list(plc_prior, vague_prior)
  prior2 <- list(plc_prior, "foo", "foo", "foo", "foo")
  expect_error(bMCTtest(dose, resp, biom, models=modlist, prior = prior1),
               "Dose and prior have non-conforming size")
  expect_error(bMCTtest(dose, resp, biom, models=modlist, prior = prior2),
               "priors need to be of class normMix")
})
