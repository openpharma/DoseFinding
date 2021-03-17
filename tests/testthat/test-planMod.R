context("planning models")

# TODO
# * various FIXMEs in test #2
# * what do we want to do with tests #3-5 (mostly plots)
# * test #4 crashes in planMod

# test 1: "validation" using fitMod and vcov.DRMod and predict.DRMod
test_that("getPredVar gives the same results as predict.DRMod", {
  n <- c(100,50,50,50,100)
  doses <- c(0,10,20,40,50)
  cf <- c(0,1,10)
  V <- DoseFinding:::aprCov(doses, "emax", cf, S=diag(1/n))
  pv1 <- DoseFinding:::getPredVar("emax", cf, V=V, pDose=50)
  # now validation using the formulas in fitMod
  doseVec <- rep(doses, n)
  respVec <- rnorm(length(doseVec))
  dd <- fitMod(doseVec, respVec, model="emax")
  # now change to achieve desired values
  dd$coefs <- cf
  dd$df <- dd$RSS <- sum(n)
  pv2 <- predict(dd, predType = "effect-curve", doseSeq=50, se.fit=TRUE)$se.fit^2
  expect_equal(pv1, pv2)
})

# test 2: "validation" using simulation
test_that("get_{TD,ED,Pred}Var gives the same result as a simulation", {
  skip_on_cran()
  # select very large sample size (to validate asymptotics)
  n <- c(100000, 50000, 50000, 50000, 100000) # FIXME: which one?
  n <- c(100, 50, 50, 50, 100)
  doses <- c(0,10,20,40,50)
  cf <- c(0,0.4,10)
  Delta <- 0.2

  mm <- Mods(emax=cf[3], doses=doses, placEff=0, maxEff=emax(50,cf[1],cf[2],cf[3]))
  true_values <- unname(c(ED(mm, p=0.5),
                          TD(mm, Delta=Delta),
                          emax(50, cf[1], cf[2], cf[3])))
  V <- DoseFinding:::aprCov(doses, "emax", cf, S=diag(0.3^2/n))
  true_variances <- unname(c(
    DoseFinding:::getEDVar("emax", cf, V=V, p=0.5, maxD=50, scale = "unrestricted"),
    DoseFinding:::getTDVar("emax", cf, V=V, Delta=Delta, scale = "unrestricted"),
    DoseFinding:::getPredVar("emax", cf, V=V, pDose=50)))

  # simulation
  mn <- emax(doses, cf[1], cf[2], cf[3])
  doseVec <- rep(doses, n)
  mnVec <- rep(mn, n)
  one_sim <- function() {
    respVec <- mnVec + rnorm(length(mnVec),0,0.3)
    ff <- fitMod(doseVec, respVec, model="emax", bnds = c(0.05, 75))
    return(c(ed = ED(ff, p=0.5),
             td = TD(ff, Delta = Delta),
             pl = predict(ff, doseSeq=50, predType = "effect-curve")))
  }
  sim <- replicate(100, one_sim()) # FIXME: originally 10000
  expect_equal(unname(rowMeans(sim)), true_values) # FIXME: tolerance?
  expect_equal(unname(apply(sim, 1, var)), true_variances)

  ## FIXME: what do we do with this?
  ## tdvart <- DoseFinding:::getTDVar("emax", cf, V=V,scale = "log", Delta=Delta)
  ## hist(td[td < 100], freq=FALSE, breaks = 101)
  ## curve(dlnorm(x, log(tdt), sqrt(tdvart)), add=TRUE)

  edt7 <- ED(mm, p=0.7)
  edt3 <- ED(mm, p=0.3)
  expect_equal(mean(sim[1,] < edt7 & sim[1,] > edt3),
               unname(pnorm(edt7, edt, sqrt(true_variances[1]))-pnorm(edt3, edt, true_variances[1])))
})


# test 5: C's example
test_that("negative values for Delta lead to an error", {
  models <- Mods(linear = NULL, linlog = NULL, emax = c(8, 10),
                 sigEmax = c(10, 2),
                 doses = c(0, 10,20, 50, 100),
                 placEff=0, maxEff=-2)
  pObj <- planMod("sigEmax",models,n=100,sigma=1.2,
                  simulation=TRUE,nSim=100, cores=4)
  expect_error(summary(pObj,Delta=-1.1),
               "\"Delta\" needs to be > 0")
})
