context("optimal designs")

# TODO
# * mixed Paper p. 1233, l. 2 (note the off and probably also the scal
#   parameter were treated as unknown in this example in the paper, hence the
#   results need not be consistent with paper)
#
# * everything from the "some other examples" section
#
# * optimizer = "exact" and "solnp" (weights vary by up to ~4 percentage points)
#
# * Example from Padmanabhan and Dragalin, Biometrical Journal 52 (2010) p. 836-852
#
# * optimial design logistic regression; compare this to Atkinson et al. (2007), p. 400

## Recreate examples from this this article
##
## @Article{dette2008,
##   author = 	 {Dette, Holger and Bretz, Frank and Pepelyshev, Andrey and Pinheiro, José},
##   title = 	 {Optimal Designs for Dose-Finding Studies},
##   journaltitle = {Journal of the American Statistical Association},
##   year = 	 2008,
##   volume = 	 103,
##   issue = 	 483,
##   pages = 	 {1225-1237},
##   doi = 	 {10.1198/016214508000000427}}

# Note: expect_equal(..., tolerance = 1e-3) in most instances, because the
# published results have three or four decimal places

test_that("the emax model (table 2, line 5) gives the same results", {
  fMod <- Mods(emax = 25, doses = c(0,150), placEff=0, maxEff=0.4)
  fMod$emax[2] <- 0.6666667
  doses <- c(0, 18.75, 150)
  probs <- 1
  deswgts1 <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD",
                        optimizer="Nelder-Mead")
  deswgts2 <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD",
                        optimizer="nlminb")
  expect_equal(deswgts1$design, deswgts2$design, tolerance = 1e-4)
  expect_equal(deswgts1$design, c(0.442, 0.5, 0.058), tolerance = 1e-3)
  ## efficiency compared to standard design (last column)
  crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                  Delta=0.2, designCrit = "TD")
  expect_equal(exp(deswgts1$crit - crt), 0.5099, tolerance = 1e-4)
})

test_that("the emax model (table 2, line 2) gives the same results", {
  fMod <- Mods(emax = 25, doses = c(0,150), placEff=0, maxEff=0.4)
  doses <- c(0, 18.75, 150)
  probs <- 1
  deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD")
  expect_equal(deswgts$design, c(0.5, 0.5, 0), tolerance = 1e-3)
})

test_that("the exponential model (table 3, line 2) gives the same results", {
  fMod <- Mods(exponential=85, doses = c(0, 150), placEff=0, maxEff=0.4)
  doses <- c(0, 50, 104.52, 150)
  probs <- 1
  deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD",
                       optimizer="Nelder-Mead")
  expect_equal(deswgts$design, c(0.5, 0, 0.5, 0), tolerance = 1e-3)
  # efficiency compared to standard design (last column)
  crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                  Delta=0.2, designCrit = "TD")
  expect_equal(exp(deswgts$crit - crt), 0.4286, tolerance = 1e-4)
})

test_that("the exponential model (table 3, line 1) gives the same results", {
  fMod <- Mods(exponential=65, doses=c(0, 150), placEff=0, maxEff=0.4)
  fMod$exponential[2] <- 0.08264711
  doses <- c(0, 101.57, 150)
  probs <- 1
  deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD")
  expect_equal(deswgts$design, c(0.440, 0.5, 0.060), tolerance = 1e-3)
})

test_that("the logistic model (table 4, line 7) gives the same results", {
  fMod <- Mods(logistic=c(50, 10.881), doses = c(0, 150), placEff=0, maxEff=0.4)
  doses <- c(0, 37.29, 64.44, 150)
  probs <- 1
  deswgts <- optDesign(fMod, probs, doses, Delta=0.05, designCrit = "TD")
  expect_equal(deswgts$design, c(0.401, 0.453, 0.099, 0.047), tolerance = 1e-3)
  ## efficiency compared to standard design (last column)
  crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                  Delta=0.05, designCrit = "TD")
  expect_equal(exp(deswgts$crit - crt), 0.1853, tolerance = 1e-4)
})

test_that("the logistic model (table 4, line 1) gives the same results", {
  fMod <- Mods(logistic=c(50, 10.881), doses = c(0, 150), placEff=0, maxEff=0.4)
  doses <- c(0, 50.22)
  probs <- 1
  deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD")
  expect_equal(deswgts$design, c(0.5, 0.5))
})

test_that("the beta model (table 5, line 5) gives the same results", {
  fMod <- Mods(betaMod = c(0.33, 2.31), doses = c(0,150), addArgs=list(scal=200),
               placEff=0, maxEff=0.4)
  doses <- c(0, 0.49, 25.2, 108.07, 150)
  probs <- 1
  deswgts <- optDesign(fMod, probs, doses, Delta=0.1,
                       control=list(maxit=1000), designCrit = "TD")
  expect_equal(deswgts$design, c(0.45, 0.48, 0.05, 0.02, 0), tolerance = 1e-2)
  ## efficiency compared to standard design (last column)
  crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                  Delta=0.1, designCrit = "TD")
  expect_equal(exp(deswgts$crit - crt), 0.130, tolerance = 1e-3)
})

test_that("the beta model (table 5, line 10) gives the same results", {
  fMod <- Mods(betaMod = c(1.39, 1.39), doses=c(0, 150), addArgs=list(scal=200),
               placEff=0, maxEff=0.4)
  doses <- c(0, 27, 94.89, 150)
  probs <- 1
  deswgts <- optDesign(fMod, probs, doses, Delta=0.1, designCrit = "TD")
  expect_equal(deswgts$design, c(0.45, 0.48, 0.05, 0.02), tolerance = 1e-2)
  ## efficiency compared to standard design (last column)
  crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                  Delta=0.1, designCrit = "TD")
  expect_equal(exp(deswgts$crit - crt), 0.501, tolerance = 1e-3)
})

test_that("the beta model (table 5, line 1) gives the same results", {
  fMod <- Mods(betaMod = c(0.23, 2.31), doses=c(0,150), addArgs=list(scal=200),
               placEff=0, maxEff=0.4)
  doses <- c(0, 0.35, 150)
  probs <- 1
  deswgts <- optDesign(fMod, probs, doses, Delta=0.2, designCrit = "TD")
  expect_equal(deswgts$design, c(0.5, 0.5, 0), tolerance = 1e-2)
  ## efficiency compared to standard design (last column)
  crt <- calcCrit(rep(1/6, 6), fMod, probs, c(0, 10, 25, 50, 100, 150),
                  Delta=0.2, designCrit = "TD")
  expect_equal(exp(deswgts$crit - crt), 0.056, tolerance = 1e-3)
})

test_that("standardized Dopt and Dopt&TD criteria work", {
  doses <- c(0, 62.5, 125, 250, 500)
  fMod1 <- Mods(sigEmax = rbind(c(25, 5), c(107.14, 2)), doses=doses, placEff=60, maxEff=280)
  fMod2 <- Mods(sigEmax = rbind(c(25, 5), c(107.14, 2)), linear = NULL,
                doses=doses, placEff=60, maxEff=280)
  w1 <- rep(0.5, 2)
  w2 <- rep(1/3, 3)
  ## des1 and des2 should be exactly the same
  des1 <- optDesign(fMod1, w1, doses, designCrit = "Dopt", standDopt = FALSE)
  des2 <- optDesign(fMod1, w1, doses, designCrit = "Dopt", standDopt = TRUE)
  expect_equal(des1$design, des2$design, tolerance =1e-6)
  ## des1 and des2 should be different (as linear and emax have different
  ## number of parameters)
  des1 <- optDesign(fMod2, w2, doses, designCrit = "Dopt", standDopt = FALSE,
                    optimizer = "solnp")
  des2 <- optDesign(fMod2, w2, doses, designCrit = "Dopt", standDopt = TRUE,
                    optimizer = "solnp")
  expect_false(all(des1$design == des2$design))
  ## same with Dopt&TD criterion: des1 and des2 will differ (due to different
  ## scaling of Dopt and TD criteria)
  des1 <- optDesign(fMod1, w1, doses, designCrit = "Dopt&TD",
                    Delta = 100, standDopt = FALSE,
                    optimizer = "solnp")
  des2 <- optDesign(fMod1, w1, doses, designCrit = "Dopt&TD",
                    Delta = 100, standDopt = TRUE,
                    optimizer = "solnp")
  expect_false(all(des1$design == des2$design))
})

## code using lower and upper bound (previous to version 0.9-6 this caused
## problems as the starting value for solnp rep(0.2, 5) was on the boundary,
## now a feasible starting values are used
test_that("feasible starting values are used when on boundary", {
  doses <- seq(0, 1, length=5)
  nold <- rep(0, times=5)
  lowbnd <- c(0.2,0.0,0.0,0.0,0.2)
  uppbnd <- c(1.0,0.3,1.0,1.0,1.0)
  trueModels <- Mods(linear=NULL, doses=doses, placEff = 0, maxEff = 1)
  des <- optDesign(models=trueModels, probs=1, doses=doses, designCrit="Dopt",
                   lowbnd=lowbnd,uppbnd=uppbnd)
  expect_equal(des$design, c(0.5, 0, 0, 0, 0.5))
})

test_that("there are no instabilities for numerical gradients", {
  mm <- Mods(betaMod=c(1.5,0.8), doses=seq(0,1,by=0.25), placEff=0, maxEff=1)
  des <- optDesign(mm, probs=1, designCrit="TD", Delta=0.5)
  expect_equal(des$design, c(0.4895, 0.3552, 0.1448, 0, 0.0105), tolerance = 1e-4)
})
