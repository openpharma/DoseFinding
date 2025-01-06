# create contrast matrix

doses <- c(0, 10, 25, 50, 100, 150)
fmodels <- Mods(linear = NULL, emax = 25, logistic = c(50, 10.88111),
                exponential = 85, betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)),
                doses = doses, addArgs = list(scal = 200), placEff = 0, maxEff = 0.4)
contMat <- optContr(fmodels, w = 1)


########################################################
## General tests: does powMCT work in various scenarios
########################################################
test_that("powMCT computes with default values", {
  power <- powMCT(contMat, altModels = fmodels, n = 50,  sigma = 1)
  expect_true(is.numeric(power))
  expect_true(all(power > 0 & power <= 1))
})

test_that("powMCT works with specified covariance matrix S", {
  doses <- c(0, 10, 25, 50, 100, 150)
  S <- 1^2 / 50 * diag(length(doses))
  power <- powMCT(contMat, altModels = fmodels, S = S, df = 50 * length(doses) - length(doses), alpha = 0.05)
  expect_true(is.numeric(power))
  expect_true(all(power > 0 & power <= 1))
})

test_that("powMCT works with placebo adjusted estimates", {
  doses <- c(10, 25, 50, 100, 150)
  fmodels <- Mods(linear = NULL, emax = 25, logistic = c(50, 10.88111),
                  exponential = 85, betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)),
                  doses = c(0, doses), addArgs = list(scal = 200), placEff = 0, maxEff = 0.4)
  contMat <- optContr(fmodels, doses = doses, w = 1, placAdj = TRUE)
  power <- powMCT(contMat, altModels = fmodels, n = 50,  sigma = 1, placAdj = TRUE)
  expect_true(is.numeric(power))
  expect_true(all(power > 0 & power <= 1))
})

test_that("powMCT works with two-sided alternative", {
  power <- powMCT(contMat, altModels = fmodels, n = 50, alpha = 0.05, sigma = 1, alternative = "two.sided")
  expect_true(is.numeric(power))
  expect_true(all(power > 0 & power <= 1))
})

########################################################
## Error testing
########################################################

test_that("powMCT errors when required arguments are missing", {
  expect_error(powMCT(contMat, altModels = fmodels), "Either S or both n and sigma need to be specified")
  expect_error(powMCT(contMat, n = 50, sigma = 1), "altModels argument needs to be specified")
})

test_that("powMCT detects invalid inputs", {
  expect_error(powMCT(list(contMat), altModels = fmodels, n = 50,  sigma = 1), "contMat needs to be a matrix")
  expect_error(powMCT(contMat, altModels = fmodels, n = rep(50, 2),  sigma = 1), "n needs to be of length nrow")
  expect_error(powMCT(contMat, altModels = fmodels, n = 1,  sigma = 1), "cannot compute power: specified \"n\" and dose vector result in df = 0")
  
  S <- matrix(1, nrow = 6, ncol = 6)
  badS1 <- matrix(1, nrow = 3, ncol = 3)
  badS2 <- matrix(1, nrow = 5, ncol = 6)
  expect_error(powMCT(contMat, altModels = fmodels, S = S), "need to specify degrees of freedom in \"df\", when specifying \"S\"")
  expect_error(powMCT(contMat, altModels = fmodels, S = S, n = 50), "Need to specify either \"S\" or both \"n\" and \"sigma\"")
  expect_error(powMCT(contMat, altModels = fmodels, S = badS1, df = 45), "S needs to have as many rows&cols as there are doses")
  expect_error(powMCT(contMat, altModels = fmodels, S = badS2, df = 45), "S needs to be a square matrix")
  
  fmodels2 <- Mods(linear = NULL, emax = 25, logistic = c(50, 10.88111),
                  exponential = 85, betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)),
                  doses = c(0, 10, 25, 50, 100), placEff = 0, maxEff = 0.4)
  expect_error(powMCT(contMat, altModels = fmodels2, n = 50,  sigma = 1), "Incompatible contMat and muMat")
  
})

########################################################
## Test power calculations
########################################################

test_that("powMCT gives same result as power.t.test", {
  doses <- c(0, 1)
  fmodels <- Mods(linear = NULL, doses = doses, placEff = 0, maxEff = 0.4)
  contMat <- optContr(fmodels, w = 1)
  power1 <- powMCT(contMat, altModels = fmodels, n = 50,  sigma = 1)
  power2 <- powMCT(contMat, altModels = fmodels, n = 50,  sigma = 1, alpha = 0.05, alternative = "two.sided")
  expect_equal(as.numeric(power1), power.t.test(n = 50, delta = 0.4, sd = 1, sig.level = 0.025, alternative = "one.sided")$power, tolerance = 0.001)
  expect_equal(as.numeric(power2), power.t.test(n = 50, delta = 0.4, sd = 1, sig.level = 0.05, alternative = "two.sided")$power, tolerance = 0.001)

})

## compare power to externally calculated power values
test_that("powMCT calculates power correctly", {
  doses <- c(0, 25, 100, 300)
  fmodels <- Mods(emax = 45, sigEmax = c(100,3 ), logistic = c(45, 15),
                  exponential = 60, quadratic = -0.0022,
                  doses = doses, placEff = 0, maxEff = 1)
  contMat <- optContr(fmodels, w = 1)
  power <- powMCT(contMat, altModels = fmodels, n = 100,  sigma = 3)
  expect_equal(as.numeric(power), c(0.672, 0.754, 0.817, 0.757, 0.635), tolerance = 0.01)
})



