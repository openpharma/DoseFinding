# create contrast matrix

doses <- c(0, 10, 25, 50, 100, 150)
fmodels <- Mods(linear = NULL, emax = 25, logistic = c(50, 10.88111),
                exponential = 85, betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)),
                doses = doses, addArgs = list(scal = 200), placEff = 0, maxEff = 0.4)
contMat <- optContr(fmodels, w = 1)

tFunc <- function(n) {
  powVals <- powMCT(contMat, altModels = fmodels, n = n, sigma = 1)
  mean(powVals)
}


#########################################################
## General tests: does sampSize work as expected
#########################################################

test_that("sampSize and sampSizeMCT work correctly", {
  skip_on_ci()

  result1 <- sampSize(upperN = 80, targFunc = tFunc, target = 0.8, alRatio = rep(1, 6))
  result2 <- sampSizeMCT(upperN = 80, power = 0.8, alRatio = rep(1, 6), contMat = contMat, 
                         altModels = fmodels, sigma = 1)
  expect_true(is.list(result1))
  expect_true(all(result1$samp.size > 0))
  expect_true(result1$target > 0)
  expect_true(is.list(result2))
  expect_true(all(result2$samp.size > 0))
  expect_true(result2$target > 0)
  expect_equal(result1$samp.size, result2$samp.size, tolerance = 1)
})

test_that("sampSize and sampSizeMCT work correctly with Ntype = total", {
  skip_on_ci()

  result1 <- sampSize(upperN = 80, targFunc = tFunc, target = 0.8, alRatio = rep(1, 6), Ntype = "total")
  result2 <- sampSizeMCT(upperN = 80, power = 0.8, alRatio = rep(1, 6), contMat = contMat, 
                         altModels = fmodels, sigma = 1, Ntype = "total")
  expect_true(is.list(result1))
  expect_true(all(result1$samp.size > 0))
  expect_true(result1$target > 0)
  expect_true(is.list(result2))
  expect_true(all(result2$samp.size > 0))
  expect_true(result2$target > 0)
  expect_equal(result1$samp.size, result2$samp.size, tolerance = 6)
})

test_that("sampSize and sampSizeMCT work correctly with sumFct = min", {
  skip_on_ci()
  
  result1 <- sampSize(upperN = 80, targFunc = tFunc, target = 0.8, alRatio = rep(1, 6), Ntype = "total")
  result2 <- sampSizeMCT(upperN = 80, power = 0.8, alRatio = rep(1, 6), contMat = contMat, 
                         altModels = fmodels, sigma = 1, Ntype = "total", sumFct = min)
  expect_true(is.list(result1))
  expect_true(all(result1$samp.size > 0))
  expect_true(result1$target > 0)
  expect_true(is.list(result2))
  expect_true(all(result2$samp.size > 0))
  expect_true(result2$target > 0)
  expect_equal(result1$samp.size, result2$samp.size, tolerance = 6)
})

test_that("sampSizeMCT handles S matrix correctly", {
  S <- 6 * diag(6)
  result <- sampSizeMCT(upperN = 500, contMat = contMat, S = S, altModels = fmodels,
                        power = 0.8, alRatio = rep(1, 6),  Ntype = "total")
  
  expect_true(is.list(result))
  expect_true(all(result$samp.size > 0))
  expect_true(result$target > 0)
})

test_that("print.sampSize prints correct output", {
  sSize_obj <- sampSize(upperN = 80, targFunc = tFunc, target = 0.8, alRatio = rep(1, 6), Ntype = "total")
  
  expect_output(print(sSize_obj), "Sample size calculation\n\n", fixed = TRUE)
  expect_output(print(sSize_obj), "alRatio: ", fixed = TRUE)
  expect_output(print(sSize_obj), paste("Total sample size:", sum(sSize_obj$samp.size)), fixed = TRUE)
  expect_output(print(sSize_obj), paste("Sample size per arm:", paste(sSize_obj$samp.size, collapse = " ")), fixed = TRUE)
  expect_output(print(sSize_obj), paste("targFunc:", sSize_obj$target), fixed = TRUE)
})


########################################################
## Error testing
########################################################

test_that("sampSize errors with invalid inputs", {
  expect_error(sampSize(upperN = 80, targFunc = NULL, target = 0.8, alRatio = rep(1, 6)),
               "targFunc")
  expect_error(sampSize(upperN = 80, targFunc = function(n) { n }, target = 0.8),
               "allocation ratios need to be specified")
  expect_error(sampSize(upperN = 80, targFunc = function(n) { n }, target = 0.8, alRatio = c(1, -1, 1)),
               "all entries of alRatio need to be positive")
})

test_that("sampSizeMCT errors with invalid inputs", {
  expect_error(sampSizeMCT(upperN = 80, contMat = contMat, altModels = fmodels,
                           power = 0.8, alRatio = rep(1, 6), alpha = 0.025, placAdj = TRUE),
               "placAdj needs to be FALSE")
  expect_error(sampSizeMCT(upperN = 80, contMat = contMat, altModels = fmodels,
                           power = 0.8, alRatio = rep(1, 6), sigma = 1, n = 50),
               "n is not allowed to be specified")
  expect_error(sampSizeMCT(upperN = 80, contMat = contMat, altModels = fmodels,
                           power = 0.8, alRatio = rep(1, 6)),
               "need sigma if S is not specified")
})


########################################################
## Compare results to powMCT
########################################################

# Tests for sampSizeMCT function
test_that("sampSizeMCT results are consistent with powMCT", {
  result <- sampSizeMCT(upperN = 80, contMat = contMat, sigma = 1, altModels = fmodels,
                        power = 0.8, alRatio = rep(1, 6), alpha = 0.025)
  power <- powMCT(contMat, altModels = fmodels, sigma = 1, n = result$samp.size[1])
  
  expect_equal(result$target, mean(power), tolerance = 0.01)
})

########################################################
## Testing targN and powN
########################################################
tFunc <- function(n) {
  powVals <- powMCT(contMat, altModels = fmodels, n = n, sigma = 1)
  powVals
}

test_that("targN calculates target function values correctly", {

  # Perform targN calculations
  result <- targN(upperN = 100, lowerN = 10, step = 10, targFunc = tFunc, alRatio = rep(1, 6))
  power <- powMCT(contMat, altModels = fmodels, n = 50, sigma = 1)
  power <- c(power, min(power), mean(power), max(power))
  expect_true(is.matrix(result))
  expect_true(nrow(result) > 0)
  expect_true(ncol(result) > 0)
  expect_true(all(result > 0 & result <= 1))
  expect_true(all(c("min", "mean", "max") %in% colnames(result)))
  expect_equal(as.numeric(power), as.numeric(result[5, ]), tolerance = 0.01)
})

test_that("powN calculates power correctly", {
  
  # Perform targN calculations
  result <- powN(upperN = 100, lowerN = 10, step = 10, alRatio = rep(1, 6), contMat = contMat, sigma = 1, altModels = fmodels)
  power <- powMCT(contMat, altModels = fmodels, n = 50, sigma = 1)
  power <- c(power, min(power), mean(power), max(power))
  expect_true(is.matrix(result))
  expect_true(nrow(result) > 0)
  expect_true(ncol(result) > 0)
  expect_true(all(result > 0 & result <= 1))
  expect_true(all(c("min", "mean", "max") %in% colnames(result)))
  expect_equal(as.numeric(power), as.numeric(result[5, ]), tolerance = 0.01)
})


test_that("targN errors with invalid inputs", {

  expect_error(targN(upperN = 100, lowerN = 10, step = 10, targFunc = tFunc),
               "allocation ratios need to be specified")
  expect_error(targN(upperN = 100, lowerN = 10, step = 10, targFunc = tFunc, alRatio = c(1, 0, 1)),
               "all entries of alRatio need to be positive")
  expect_error(targN(upperN = 100, lowerN = 10, step = 10, targFunc = tFunc, alRatio = rep(1, 6), sumFct = mean),
               "sumFct needs to be a character vector")
})


test_that("powN errors with invalid inputs", {
  
  expect_error(powN(upperN = 100, lowerN = 10, step = 10, alRatio = rep(1, 6), contMat = contMat, sigma = 1, altModels = fmodels, placAdj = TRUE),
               "placAdj needs to be FALSE for powN")
  expect_error(powN(upperN = 100, lowerN = 10, step = 10, alRatio = rep(1, 6), contMat = contMat, sigma = 1, n = 50, altModels = fmodels),
               "n is not allowed to be specified for sample size calculation")
  expect_error(powN(upperN = 100, lowerN = 10, step = 10, alRatio = rep(1, 6), contMat = contMat, altModels = fmodels),
               "need sigma if S is not specified")
})

test_that("powN handles S matrix correctly", {
  S <- 6 * diag(6)
  result <- powN(upperN = 100, lowerN = 10, step = 10,  contMat = contMat, S = S, altModels = fmodels,
                        alRatio = rep(1, 6), Ntype = "total")
  expect_true(is.matrix(result))
  expect_true(nrow(result) > 0)
  expect_true(ncol(result) > 0)
  expect_true(all(result > 0 & result <= 1))
  expect_true(all(c("min", "mean", "max") %in% colnames(result)))
})

# Tests for plot.targN function
tn <- targN(upperN = 100, lowerN = 10, step = 10, targFunc = tFunc, alRatio = rep(1, 6))

test_that("plot.targN works as expected", {
  expect_error(plot(tn), NA)
  # Test with superpose = TRUE
  expect_error(plot(tn, superpose = TRUE), NA)
  # Test with superpose = FALSE
  expect_error(plot(tn, superpose = FALSE), NA)
  # Test with line.at specified
  expect_error(plot(tn, line.at = 0.8), NA)
  # Test with line.at as NULL
  expect_error(plot(tn, line.at = NULL), NA)
  # Test with custom xlab and ylab
  expect_error(plot(tn, xlab = "Sample Size", ylab = "Power"), NA)
  # Test with default xlab and ylab
  expect_error(plot(tn, xlab = NULL, ylab = NULL), NA)
})
