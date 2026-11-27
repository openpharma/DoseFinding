# Setting up some test data
doses <- c(0, 0.5, 1, 2, 4)
drFit <- c(1, 2, 3, 4, 5)  # Example response
S <- diag(5)  # Covariance matrix for simplicity

test_that("bFitMod errors with invalid inputs", {
  expect_error(bFitMod(dose = doses, resp = drFit, model = "invalidModel", S = S), 
               "invalid model selected")
  expect_error(bFitMod(dose = doses, resp = drFit[1:4], model = "linear", S = S), 
               "dose and resp need to be of the same size")
  expect_error(bFitMod(dose = doses, resp = drFit, model = "linear", S = diag(4)), 
               "S and dose have non-conforming size")
})


test_that("bFitMod correctly fits a 'linear' model with Bayes", {
  prior <- list(norm = c(0, 10), norm = c(0, 100))
  fit <- bFitMod(dose = doses, resp = drFit, model = "linear", S = S, 
                 type = "Bayes", nSim = 100, prior = prior)
  expect_s3_class(fit, "bFitMod")
  expect_equal(attr(fit, "model"), "linear")
  expect_equal(attr(fit, "type"), "Bayes")
  expect_true(!is.null(fit$samples))
})

test_that("bFitMod correctly fits a 'linear' model with bootstrap", {
  fit <- bFitMod(dose = doses, resp = drFit, model = "linear", S = S, 
                 type = "bootstrap", nSim = 100)
  expect_s3_class(fit, "bFitMod")
  expect_equal(attr(fit, "model"), "linear")
  expect_equal(attr(fit, "type"), "bootstrap")
  expect_true(!is.null(fit$samples))
})

test_that("print.bFitMod does not throw an error", {
  prior <- list(norm = c(0, 10), norm = c(0, 100))
  fit <- bFitMod(dose = doses, resp = drFit, model = "linear", S = S, 
                 type = "Bayes", nSim = 100, prior = prior)
  
  expect_output(print(fit), regexp = "Dose Response Model")
  expect_output(print(fit), regexp = "Summary of posterior draws")
})

test_that("bFitMod handles placebo adjustment appropriately", {
  prior <- list(norm = c(0, 10), norm = c(0, 100))
  expect_error(bFitMod(dose = doses, resp = drFit, model = "linlog", S = S, 
                       placAdj = TRUE, type = "Bayes", nSim = 100, prior = prior),
               "logistic and linlog models cannot be fitted to placebo adjusted data")
})

test_that("bFitMod correctly handles 'linInt' model", {
  fit <- bFitMod(dose = doses, resp = drFit, model = "linInt", S = S, 
                 type = "bootstrap", nSim = 100)
  expect_s3_class(fit, "bFitMod")
  expect_equal(attr(fit, "model"), "linInt")
  expect_true(!is.null(attr(fit, "nodes")))
  expect_true(!is.null(fit$samples))
})

test_that("bFitMod correctly handles additional arguments", {
  prior <- list(norm = c(0, 10), norm = c(0, 100), beta=c(0, 1.5, 0.45, 1.7), beta=c(0, 1.5, 0.45, 1.7))
  fit <- bFitMod(dose = doses, resp = drFit, model = "betaMod", S = S, 
                 type = "Bayes", nSim = 100, prior = prior, 
                 addArgs = list(scal = 1.2*max(doses)))
  expect_s3_class(fit, "bFitMod")
  expect_equal(attr(fit, "model"), "betaMod")
  expect_equal(attr(fit, "scal"), 1.2 * max(doses))
  expect_true(!is.null(fit$samples))
})

# Assuming the `biom` dataset is available in the environment for examples
data(biom)
anMod <- lm(resp ~ factor(dose) - 1, data = biom)
drFit <- coef(anMod)
S <- vcov(anMod)
dose <- sort(unique(biom$dose))

# Assuming normal priors for test example
prior <- list(norm = c(0, 10), norm = c(0, 100), beta = c(0, 1.5, 0.45, 1.7))

# Fit a model
gsample <- bFitMod(dose, drFit, S, model = "emax", start = c(0, 1, 0.1), nSim = 1000, prior = prior)

test_that("predict.bFitMod returns correct quantiles", {
  doseSeq <- c(0, 0.5, 1)
  pred <- predict(gsample, doseSeq = doseSeq)
  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), 5)  # Expecting rows for different quantiles
  expect_equal(length(unique(doseSeq)), ncol(pred))  # One column per dose in doseSeq
})

test_that("plot.bFitMod generates a plot", {
  expect_error(plot(gsample), NA)
  # Check for plotting is a little tricky, one way to check if some plot is generated
  expect_true(is.null(dev.list()) || length(dev.list()) > 0)
})

test_that("coef.bFitMod returns model coefficients", {
  coefs <- coef(gsample)
  expect_true(is.numeric(coefs))
  expect_equal(length(coefs), length(gsample$samples))
})

# To ensure the appropriate methods are defined, use methods(...) to list them:
test_that("appropriate methods for bFitMod are defined", {
  expect_true("predict.bFitMod" %in% methods("predict"))
  expect_true("plot.bFitMod" %in% methods("plot"))
  expect_true("coef.bFitMod" %in% methods("coef"))
})

test_that("bFitMod keeps covariance aligned when sorting doses", {
  dose <- c(0, 1, 3)
  resp <- c(.1, .5, .8)
  S <- diag(c(.02, .05, .09)) + .01
  perm <- c(3, 1, 2)
  set.seed(19)
  a <- bFitMod(dose, resp, "linear", S, type = "bootstrap", nSim = 8)
  set.seed(19)
  b <- bFitMod(dose[perm], resp[perm], "linear", S[perm, perm],
               type = "bootstrap", nSim = 8)
  expect_equal(attr(b, "data")$S, S)
  expect_equal(coef(a), coef(b))
})

test_that("legacy effect predictions subtract each draw's placebo response", {
  object <- list(samples = rbind(c(1, 2, 15), c(-1, 3, 15)))
  class(object) <- "bFitMod"
  attr(object, "model") <- "emax"
  attr(object, "placAdj") <- FALSE
  pred <- predict(object, doseSeq = c(15, 60), predType = "effect-curve",
                  summaryFct = NULL)
  expect_equal(unname(pred), rbind(c(1, 1.6), c(1.5, 2.4)))
  expect_equal(dim(predict(object, doseSeq = 15, summaryFct = NULL)), c(2L, 1L))
})
