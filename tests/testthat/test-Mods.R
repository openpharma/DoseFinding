test_that("Mods function requires dose levels", {
  expect_error(Mods(linear = NULL), "Need to specify dose levels")
})

test_that("Mods function ensures dose levels include placebo and are non-negative", {
  expect_error(Mods(linear = NULL, doses = c(0.05, 0.2)), "Need to include placebo dose")
  expect_error(Mods(linear = NULL, doses = c(-0.05, 0, 0.2)), "Only dose-levels >= 0 allowed")
})

test_that("Mods function checks addArgs parameters for validity", {
  expect_error(Mods(linear = NULL, doses = c(0, 0.05, 0.2), addArgs = list(scal = 0.1, off = 0.01)), 
               "\"scal\" parameter needs to be ")
  expect_error(Mods(linear = NULL, doses = c(0, 0.05, 0.2), addArgs = list(scal = 1.2, off = -0.1)), 
               "\"off\" parameter needs to be positive")
})

test_that("Mods function generates an object of class Mods", {
  models <- Mods(linear = NULL, emax = 0.05, doses = c(0, 0.05, 0.2, 0.6, 1), addArgs = list(scal = 1.2, off = 0.1))
  expect_s3_class(models, "Mods")
  expect_true(!is.null(attr(models, "placEff")))
  expect_true(!is.null(attr(models, "maxEff")))
  expect_true(!is.null(attr(models, "direction")))
  expect_true(!is.null(attr(models, "doses")))
  expect_true(!is.null(attr(models, "scal")))
  expect_true(!is.null(attr(models, "off")))
})

test_that("Mods function calculates responses correctly", {
  doses <- c(0, 10, 25, 50, 100, 150)
  fmodels <- Mods(linear = NULL, emax = 25, 
                  logistic = c(50, 10.88111), exponential = 85, 
                  betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)), 
                  linInt = rbind(c(0, 1, 1, 1, 1), c(0, 0, 1, 1, 0.8)), 
                  doses = doses, placEff = 0.5, maxEff = -0.4, 
                  addArgs = list(scal = 200))
  responses <- getResp(fmodels, doses)
  expect_equal(nrow(responses), length(doses))
})

test_that("Mods function can specify all model parameters (fullMod = TRUE)", {
  fmods <- Mods(emax = c(0, 1, 0.1), linear = cbind(c(-0.4, 0), c(0.2, 0.1)), 
                sigEmax = c(0, 1.1, 0.5, 3), 
                doses = 0:4, fullMod = TRUE)
  responses <- getResp(fmods, doses = seq(0, 4, length = 11))
  expect_equal(nrow(responses), 11)
  expect_equal(ncol(responses), length(attr(fmods, "maxEff")))
})


## test plotting functions
test_that("plotMods function basic functionality", {
  models <- Mods(linear = NULL, emax = 0.05, doses = c(0, 0.05, 0.2, 0.6, 1), 
                 addArgs = list(scal = 1.2, off = 0.1))
  p <- plotMods(models)
  
  expect_s3_class(p, "ggplot")
  expect_true("GeomLine" %in% sapply(p$layers, function(layer) class(layer$geom)[1]))
  expect_true("GeomPoint" %in% sapply(p$layers, function(layer) class(layer$geom)[1]))
  
  p_superpose <- plotMods(models, superpose = TRUE)
  expect_s3_class(p_superpose, "ggplot")
  expect_true("GeomLine" %in% sapply(p_superpose$layers, function(layer) class(layer$geom)[1]))
})

test_that("plot.Mods function basic functionality", {
  models <- Mods(linear = NULL, emax = 0.05, doses = c(0, 0.05, 0.2, 0.6, 1), 
                 addArgs = list(scal = 1.2, off = 0.1))
  
  p <- plot(models)
  
  expect_s3_class(p, "trellis")
})

test_that("plotMods handles customizations correctly", {
  models <- Mods(linear = NULL, emax = 0.05, doses = c(0, 0.05, 0.2, 0.6, 1), 
                 addArgs = list(scal = 1.2, off = 0.1))
  
  p_custom <- plotMods(models, xlab = "Custom X Label", ylab = "Custom Y Label")
  
  expect_s3_class(p_custom, "ggplot")
  expect_equal(p_custom$labels$x, "Custom X Label")
  expect_equal(p_custom$labels$y, "Custom Y Label")
})

test_that("plot.Mods handles customizations correctly", {
  models <- Mods(linear = NULL, emax = 0.05, doses = c(0, 0.05, 0.2, 0.6, 1), 
                 addArgs = list(scal = 1.2, off = 0.1))
  
  p_custom <- plot(models, lwd = 3, pch = 3, cex = 1.2, col = "red")
  
  expect_s3_class(p_custom, "trellis")
})