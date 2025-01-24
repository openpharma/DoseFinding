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

#########################
## tests for ED and TD
#########################
data(biom)
modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                linInt = c(0, 0.5, 1, 1), doses = c(0, 0.05, 0.2, 0.6, 1))
## produce first stage fit (using dose as factor)
anMod <- lm(resp~factor(dose)-1, data=biom)
drFit <- coef(anMod)
S <- vcov(anMod)
dose <- sort(unique(biom$dose))
mod_dr <- fitMod(dose, drFit, S = S, type = "general", model = "emax",  bnds = c(0.01, 4))
prior <- list(norm = c(0, 10), norm = c(0,100), beta=c(0,1.5,0.45,1.7))
mod_bfit <- bFitMod(dose, drFit, S, model = "emax", 
                    start = c(0, 1, 0.1), nSim = 1000, prior = prior)
mod_maFit <- maFitMod(dose, drFit, S, model = c("emax", "sigEmax"), nSim = 10)

test_that("TD errors with type discrete if incorrect dose-range supplied", {
  expect_error(TD(modlist, Delta=0.3, TDtype = "discrete", doses=dose[-1]), "need placebo dose for TD calculation")
  expect_error(TD(mod_dr, Delta=0.3, TDtype = "discrete", doses=dose[-1]), "need placebo dose for TD calculation")
  expect_error(TD(mod_bfit, Delta=0.3, TDtype = "discrete", doses=dose[-1]), "need placebo dose for TD calculation")
  expect_error(TD(mod_maFit, Delta=0.3, TDtype = "discrete", doses=dose[-1]), "need placebo dose for TD calculation")
  
  expect_error(TD(modlist, Delta=0.3, TDtype = "discrete", doses=c(dose, 2)), "Doses provided may not exceed the observed dose range")
  expect_error(TD(mod_dr, Delta=0.3, TDtype = "discrete", doses=c(dose, 2)), "Doses provided may not exceed the observed dose range")
  expect_error(TD(mod_bfit, Delta=0.3, TDtype = "discrete", doses=c(dose, 2)), "Doses provided may not exceed the observed dose range")
  expect_error(TD(mod_maFit, Delta=0.3, TDtype = "discrete", doses=c(dose, 2)), "Doses provided may not exceed the observed dose range")
  
  
})

test_that("TD gives consistent results for discrete and continuous type", {
  td1a <- TD(modlist, Delta=0.3, TDtype = "discrete", doses=seq(0, max(dose), 0.002))
  td1b <- TD(modlist, Delta=0.3, TDtype = "discrete", doses=seq(0, max(dose) - 0.1, 0.002))
  td2 <- TD(modlist, Delta=0.3, TDtype = "continuous")
  
  expect_equal(td1a, td2, tolerance = 0.01)
  expect_equal(td1b, td2, tolerance = 0.01)
  
  td1a <- TD(mod_dr, Delta=0.3, TDtype = "discrete", doses=seq(0, max(dose), 0.002))
  td1b <- TD(mod_dr, Delta=0.3, TDtype = "discrete", doses=seq(0, max(dose) - 0.1, 0.002))
  td2 <- TD(mod_dr, Delta=0.3, TDtype = "continuous")
  expect_equal(td1a, td2, tolerance = 0.01)
  expect_equal(td1b, td2, tolerance = 0.01)
  
  td1a <- median(TD(mod_bfit, Delta=0.3, TDtype = "discrete", doses=seq(0, max(dose), 0.002)))
  td1b <- median(TD(mod_bfit, Delta=0.3, TDtype = "discrete", doses=seq(0, max(dose) - 0.1, 0.002)))
  td2 <- median(TD(mod_bfit, Delta=0.3, TDtype = "continuous"))
  expect_equal(td1a, td2, tolerance = 0.01)
  expect_equal(td1b, td2, tolerance = 0.01)
  
  td1a <- TD(mod_maFit, Delta=0.3, TDtype = "discrete", doses=seq(0, max(dose), 0.002))
  td1b <- TD(mod_maFit, Delta=0.3, TDtype = "discrete", doses=seq(0, max(dose) - 0.1, 0.002))
  td2 <- TD(mod_maFit, Delta=0.3, TDtype = "continuous")
  expect_equal(td1a, td2, tolerance = 0.01)
  expect_equal(td1b, td2, tolerance = 0.01)
  
})


test_that("ED errors with type discrete if incorrect dose-range supplied", {
  expect_error(ED(modlist, p=0.9, EDtype = "discrete", doses=dose[-1]), "need placebo dose for ED calculation")
  expect_error(ED(mod_dr, p=0.9, EDtype = "discrete", doses=dose[-1]), "need placebo dose for ED calculation")
  expect_error(ED(mod_bfit, p=0.9, EDtype = "discrete", doses=dose[-1]), "need placebo dose for ED calculation")
  expect_error(ED(mod_maFit, p=0.9, EDtype = "discrete", doses=dose[-1]), "need placebo dose for ED calculation")
  
  expect_error(ED(modlist, p=0.9, EDtype = "discrete", doses=c(dose, 2)), "Doses provided may not exceed the observed dose range")
  expect_error(ED(mod_dr, p=0.9, EDtype = "discrete", doses=c(dose, 2)), "Doses provided may not exceed the observed dose range")
  expect_error(ED(mod_bfit, p=0.9, EDtype = "discrete", doses=c(dose, 2)), "Doses provided may not exceed the observed dose range")
  expect_error(ED(mod_maFit, p=0.9, EDtype = "discrete", doses=c(dose, 2)), "Doses provided may not exceed the observed dose range")
  
  
})

test_that("ED gives consistent results for discrete and continuous type", {
  ed1a <- ED(modlist, p=0.9, EDtype = "discrete", doses=seq(0, max(dose), 0.002))
  ed1b <- ED(modlist, p=0.9, EDtype = "discrete", doses=seq(0, max(dose) - 0.05, 0.002))
  ed2 <- ED(modlist, p=0.9, EDtype = "continuous")
  expect_equal(ed1a, ed1b)
  expect_equal(ed1a, ed2, tolerance = 0.01)
  expect_equal(ed1b, ed2, tolerance = 0.01)
  
  ed1a <- ED(mod_dr, p=0.9, EDtype = "discrete", doses=seq(0, max(dose), 0.002))
  ed1b <- ED(mod_dr, p=0.9, EDtype = "discrete", doses=seq(0, max(dose) - 0.05, 0.002))
  ed2 <- ED(mod_dr, p=0.9, EDtype = "continuous")
  expect_equal(ed1a, ed1b)
  expect_equal(ed1a, ed2, tolerance = 0.01)
  expect_equal(ed1b, ed2, tolerance = 0.01)
  
  ed1a <- median(ED(mod_bfit, p=0.9, EDtype = "discrete", doses=seq(0, max(dose), 0.002)))
  ed1b <- median(ED(mod_bfit, p=0.9, EDtype = "discrete", doses=seq(0, max(dose) - 0.05, 0.002)))
  ed2 <- median(ED(mod_bfit, p=0.9, EDtype = "continuous"))
  expect_equal(ed1a, ed1b)
  expect_equal(ed1a, ed2, tolerance = 0.01)
  expect_equal(ed1b, ed2, tolerance = 0.01)
  
  ed1a <- ED(mod_maFit, p=0.9, EDtype = "discrete", doses=seq(0, max(dose), 0.002), direction = "increasing")
  ed1b <- ED(mod_maFit, p=0.9, EDtype = "discrete", doses=seq(0, max(dose) - 0.05, 0.002), direction = "increasing")
  ed2 <- ED(mod_maFit, p=0.9, EDtype = "continuous", direction = "increasing")
  expect_equal(ed1a, ed1b)
  expect_equal(ed1a, ed2, tolerance = 0.01)
  expect_equal(ed1b, ed2, tolerance = 0.01)
  
})