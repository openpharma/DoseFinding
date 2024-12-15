context("maFitMod")

data(biom)
anMod <- lm(resp~factor(dose)-1, data=biom)
drFit <- coef(anMod)
S <- vcov(anMod)
doses <- c(0,0.05,0.2,0.6,1)


test_that("Test input parameters of maFitMod", {
  ## models
  expect_error(
    maFitMod(doses, drFit, S),
    "Need to specify the models that should be fitted")
  expect_error(
    maFitMod(doses, drFit, S, models = c("linear", "emax", "XYZ")),
    "Invalid dose-response model specified: XYZ")
  expect_error(
    maFitMod(doses, drFit, S, models = c("ABC", "emax", "XYZ")),
    "Invalid dose-response model specified: ABC, XYZ")
  ## bnds
  bnds <- c(10,20)
  expect_error(
    maFitMod(doses, drFit, S, models = c("emax"), nSim = 10, bnds = bnds),
    "bnds needs to be a list")

  bnds <- defBnds(1)
  bnds$emax <- NULL
  expect_message(
    maFitMod(doses, drFit, S, models = c("emax"), nSim = 10, bnds = bnds))
  
  bnds$emax <- c(10,20)
  tmp <- maFitMod(doses, drFit, S, models = c("emax"), nSim = 10,
                  bnds = bnds)
  out <- sapply(tmp$fits, function(x) coef(x)[3]) # all should be in [10,20]
  eps <- 0.0001
  for(i in 1:10){
    expect_lt(out[i], 20+eps)
    expect_gt(out[i], 10-eps)
  }
  
  bnds$betaMod <- rbind(c(1,2), c(1,2))
  tmp <- maFitMod(doses, drFit, S, models = c("betaMod"), nSim = 10,
                  bnds = bnds)
  out <- sapply(tmp$fits, function(x) coef(x)[3:4]) # all should be in [1,2]
  eps <- 0.0001
  for(i in 1:10){
    expect_lt(out[1,i], 2+eps)
    expect_gt(out[1,i], 1-eps)
    expect_lt(out[2,i], 2+eps)
    expect_gt(out[2,i], 1-eps)
  }

  ## addArgs
  tmp <- maFitMod(doses, drFit, S, models = c("betaMod"), nSim = 10)
  out <- sapply(tmp$fits, function(x) attr(x, "scal")) # should be 1.2
  for(i in 1:10){
    expect_equal(out[i], 1.2)
  }
  tmp <- maFitMod(doses, drFit, S, models = c("betaMod"), nSim = 10,
                  addArgs = list(scal = 200)) 
  out <- sapply(tmp$fits, function(x) attr(x, "scal")) # should be 200
  for(i in 1:10){
    expect_equal(out[i], 200)
  }
  tmp <- maFitMod(doses, drFit, S, models = c("linlog"), nSim = 10,
                  addArgs = list(off = 123)) 
  out <- sapply(tmp$fits, function(x) attr(x, "off")) # should be 123
  for(i in 1:10){
    expect_equal(out[i], 123)
  }
})

test_that("test model fitting", {
  set.seed(295)
  bnds <- defBnds(max(doses))
  expect_silent(fits1 <- maFitMod(doses, drFit, S, models = c("linear"), nSim = 10))
  expect_silent(fits2 <- maFitMod(doses, drFit, S, models = c("linInt"), nSim = 10))
  expect_silent(fits3 <- maFitMod(doses, drFit, S, models = c("emax", "sigEmax"), nSim = 10, bnds = bnds))
  expect_silent(fits4 <- maFitMod(doses, drFit, S, models = c("linear", "emax", "betaMod"), nSim = 10, bnds = bnds))
  expect_silent(fits5 <- maFitMod(doses, drFit, S, models = c("linear", "emax", "betaMod"), nSim = 10, bnds = bnds))
  builtin <- c("linlog", "linear", "quadratic", "linInt", "emax",
               "exponential", "logistic", "betaMod", "sigEmax")
  expect_silent(fits6 <- maFitMod(doses, drFit, S, models = builtin, nSim = 10, bnds = bnds))
  expect_true(class(fits6) == "maFit")
  expect_equal(length(fits6$fits), 10)
  expect_equal(length(fits6$selModels), 10)
  expect_named(fits6, c("fits", "selModels", "args"))
  
  ## test print method (HOW TO TEST PRINT, WITHOUT PRINTING TO CONSOLE)
  ## expect_no_condition(print(fits1))
  ## expect_no_condition(print(fits2, digits=10))
  ## expect_no_condition(print(fits3))
  ## expect_no_condition(print(fits4)
  ## expect_no_condition(print(fits5))
  ## expect_no_condition(print(fits6))
  
  ## test prediction
  expect_error(
    predict(fits1),
    "Need to provide doseSeq argument")
  dsq <- seq(0,1,length=101)
  expect_silent(p1 <- predict(fits1, doseSeq = dsq))
  expect_silent(p2 <- predict(fits2, doseSeq = dsq))
  expect_silent(p3 <- predict(fits3, doseSeq = dsq))
  expect_silent(p4 <- predict(fits4, doseSeq = dsq))
  expect_silent(p5 <- predict(fits5, doseSeq = dsq))
  expect_silent(p6 <- predict(fits5, doseSeq = dsq))
  expect_equal(dim(p6), c(10,101))
  
  expect_silent(plot(fits1))
  expect_silent(plot(fits2, xlab = "ABC"))
  expect_silent(plot(fits3, ylab = "XYZ"))
  expect_silent(plot(fits4, title = "123"))
  expect_silent(plot(fits5))
  expect_silent(plot(fits6))
  
  expect_silent(plot(fits1, plotData = "none"))
  expect_silent(plot(fits2, plotData = "none"))
  expect_silent(plot(fits3, plotData = "none"))
  expect_silent(plot(fits4, plotData = "none"))
  expect_silent(plot(fits5, plotData = "none"))
  expect_silent(plot(fits6, plotData = "none"))
  
  expect_silent(plot(fits1, plotData = "meansCI"))
  expect_silent(plot(fits2, plotData = "meansCI"))
  expect_silent(plot(fits3, plotData = "meansCI"))
  expect_silent(plot(fits4, plotData = "meansCI"))
  expect_silent(plot(fits5, plotData = "meansCI"))
  expect_silent(plot(fits6, plotData = "meansCI"))
  expect_error(plot(fits6, title = 23),
               "title needs to be a character")
  expect_error(plot(fits6, trafo = 23),
               "trafo needs to be a function")
  
  expect_error(
    ED.maFit(fits5, p = 0.9)) # should fail (direction not specified)
  expect_silent(ED.maFit(fits5, p = 0.9, direction = "increasing"))
  expect_silent(ED.maFit(fits5, p = 0.5, direction = "increasing"))
  
  ## check decreasing direction
  drFit2 <- 1-drFit
  fits6 <- maFitMod(doses, drFit2, S, models = builtin, nSim = 10, bnds = bnds)
  expect_true(
    is.na(ED.maFit(fits6, p = 0.9, direction = "increasing")))
  expect_silent(ED.maFit(fits6, p = 0.9, direction = "decreasing"))
  expect_silent(ED.maFit(fits6, p = 0.5, direction = "decreasing"))
})

## compare model fitting to bFitMod
## Helper function to create example data
make_example_data <- function() {
  dose <- c(0, 5, 10, 20, 40)
  resp <- sort(rnorm(5,0,0.1))
  S <- diag(5)*rexp(5,2) + crossprod(matrix(rexp(25,10),5,5))
  list(ds = dose, means = resp, vc = S)
}
test_that("Compare fits of maFitMod to bFitMod", {
  bnds <- defBnds(40)
  set.seed(726)
  dat <- make_example_data()
  builtin <- c("linlog", "linear", "quadratic", "linInt", "emax",
               "exponential", "logistic", "betaMod", "sigEmax")
  for(i in 1:length(builtin)){
    set.seed(376)
    res1 <- maFitMod(dat$ds, dat$means, S=dat$vc, models = builtin[i], nSim = 10,
                     bnds = bnds)
    set.seed(376)
    res2 <- bFitMod(dat$ds, dat$means, S=dat$vc, model = builtin[i], type = "bootstrap",
                    nSim = 10, bnds = bnds[[builtin[i]]])
    
    ds <- seq(0,40,by=1)
    pred0 <- predict(res1, doseSeq = ds)
    pred1 <- apply(pred0, 2, function(x){
      quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))
    })
    pred2 <- predict(res2, doseSeq = ds)
    rownames(pred1) <- NULL
    expect_equal(pred1, pred2)
  }
})
  
