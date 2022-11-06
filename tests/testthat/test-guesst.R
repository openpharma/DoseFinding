context("guesstimates")

test_that("emax", {
  emx1 <- guesst(d=0.3, p=0.8, model="emax")
  expect_equal(unname(emax(0.3,0,1,emx1)),
               0.8,
               tolerance = 0.001)
  
})

test_that("emax local", {
  emx2 <- guesst(d=0.3, p=0.8, model="emax", local = TRUE, Maxd = 1)
  expect_equal(unname(emax(0.3,0,1,emx2)/emax(1,0,1,emx2)),
               0.8,
               tolerance = 0.001)
})

test_that("betaMod", {
  bta <- guesst(d=0.4, p=0.8, model="betaMod", dMax=0.8, scal=1.2, Maxd=1)
  expect_equal(betaMod(c(0.4,0.8), 0, 1, bta[1], bta[2], scal=1.2),
               c(0.8, 1.0),
               tolerance = 0.001)
})

test_that("exponential", {
  expo <- guesst(d = 0.8, p = 0.5, "exponential", Maxd=1)
  expect_equal(unname(exponential(0.8,0,1,expo)/exponential(1,0,1,expo)),
               0.5,
               tolerance = 0.001)
})

test_that("quadratic", {
  quad <- guesst(d = 0.7, p = 1, "quadratic")
  mm <- Mods(quadratic=quad, doses=c(0,0.7,1))
  expect_equal(getResp(mm)[2],
               1,
               tolerance = 0.001)
})

test_that("logistic", {
  lgc1 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "logistic")
  expect_equal(logistic(c(0.2,0.6), 0, 1, lgc1[1], lgc1[2]),
               c(0.2, 0.95),
               tolerance = 0.001)
})

test_that("logistic local", {
  lgc2 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "logistic", 
                 local = TRUE, Maxd = 1)
  r0 <- logistic(0, 0, 1, lgc2[1], lgc2[2])
  r1 <- logistic(1, 0, 1, lgc2[1], lgc2[2])
  expect_equal((logistic(c(0.2,0.6), 0, 1, lgc2[1], lgc2[2])-r0)/(r1-r0),
               c(0.2, 0.95),
               tolerance = 0.001)
})

test_that("sigEmax", {
  sgE1 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "sigEmax")
  expect_equal(sigEmax(c(0.2,0.6), 0, 1, sgE1[1], sgE1[2]),
               c(0.2, 0.95),
               tolerance = 0.001)
})

test_that("sigEmax local", {
  sgE2 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "sigEmax",
                 local = TRUE, Maxd = 1)
  r1 <- sigEmax(1, 0, 1, sgE2[1], sgE2[2])
  expect_equal(sigEmax(c(0.2,0.6), 0, 1, sgE2[1], sgE2[2])/r1,
               c(0.2,0.95),
               tolerance = 0.001)
})

