context("dose response model functions")

ud <- function(x) unname(drop(x))

test_that("betaMod does not produce NaN for large delta1, delta2", {
  expect_equal(betaMod(100, 1, 2, 10, 10, 200), 3)
  expect_equal(betaMod(100, 1, 2, 150, 150, 200), 3)
  expect_equal(betaMod(100, 1, 2, 100, 50, 200), 1.000409)
  expect_equal(betaMod(0, 1, 2, 50, 50, 200), 1)
  expect_equal(betaMod(0, 1, 2, 75, 75, 200), 1)
  expect_equal(ud(betaModGrad(100, 2, 50, 50, 200)), c(1, 1, 0, 0))
  expect_equal(ud(betaModGrad(100, 2, 150, 150, 200)), c(1, 1, 0, 0))
  expect_equal(ud(betaModGrad(0, 2, 50, 50, 200)), c(1, 0, 0, 0))
  expect_equal(ud(betaModGrad(0, 2, 100, 100, 200)), c(1, 0, 0, 0))
})

test_that("sigEmax does not produce NaN for large dose and large h", {
  expect_equal(sigEmax(100, 1, 1, 50, 2), 1.8)
  expect_equal(sigEmax(100, 1, 1, 50, 150), 2)
  expect_equal(sigEmax(150, 1, 1, 50, 150), 2)
  expect_equal(sigEmax(0, 1, 1, 50, 10), 1)
  expect_equal(sigEmax(0, 1, 1, 50, 400), 1)
  expect_equal(sigEmax(c(50, 150), 1, 1, 50, 0), c(1.5, 1.5))
  expect_equal(ud(sigEmaxGrad(100, 1, 50, 10)), c(1, 0.999024390243902, -0.000194931588340274, 0.000675581404300663))
  expect_equal(ud(sigEmaxGrad(100, 1, 50, 150)), c(1, 1, 0, 0))
  expect_equal(ud(sigEmaxGrad(150, 1, 50, 150)), c(1, 1, 0, 0))
  expect_equal(ud(sigEmaxGrad(0, 1, 50, 0)), c(1, 0.5, 0, 0))
  expect_equal(ud(sigEmaxGrad(0, 1, 50, 150)), c(1, 0, 0, 0))
  # this is the only NaN we can't get rid off, as the function
  #   (a,b,x) ↦ a^x/(a^x+b^x)
  # has a non-removable discontinuity at (0, 0, x) for all x > 0
  # fortunately an ed50=0 does not make much sense from a modeling perspective
  expect_equal(sigEmax(0, 1, 1, 0, 5), NaN)
})
