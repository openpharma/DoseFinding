get_pow_mct_interim_example <- function() {
  doses <- c(0, 0.5, 1, 2, 4, 8)
  mods <- Mods(
    emax = c(0.5, 1, 2, 4),
    sigEmax = rbind(c(0.5, 3), c(1, 3), c(2, 3), c(4, 3)),
    quadratic = -0.1,
    doses = doses
  )
  w <- c(1, 0.5, 0.5, 0.5, 1, 1)
  contMat <- optContr(models = mods, w = w)$contMat
  sigma <- 0.3
  n_final <- round(531 * w / sum(w))
  n <- floor(n_final / 2)
  S0t <- diag(sigma^2 / n)
  S_end <- diag(sigma^2 / n_final)
  ## assumed interim estimate
  mu0t <- 0.05 * doses / (doses + 1) + rnorm(6, 0, 0.382 / sqrt(n))
  ## assumed mu (needed for conditional power)
  mu_assumed <- 0.135 * doses / (doses + 1)
  list(
    contMat = contMat,
    S0t = S0t,
    S_end = S_end,
    mu0t = mu0t,
    mu_assumed = mu_assumed
  )
}

test_that("powMCTInterim works as expected with predictive power", {
  example_data <- get_pow_mct_interim_example()
  result <- powMCTInterim(
    contMat = example_data$contMat,
    S0t = example_data$S0t,
    S_end = example_data$S_end,
    mu0t = example_data$mu0t,
    type = "predictive"
  )
  expect_snapshot_value(result, style = "deparse")
})

test_that("powMCTInterim works as expected with conditional power", {
  example_data <- get_pow_mct_interim_example()
  result <- powMCTInterim(
    contMat = example_data$contMat,
    S0t = example_data$S0t,
    S_end = example_data$S_end,
    mu0t = example_data$mu0t,
    mu = example_data$mu_assumed,
    type = "conditional"
  )
  expect_snapshot_value(result, style = "deparse")
})
