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
  S_0t <- diag(sigma^2 / n)
  S_01 <- diag(sigma^2 / n_final)
  ## assumed interim estimate
  mu_0t <- 0.05 * doses / (doses + 1) + rnorm(6, 0, 0.382 / sqrt(n))
  ## assumed mu (needed for conditional power)
  mu_assumed <- 0.135 * doses / (doses + 1)
  list(
    contMat = contMat,
    S_0t = S_0t,
    S_01 = S_01,
    mu_0t = mu_0t,
    mu_assumed = mu_assumed
  )
}

test_that("powMCTInterim works as expected with predictive power", {
  set.seed(123)
  example_data <- get_pow_mct_interim_example()
  result <- powMCTInterim(
    contMat = example_data$contMat,
    S_0t = example_data$S_0t,
    S_01 = example_data$S_01,
    mu_0t = example_data$mu_0t,
    type = "predictive"
  )
  expect_equal(as.numeric(result), 0.9473, tolerance = 1e-4)
})

test_that("powMCTInterim works as expected with conditional power", {
  set.seed(321)
  example_data <- get_pow_mct_interim_example()
  result <- powMCTInterim(
    contMat = example_data$contMat,
    S_0t = example_data$S_0t,
    S_01 = example_data$S_01,
    mu_0t = example_data$mu_0t,
    mu_assumed = example_data$mu_assumed,
    type = "conditional"
  )
  expect_equal(as.numeric(result), 0.2739, tolerance = 5e-4)
})

test_that("powMCTInterim works as expected when using NULL explicitly for mu_assumed", {
  set.seed(321)
  example_data <- get_pow_mct_interim_example()
  expect_message(
    result <- powMCTInterim(
      contMat = example_data$contMat,
      S_0t = example_data$S_0t,
      S_01 = example_data$S_01,
      mu_0t = example_data$mu_0t,
      mu_assumed = NULL,
      type = "conditional"
    ),
    "mu_assumed not supplied, setting mu_assumed = mu_0t"
  )
  expect_equal(as.numeric(result), 0.00933, tolerance = 1e-4)
})

test_that("powMCTInterim can specify control via list", {
  set.seed(123)
  example_data <- get_pow_mct_interim_example()
  result <- powMCTInterim(
    contMat = example_data$contMat,
    S_0t = example_data$S_0t,
    S_01 = example_data$S_01,
    mu_0t = example_data$mu_0t,
    type = "predictive",
    control = list(maxpts = 1e5)
  )
  expect_equal(as.numeric(result), 0.9473, tolerance = 1e-4)
})

# Simulation based calculation, used for the integration tests below.
interim_power <- function(
  nSim,
  contMat,
  S0t,
  S_end,
  mu0t,
  type = c("predictive", "conditional"),
  mu_assumed,
  alpha
) {
  S0t_inv <- solve(S0t)
  St1_inv <- solve(S_end) - S0t_inv
  St1 <- solve(St1_inv)

  ## simulate incremental information for second stage data
  if (type == "predictive") {
    mu_stage2 <- rmvnorm(nSim, mu0t, S0t + St1)
  }
  if (type == "conditional") {
    mu_stage2 <- rmvnorm(nSim, mu_assumed, St1)
  }

  ## pre-calculate critical value
  covMat <- t(contMat) %*% S_end %*% contMat
  corMat <- cov2cor(covMat)
  critV <- critVal(corMat, alpha = alpha, df = Inf)
  den <- sqrt(diag(covMat)) ## numerator of t-statistics
  ## now determine for each sampled second stage estimate, whether final test
  ## would be significant
  maxT <- numeric(nSim)
  for (i in 1:nSim) {
    ## combine estimates from part 1 and 2
    mu_end <- as.numeric(
      S_end %*% (S0t_inv %*% mu0t + St1_inv %*% mu_stage2[i, ])
    )
    ## calculate maximum statistic
    ct <- as.vector(mu_end %*% contMat)
    tStat <- ct / den
    maxT[i] <- max(tStat)
  }
  ## final power value
  mean(maxT > critV)
}

test_that("powMCTInterim gives same conditional power result as with simulation based approach", {
  skip_on_cran()
  skip_on_ci()

  set.seed(245)
  example_data <- get_pow_mct_interim_example()

  result <- powMCTInterim(
    contMat = example_data$contMat,
    S_0t = example_data$S_0t,
    S_01 = example_data$S_01,
    mu_0t = example_data$mu_0t,
    mu_assumed = example_data$mu_assumed,
    type = "conditional"
  )

  expected <- interim_power(
    nSim = 1e6,
    contMat = example_data$contMat,
    S0t = example_data$S_0t,
    S_end = example_data$S_01,
    mu0t = example_data$mu_0t,
    type = "conditional",
    mu_assumed = example_data$mu_assumed,
    alpha = 0.025
  )

  expect_equal(result, expected, tolerance = 1e-2, ignore_attr = TRUE)
})

test_that("powMCTInterim gives same predictive power result as with simulation based approach", {
  skip_on_cran()
  skip_on_ci()

  set.seed(245)
  example_data <- get_pow_mct_interim_example()

  result <- powMCTInterim(
    contMat = example_data$contMat,
    S_0t = example_data$S_0t,
    S_01 = example_data$S_01,
    mu_0t = example_data$mu_0t,
    type = "predictive"
  )

  expected <- interim_power(
    nSim = 1e6,
    contMat = example_data$contMat,
    S0t = example_data$S_0t,
    S_end = example_data$S_01,
    mu0t = example_data$mu_0t,
    type = "predictive",
    alpha = 0.025
  )

  expect_equal(result, expected, tolerance = 1e-2, ignore_attr = TRUE)
})
