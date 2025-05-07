context("Optimal Contrasts")

require_extra_packages <- function() {
  if (!(require("quadprog") && require("Rsolnp"))) {
    skip("packages quadprog and Rsolnp not available")
  }
}

# calculation of optimal contrast by enumerating all active sets
allActiveSets <- function(S, mu, mult){
  k <- length(mu)
  CC <- cbind(-1, diag(k - 1))
  SPa <- CC %*% S %*% t(CC)
  muPa <- as.numeric(CC %*% mu)
  # generate all possible active sets
  mat <- matrix(nrow = 2^(k-1), ncol = (k-1))
  for(i in 1:(k-1))
    mat[,i] <- rep(rep(c(FALSE,TRUE), each=2^(i-1)), 2^((k-1)-i))
  val <- numeric(2^(k-1))
  feasible <- logical(2^(k-1))
  cont <- matrix(nrow = 2^(k-1), ncol = (k-1))
  for(i in 1:(2^(k-1))){
    nonzero <- mat[i,]
    if(sum(nonzero) > 0){
      cont[i,!nonzero] <- 0
      cont[i,nonzero] <- solve(SPa[nonzero, nonzero]) %*% muPa[nonzero]
      feasible[i] <- all(mult*cont[i,] >= 0)
      contrast <- c(-sum(cont[i,]), cont[i,])
      val[i] <- as.numeric(t(contrast)%*%mu/sqrt(t(contrast)%*%S%*%contrast))
    }
  }
  if(!any(feasible))
    return(rep(NA, k))
  mm <- max(val[which(feasible)])
  c(-sum(cont[val == mm,]), cont[val == mm,])
}

# helper functions
getStand <- function(x) x/sqrt(sum(x^2))
getNCP <- function(cont, mu, S) {
  as.numeric(t(cont)%*%mu/sqrt(t(cont)%*%S%*%cont))
}

one_sim <- function() {
  cont <- vector("list", 5)
  # simulate mean and covariance matrix
  kk <- round(runif(1, 4, 10))
  A <- matrix(runif(kk^2, -1, 1), kk, kk)
  S <- crossprod(A)+diag(kk)
  S_inv <- solve(S)
  mult <- sign(rnorm(1))
  mu <- mult*sort(rnorm(kk, 1:kk, 1))
  # unconstrained solution
  ones <- rep(1, kk)
  unConst <- S_inv%*%(mu - c(t(mu)%*%S_inv%*%ones/(t(ones)%*%S_inv%*%ones)))
  cont[[1]] <- getStand(unConst)
  # function from DoseFinding package
  cont[[2]] <- DoseFinding:::constOptC(mu, S_inv, placAdj=FALSE,
                                   ifelse(mult == 1, "increasing", "decreasing"))
  # alternative solution using quadratic programming
  A <- t(rbind(rep(1, kk), mu,
               mult * diag(kk) * c(-1, rep(1, kk - 1))))
  bvec <- c(0, 1, rep(0, kk))
  rr <- solve.QP(S, rep(0, kk), A, bvec, meq = 2)
  cont[[3]] <- getStand(rr$solution)
  # using solnp
  mgetNCP <- function(x, ...){
    cont <- c(-sum(x), x)
    -getNCP(cont, ...)
  }
  res <- solnp(rep(1, kk-1), mgetNCP, mu=mu, S=S,
               LB=rep(0, kk-1), UB=rep(20, kk-1),
               control = list(trace = 0))
  cont[[4]] <- getStand(c(-sum(res$pars), res$pars))
  # using enumeration
  cont[[5]] <- allActiveSets(S=S, mu=mu, mult=mult)
  return(sapply(cont, getNCP, mu = mu, S = S))
}

test_that("calculation of contrasts works", {
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
  require_extra_packages()
  ncps <- replicate(1000, one_sim())
  ## calculate best result among alternative methods (solnp sometimes fails)
  best_ncp <- apply(ncps[c(3,4,5),], 2, max)
  ## compare to DoseFinding::constOptC
  expect_equal(ncps[2,], best_ncp)
})

test_that("constant shapes are handled correctly", {
  data(biom)
  # define shapes for which to calculate optimal contrasts
  modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                  linInt = rbind(c(0, 0, 0, 1), c(0, 1, 1, 1)),
                  doses = c(0, 0.05, 0.2, 0.6, 1), placEff = 1)
  cont_mat <- function(doses, placAdj, type) {
    optContr(modlist, w=1, doses=doses, placAdj=placAdj, type = type)$contMat
  }
  ## code should notice that linInt shapes are constant over specified dose rng (no contrast can be calculated)
  expect_message(cont_mat(0.05, TRUE, "u"), "The linInt1, linInt2 models have a constant shape, cannot
calculate optimal contrasts for these shapes.")
  expect_message(cont_mat(0.05, TRUE, "c"), "The linInt1, linInt2 models have a constant shape, cannot
calculate optimal contrasts for these shapes.")
  expect_message(cont_mat(c(0.05, 0.5), TRUE, "u"), "The linInt1 model has a constant shape, cannot
calculate optimal contrasts for this shape.")
  expect_message(cont_mat(c(0.05, 0.5), TRUE, "c"), "The linInt1 model has a constant shape, cannot
calculate optimal contrasts for this shape.")
  expect_message(cont_mat(c(0, 0.05), FALSE, "u"), "The linInt1, linInt2 models have a constant shape, cannot
calculate optimal contrasts for these shapes.")
  expect_message(cont_mat(c(0, 0.05), FALSE, "c"), "The linInt1, linInt2 models have a constant shape, cannot
calculate optimal contrasts for these shapes.")
  expect_message(cont_mat(c(0, 0.05, 0.5), FALSE, "u"), "The linInt1 model has a constant shape, cannot
calculate optimal contrasts for this shape.")
  expect_message(cont_mat(c(0, 0.05, 0.5), FALSE, "c"), "The linInt1 model has a constant shape, cannot
calculate optimal contrasts for this shape.")
  ## in case of all constant shapes stop with error
  modlist2 <- Mods(linInt = rbind(c(0, 1, 1, 1), c(0, 0, 0, 1)),
                   doses = c(0, 0.05, 0.2, 0.6, 1), placEff = 1)
  expect_error(optContr(modlist2, w=1, doses=c(0.05), placAdj=TRUE, type = "u"),
               "All models correspond to a constant shape, no optimal contrasts calculated.")
  expect_error(optContr(modlist2, w=1, doses=c(0.05), placAdj=TRUE, type = "c"),
               "All models correspond to a constant shape, no optimal contrasts calculated.")
  expect_error(optContr(modlist2, w=1, doses=c(0, 0.05), placAdj=FALSE, type = "u"),
               "All models correspond to a constant shape, no optimal contrasts calculated.")
  expect_error(optContr(modlist2, w=1, doses=c(0, 0.05), placAdj=FALSE, type = "c"),
               "All models correspond to a constant shape, no optimal contrasts calculated.")
  ## mixed cases where some linInt models are non-constant
  expect_message(optContr(modlist2, w=1, doses=c(0.05, 0.5), placAdj=TRUE, type = "u"), "The linInt2 model has a constant shape, cannot
calculate optimal contrasts for this shape.")
  expect_message(optContr(modlist2, w=1, doses=c(0.05, 0.5), placAdj=TRUE, type = "c"), "The linInt2 model has a constant shape, cannot
calculate optimal contrasts for this shape.")
  expect_message(optContr(modlist2, w=1, doses=c(0, 0.05, 0.5), placAdj=FALSE, type = "u"), "The linInt2 model has a constant shape, cannot
calculate optimal contrasts for this shape.") 
  expect_message(optContr(modlist2, w=1, doses=c(0, 0.05, 0.5), placAdj=FALSE, type = "c"), "The linInt2 model has a constant shape, cannot
calculate optimal contrasts for this shape.")
})

test_that("optContr errors when invalid inputs are provided", {
  expect_error(optContr(models = list(), doses = c(0, 10), w = c(1, 1)),
               "models needs to be of class Mods")
  models <- Mods(linear = NULL, emax = 25, direction = c("increasing", "decreasing"), doses = c(0, 10))
  models <- Mods(linear = NULL, doses = c(0, 10))
  expect_error(optContr(models, doses = c(0, 10)),
               "Need to specify exactly one of \"w\" or \"S\"")
  expect_error(optContr(models, doses = c(0, 10), w = c(1, 1), S = diag(2)),
               "Need to specify exactly one of \"w\" or \"S\"")
  expect_error(optContr(models, doses = c(0, 10), w = c(1, 1), placAdj = TRUE),
               "If placAdj == TRUE there should be no placebo group in \"doses\"")
  expect_error(optContr(models, doses = c(0, 10), w = c(1, 1, 1)),
               "w needs to be of length 1 or of the same length as doses")
  expect_error(optContr(models, doses = c(0, 10), S = c(1, 1)),
               "S needs to be a matrix")
})

models <- Mods(linear = NULL, doses = c(0, 10))

test_that("print.optContr prints contrast matrix", {
  contMat <- optContr(models, doses = c(0, 10), w = c(1, 1))
  expect_output(print(contMat), "Optimal contrasts\n.*")
})

test_that("summary.optContr summarizes and prints an optContr object", {
  contMat <- optContr(models, doses = c(0, 10), w = c(1, 1))
  expect_output(summary(contMat), "Optimal contrasts\n.*")
  expect_output(summary(contMat), "Contrast Correlation Matrix:.*")
})

test_that("plot.optContr plots contrast coefficients", {
  contMat <- optContr(models, doses = c(0, 10), w = c(1, 1))
  expect_silent(plot(contMat, plotType = "contrasts"))
  expect_silent(plot(contMat, plotType = "means"))
})

test_that("plotContr creates a ggplot object for the contrast coefficients", {
  contMat <- optContr(models, doses = c(0, 10), w = c(1, 1))
  expect_s3_class(plotContr(contMat), "ggplot")
})

test_that("plotContr creates a ggplot object with the correct data", {
  contMat <- optContr(models, doses = c(0, 10), w = c(1, 1))
  plot <- plotContr(contMat)
  
  # Ensure all dose levels are present in the plot
  expect_true(all(levels(as.factor(plot$data$dose)) %in% c(0, 10)))
  # Ensure all models are present in the plot
  expect_true(all(levels(as.factor(plot$data$model)) %in% c("linear")))
  # Check y-axis label
  expect_equal(plot$labels$y, "Contrast coefficients")
  # Check x-axis label
  expect_equal(plot$labels$x, "Dose")
})

test_that("lattice plot for optContr with superpose options works correctly", {
  contMat <- optContr(models, doses = c(0, 10), w = c(1, 1))
  expect_no_error(plot(contMat, plotType = "contrasts", superpose = TRUE))
})

test_that("lattice plot for optContr without superpose options works correctly", {
  contMat <- optContr(models, doses = c(0, 10), w = c(1, 1))
  expect_no_error(plot(contMat, plotType = "contrasts", superpose = FALSE))
})

# Additional test to ensure plotContr produces the correct ggplot2 plot
test_that("plotContr returns a ggplot2 plot with correct elements", {
  models <- Mods(linear = NULL, doses = c(0, 10, 25, 50, 100, 150))
  contMat <- optContr(models, doses = c(0, 10, 25, 50, 100, 150), w = rep(50, 6))
  p <- plotContr(contMat)
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "top")
})

# Additional test to ensure plot.optContr correctly sets y-axis labels
test_that("plot.optContr sets correct y-axis labels", {
  contMat <- optContr(models, doses = c(0, 10), w = c(1, 1))
  
  p1 <- plot(contMat, plotType = "contrasts", ylab = "Contrast coefficients")
  expect_equal(p1$ylab, "Contrast coefficients")
  
  p2 <- plot(contMat, plotType = "means", ylab = "Normalized model means")
  expect_equal(p2$ylab, "Normalized model means")
})
