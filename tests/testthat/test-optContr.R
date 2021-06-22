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
