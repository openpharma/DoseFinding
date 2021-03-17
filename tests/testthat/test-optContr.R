context("Optimal Contrasts")

# TODO
# * some test cases for constant shapes (what is expected?)
# * solnp gives other results than DoseFinding, quadratic programming or enumeration

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
  q <- apply(ncps, 1, quantile)
  expect_equal(q[, 5], q[, 2])
  expect_equal(q[, 5], q[, 3])
  expect_equal(q[, 5], q[, 4])
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
  expect_equal(cont_mat(0.05, TRUE, "u"), cont_mat(0.05, TRUE, "c"))
  expect_equal(cont_mat(c(0.05, 0.5), TRUE, "u"), cont_mat(c(0.05, 0.5), TRUE, "c"))
  expect_equal(cont_mat(c(0, 0.05), FALSE, "u"), cont_mat(c(0, 0.05), FALSE, "c"))
  # FIXME: what do we expect here?
  # expect_equal(cont_mat(c(0, 0.05, 0.5), FALSE, "u"), cont_mat(c(0, 0.05, 0.5), FALSE, "c"))
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
  # FIXME: what about these?
  ## optContr(modlist2, w=1, doses=c(0.05, 0.5), placAdj=TRUE, type = "u")
  ## optContr(modlist2, w=1, doses=c(0.05, 0.5), placAdj=TRUE, type = "c")
  ## optContr(modlist2, w=1, doses=c(0, 0.05, 0.5), placAdj=FALSE, type = "u")
  ## optContr(modlist2, w=1, doses=c(0, 0.05, 0.5), placAdj=FALSE, type = "c")
})
