#' Calculate Conditional or Predictive Power for Multiple Contrast Test
#'
#' @description
#' Calculates the predictive or conditional power for a multiple contrast test based on
#' interim data, e.g. for a futility interim analysis. This function can also be applied
#' to longitudinal endpoints, where at the time of interim analysis incomplete data is
#' available.
#'
#' @inheritParams powMCT
#' @param S_0t The covariance matrix for the first stage estimates
#' @param S_01 The covariance matrix anticipated for the estimates at
#'   study end
#' @param mu_0t The first stage estimates
#' @param type Whether predictive power (for a flat prior) or
#'   conditional power should be calculated. For conditional power
#'   mu_assumed needs to be specified.
#' @param mu_assumed Mean vector to assume for the second stage (only used when type is
#'   \samp{conditional}). If `NULL` (default), the first stage estimates `mu_0t` are used.
#' @param control A list specifying additional control parameters for the \samp{pmvnorm} calls in the code, see also
#' \samp{mvtnorm.control} for details.
#' @return Numeric containing the calculated power values
#' @seealso \code{\link{powMCT}} \code{\link{MCTtest}}, \code{\link{optContr}}
#' @references Bornkamp, B., Zhou, J., Xi, D. and Cao W. (2025). Futility analyses for the MCP-Mod methodology based
#' on longitudinal models, \emph{arXiv:2406.19965}
#' @examples
#'
#' doses <- c(0, 0.5, 1, 2, 4, 8)
#' mods <- Mods(emax = c(0.5, 1, 2, 4), sigEmax = rbind(c(0.5, 3), c(1, 3), c(2, 3), c(4, 3)), quadratic = -0.1, doses = doses)
#' w <- c(1, 0.5, 0.5, 0.5, 1, 1)
#' contMat <- optContr(models = mods, w = w)$contMat
#' sigma <- 0.3
#' n_final <- round(531 * w / sum(w))
#' n <- floor(n_final / 2)
#' S_0t <- diag(sigma^2 / n)
#' S_01 <- diag(sigma^2 / n_final)
#' ## assumed interim estimate
#' mu_0t <- 0.05 * doses / (doses + 1) + rnorm(6, 0, 0.382 / sqrt(n))
#' ## assumed mu (needed for conditional power)
#' mu_assumed <- 0.135 * doses / (doses + 1)
#' ## compare simulation based and numerical integration approach
#' powMCTInterim(contMat = contMat, S_0t = S_0t, S_01 = S_01, mu_0t = mu_0t, type = "predictive")
#' powMCTInterim(contMat = contMat, S_0t = S_0t, S_01 = S_01, mu_0t = mu_0t, type = "conditional", mu_assumed = mu_assumed)
#' powMCTInterim(contMat = contMat, S_0t = S_0t, S_01 = S_01, mu_0t = mu_0t, type = "predictive", alternative = "two.sided")
#' powMCTInterim(contMat = contMat, S_0t = S_0t, S_01 = S_01, mu_0t = mu_0t, type = "predictive", control = mvtnorm.control(maxpts = 1e5))
#' @export
powMCTInterim <- function(
  contMat,
  mu_0t,
  S_0t,
  S_01,
  alpha = 0.025,
  type = c("predictive", "conditional"),
  mu_assumed = NULL,
  alternative = c("one.sided", "two.sided"),
  control = mvtnorm.control()
) {
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  if (inherits(contMat, "optContr")) {
    contMat <- contMat$contMat
  }

  nDoses <- nrow(contMat)
  nContrasts <- ncol(contMat)

  if (!is.matrix(contMat)) {
    stop("contMat needs to be a matrix")
  }
  if ((nrow(S_0t) != ncol(S_0t)) || (nrow(S_01) != ncol(S_01))) {
    stop("S_0t and S_01 need to be square matrices")
  }
  if ((nrow(S_0t) != nDoses) || (nrow(S_01) != nDoses)) {
    stop(
      "Number of rows & cols of S_0t and S_01 need to match number of doses (i.e. number of rows of contMat)"
    )
  }
  if (type == "conditional") {
    if (is.null(mu_assumed)) {
      message("mu_assumed not supplied, setting mu_assumed = mu_0t")
      mu_assumed <- mu_0t
    }
    if (length(mu_assumed) != nDoses) {
      stop(
        "Length mu_assumed needs to match number of doses (i.e. number of rows of contMat)"
      )
    }
  }
  if ((length(mu_0t) != nDoses)) {
    stop(
      "Length of mu_0t needs to match number of doses (i.e. number of rows of contMat)"
    )
  }
  algorithm <- if (!inherits(control, "GenzBretz")) {
    if (!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    do.call(mvtnorm.control, control)
  } else {
    control
  }

  S_0t_inv <- solve(S_0t)
  St1_inv <- solve(S_01) - S_0t_inv
  St1 <- solve(St1_inv)

  # Pre-calculate the critical value.
  covMat <- t(contMat) %*% S_01 %*% contMat
  corMat <- cov2cor(covMat)
  criticalValue <- critVal(corMat, alpha = alpha, df = Inf)
  tTestDenominator <- sqrt(diag(covMat))

  P <- diag(1 / tTestDenominator)

  # Calculate the mean vector and covariance matrix for the predictive
  # or conditional distribution at the second stage.
  if (type == "predictive") {
    meanVector <- P %*% t(contMat) %*% mu_0t
    V0 <- S_0t + St1
    tmp <- P %*% t(contMat) %*% S_01 %*% St1_inv
    covMatrix <- tmp %*% V0 %*% t(tmp)
  }
  if (type == "conditional") {
    m0 <- S_01 %*% (S_0t_inv %*% mu_0t + St1_inv %*% mu_assumed)
    meanVector <- P %*% t(contMat) %*% m0
    tmp <- P %*% t(contMat) %*% S_01
    covMatrix <- tmp %*% St1_inv %*% t(tmp)
  }

  # Define integration boundaries.
  lower <- if (alternative == "two.sided") {
    rep(-criticalValue, nContrasts)
  } else {
    rep(-Inf, nContrasts)
  }
  upper <- rep(criticalValue, nContrasts)

  # Perform integration to obtain power value.
  intResult <- mvtnorm::pmvnorm(
    lower = lower,
    upper = upper,
    mean = as.numeric(meanVector),
    sigma = covMatrix,
    algorithm = algorithm
  )
  1 - intResult
}
