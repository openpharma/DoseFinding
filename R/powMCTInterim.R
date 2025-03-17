#' Calculate conditional or predictive power for multiple contrast test
#' 
#' Calculate predictive or conditional power for a multiple contrast test based on interim data, e.g. for a futility
#' interim analysis. This function can also be applied to longitudinal endpoints, where at the time of interim analysis
#' incomplete data is available.
#' 
#' 
#' @inheritParams powMCT 
#' @param S0t The covariance matrix for the first stage estimates
#' @param S_end The covariance matrix anticipated for the estimates at
#'   study end
#' @param mu0t The first stage estimates
#' @param type Whether predictive power (for a flat prior) or
#'   conditional power should be calculated. For conditional power
#'   mu_assumed needs to be specified.
#' @param mu_assumed Mean vector to assume for the second stage
#' @param control A list specifying additional control parameters for the \samp{pmvnorm} calls in the code, see also
#' \samp{mvtnorm.control} for details.
#' @return Numeric containing the calculated power values
#' @seealso \code{\link{powMCT}} \code{\link{MCTtest}}, \code{\link{optContr}}
#' @references Bornkamp, B., Zhou, J., Xi, D. and Cao W. (2025). Futility analyses for the MCP-Mod methodology based
#' on longitudinal models, \emph{arXiv:2406.19965}
#' @examples
#' 
#' doses <- c(0, 0.5, 1, 2, 4, 8)
#' mods <- Mods(emax = c(0.5,1,2,4), sigEmax = rbind(c(0.5, 3), c(1, 3), c(2, 3), c(4,3)), quadratic = -0.1, doses = doses)
#' w <- c(1,0.5,0.5,0.5,1,1)
#' contMat <- optContr(models=mods,w=w)$contMat
#' sigma <- 0.3
#' n_final <- round(531*w/sum(w))
#' n <- floor(n_final/2)
#' S0t <- diag(sigma^2/n)
#' S_end <- diag(sigma^2/n_final)
#' ## assumed interim estimate
#' mu0t <- 0.05*doses/(doses+1) + rnorm(6,0,0.382/sqrt(n))
#' ## assumed mu (needed for conditional power)
#' mu_assumed <- 0.135*doses/(doses+1)
#' ## compare simulation based and numerical integration approach
#' powMCTInterim(contMat = contMat, S0t = S0t, S_end = S_end, mu0t = mu0t, type = 'predictive')
#' powMCTInterim(contMat = contMat, S0t = S0t, S_end = S_end, mu0t = mu0t, type = 'conditional', mu_assumed = mu_assumed)

#' @export
powMCTInterim <- function(contMat, 
                          mu0t, S0t, S_end,  
                          alpha = 0.025,
                          type = c("predictive", "conditional"),
                          mu_assumed = NULL,
                          alternative = c("one.sided", "two.sided"),
                          control = mvtnorm.control()){
  
  alternative <- match.arg(alternative)
  if(inherits(contMat, "optContr"))
    contMat <- contMat$contMat
  
  nD <- nrow(contMat) # nr of doses
  nC <- ncol(contMat) # nr of contrasts
  
  if(!is.matrix(contMat))
    stop("contMat needs to be a matrix")
  if((nrow(S0t) != ncol(S0t)) | (nrow(S_end) != ncol(S_end)))
    stop("S0t and S_end need to be square matrices")
  if((nrow(S0t) != nD) | (nrow(S_end) != nD)) 
    stop("Number of rows & cols of S0t and S_end need to match number of doses (i.e. number of rows of contMat)")
  if(type == "conditional"){
    if(missing(mu_assumed)){
      message("mu_assumed not supplied, setting mu_assumed = mu0t")
      mu_assumed <- mu0t
    }
    if(length(mu_assumed) != nD)
      stop("Length mu_assumed needs to match number of doses (i.e. number of rows of contMat)")
  }
  if((length(mu0t) != nD)) 
    stop("Length of mu0t needs to match number of doses (i.e. number of rows of contMat)")

  S0t_inv <- solve(S0t)
  St1_inv <- solve(S_end) - S0t_inv
  St1 <- solve(St1_inv)
  
  ## pre-calculate critical value
  covMat <- t(contMat) %*% S_end %*% contMat
  corMat <- cov2cor(covMat)
  critV <- critVal(corMat, alpha = alpha, df = Inf)
  den <- sqrt(diag(covMat))  ## numerator of t-statistics
  
  P <- diag(1/den)
  ## simulate incremental information for second stage data
  if (type == "predictive") {
    mnV <- P%*%t(contMat)%*%mu0t
    V0 <- S0t + St1
    tmp <- P%*%t(contMat)%*%S_end%*%St1_inv
    V <- tmp%*%V0%*%t(tmp)
  }
  if (type == "conditional") {
    m0 <- S_end%*%(S0t_inv%*%mu0t + St1_inv%*%mu_assumed)
    mnV <- P%*%t(contMat)%*%m0
    tmp <- P%*%t(contMat)%*%S_end
    V <- tmp%*%St1_inv%*%t(tmp)
  }
  
  if(alternative[1] == "two.sided"){
    lower <- rep(-critV, nC)
  } else {
    lower <- rep(-Inf, nC)
  }
  upper <- rep(critV, nC)
  if (!missing(control)) {
    if(!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  
  pmvnormCall <- c(list(lower, upper,  mean = as.numeric(mnV), sigma = V,  algorithm = ctrl))
  1 - do.call(mvtnorm::pmvnorm, pmvnormCall)
}