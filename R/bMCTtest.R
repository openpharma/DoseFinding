#' Performs Bayesian multiple contrast test
#'
#' This function performs a Bayesian multiple contrast test using normal mixture priors for the response on each dose,
#' as proposed in Fleischer et al. (2022). For a general description of the multiple contrast test see
#' \code{\link{MCTtest}}.
#'
#' If \samp{type = "normal"}, an ANCOVA model based on a homoscedastic normality assumption is fitted and posteriors for
#' dose-response and contrast vectors are obtained assuming a known variance.
#'
#' For \samp{type = "general"} it is assumed multivariate normally distributed estimates are specified in \samp{resp}
#' with covariance given by \samp{S}, which define the likelihood.  Posteriors for dose-response and contrast vectors
#' are then obtained assuming a known covariance matrix S
#'
#' The multiple contrast test decision is based on the maximum posterior probability of a contrast being greater than
#' zero. Thresholds for the posterior probability can either be supplied or will be derived from frequentist critical
#' values. In the latter case the Bayesian test will give approximately the same results as the frequentist multiple
#' contrast test if uninformative priors are used.
#'
#' For the default calculation of optimal contrasts the prior information is ignored (i.e. contrasts are calculated in
#' the same way as in \code{\link{MCTtest}}).  Fleischer et al. (2022) discuss using contrasts that take the prior
#' effective sample sizes into account, which can be slightly more favourable for the Bayesian MCT test. Such
#' alternative contrasts can be directly handed over via the \samp{contMat} argument.
#'
#' For analysis with covariate adjustment, covariate-adjusted \samp{resp} and \samp{S} can be supplied together with
#' using \samp{type = "general"}. See  `vignette("binary_data")` vignette "Design and analysis template MCP-Mod for binary data" for an example
#' on how to obtain covariate adjusted estimates.
#'
#' @inheritParams MCTtest
#' @param prior List of length equal to the number of doses with the prior for each arm.  Each element needs to be of
#'   class \samp{normMix} (See \samp{RBesT} package documentation). It is assumed that the i-th component of the prior
#'   list corresponds to the i-th largest dose. For example the first entry in the list is the prior for the placebo
#'   group, the second entry the prior for the second lowest dose and so on.  Internally the priors across the different
#'   arms are combined (densities multiplied) assuming independence. The resulting multivariate normal mixture prior
#'   will have as many components as the product of the number of components of the individual mixture priors. The
#'   posterior mixture is part of the result object under "posterior".
#' @param alpha Significance level for the frequentist multiple contrast test. If no critical values are supplied via
#'   \samp{critV} this is used to derive critical values for Bayesian decision rule.
#' @param contMat Contrast matrix to apply to the posterior dose-response estimates. The contrasts need to be in the
#'   columns of the matrix (i.e. the column sums need to be 0). If not specified optimal contrasts are calculated using
#'   \code{\link{optContr}}.
#' @param critV Supply a critical value for the maximum posterior probability of the contrasts being greater than zero
#'   that needs to be surpassed to establish a non-flat dose-response. If this argument is NULL, this will be derived
#'   from critical values for frequentist MCP-Mod using the provided \samp{alpha}.
#' @return An object of class bMCTtest, a list containing the output.
#' @author Marius Thomas
#' @export
#' @seealso \code{\link{MCTtest}}, \code{\link{optContr}}
#' @references Fleischer, F., Bossert, S., Deng, Q., Loley, C. and Gierse, J. (2022).  Bayesian MCP-Mod,
#'   \emph{Pharmaceutical Statistics}, \bold{21}, 654--670
#' @examples
#'
#'
#' if (require("RBesT")) {
#' 
#' ###############################
#' ## Normal outcome
#' ###############################
#'
#' data(biom)
#' ## define shapes for which to calculate optimal contrasts
#' doses <- c(0, 0.05, 0.2, 0.6, 1)
#' modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
#'                 linInt = c(0, 1, 1, 1), doses = doses)
#' ## specify an informative prior for placebo, weakly informative for other arms
#' plc_prior <- mixnorm(inf = c(0.8, 0.4, 0.1), rob = c(0.2, 0.4, 10))
#' vague_prior <- mixnorm(c(1, 0, 10))
#' ## i-th component of the prior list corresponds to the i-th largest dose
#' ## (e.g. 1st component -> placebo prior; last component prior for top dose)
#' prior <- list(plc_prior, vague_prior, vague_prior, vague_prior, vague_prior)
#'
#' m1 <- bMCTtest(dose, resp, biom, models=modlist, prior = prior)
#' ## now supply a critical value (= threshold for maxmimum posterior probability)
#' m2 <- bMCTtest(dose, resp, biom, models=modlist, prior = prior, critV = 0.99)
#'
#' ####################################
#' ## Binary outcome with covariates
#' ####################################
#'\dontrun{
#' ## generate data
#' logit <- function(p) log(p / (1 - p))
#' inv_logit <- function(y) 1 / (1 + exp(-y))
#' doses <- c(0, 0.5, 1.5, 2.5, 4)
#'
#' ## set seed and ensure reproducibility across R versions
#' set.seed(1, kind = "Mersenne-Twister", sample.kind = "Rejection", normal.kind = "Inversion")
#' group_size <- 100
#' dose_vector <- rep(doses, each = group_size)
#' N <- length(dose_vector)
#' ## generate covariates
#' x1 <- rnorm(N, 0, 1)
#' x2 <- factor(sample(c("A", "B"), N, replace = TRUE, prob = c(0.6, 0.4)))
#' ## assume approximately logit(10%) placebo and logit(35%) asymptotic response with ED50=0.5
#' prob <- inv_logit(emax(dose_vector, -2.2, 1.6, 0.5) + 0.3 * x1 + 0.3 * (x2 == "B"))
#' dat <- data.frame(y = rbinom(N, 1, prob),
#'                   dose = dose_vector, x1 = x1, x2 = x2)
#'
#' ## specify an informative prior for placebo (on logit scale), weakly informative for other arms
#' plc_prior <- mixnorm(inf = c(0.8, -2, 0.5), rob = c(0.2, -2, 10))
#' vague_prior <- mixnorm(c(1, 0, 10))
#' prior <- list(plc_prior, vague_prior, vague_prior, vague_prior, vague_prior)
#'
#' ## candidate models
#' mods <- Mods(emax = c(0.25, 1), sigEmax = rbind(c(1, 3), c(2.5, 4)), betaMod = c(1.1, 1.1),
#'              placEff = logit(0.1), maxEff = logit(0.35)-logit(0.1),
#'              doses = doses)
#'
#' fit_cov <- glm(y~factor(dose) + 0 + x1 + x2, data = dat, family = binomial)
#'
#' covariate_adjusted_estimates <- function(mu_hat, S_hat, formula_rhs,
#'                                          doses, other_covariates, n_sim) {
#'   ## predict every patient under *every* dose
#'   oc_rep <- as.data.frame(lapply(other_covariates, function(col) rep(col, times = length(doses))))
#'   d_rep <- rep(doses, each = nrow(other_covariates))
#'   pdat <- cbind(oc_rep, dose = d_rep)
#'   X <- model.matrix(formula_rhs, pdat)
#'   ## average on probability scale then backtransform to logit scale
#'   mu_star <- logit(tapply(inv_logit(X %*% mu_hat), pdat$dose, mean))
#'   ## estimate covariance matrix of mu_star
#'   pred <- replicate(n_sim, logit(tapply(inv_logit(X %*% drop(mvtnorm::rmvnorm(1, mu_hat, S_hat))),
#'                                         pdat$dose, mean)))
#'   return(list(mu_star = as.numeric(mu_star), S_star = cov(t(pred))))
#' }
#'
#' ca <- covariate_adjusted_estimates(coef(fit_cov), vcov(fit_cov), ~factor(dose)+0+x1+x2,
#'                                    doses, dat[, c("x1", "x2")], 1000)
#' bMCTtest(doses, ca$mu_star, S = ca$S_star, type = "general", models = mods, prior = prior)
#'}
#' ################################################
#' ## example with contrasts handed over
#' ################################################
#'
#' data(biom)
#' ## define shapes for which to calculate optimal contrasts
#' doses <- c(0, 0.05, 0.2, 0.6, 1)
#' modlist <- Mods(emax = 0.05, linear = NULL, sigEmax = c(0.5, 5),
#'                 linInt = c(0, 1, 1, 1), doses = doses)
#'
#' ## specify an informative prior for placebo, weakly informative for other arms
#' plc_prior <- mixnorm(inf = c(0.8, 0.4, 0.1), rob = c(0.2, 0.4, 10), sigma = 0.7)
#' vague_prior <- mixnorm(c(1, 0, 10), sigma = 0.7)
#' prior <- list(plc_prior, vague_prior, vague_prior, vague_prior, vague_prior)
#'
#' ## use prior effective sample sizes to calculate optimal contrasts
#' prior_ess <- unlist(lapply(prior, ess))
#' n_grp <- as.numeric(table(biom$dose))
#' weights <- n_grp + prior_ess
#' cmat <- optContr(modlist, w = weights)
#'
#' bMCTtest(dose, resp, biom, models=modlist, prior = prior, contMat = cmat)
#' }
#' 
bMCTtest <- function (dose, resp, data = NULL, models, S = NULL, type = c("normal", "general"), 
                      prior, alpha = 0.025, na.action = na.fail, mvtcontrol = mvtnorm.control(),
                      contMat = NULL, critV = NULL) 
{
  type <- match.arg(type)
  cal <- as.character(match.call())
  
  lst <- checkAnalyArgs_bMCP(dose, resp, data, S, type, prior, na.action, cal)
  dd <- lst$dd
  type <- lst$type
  S <- lst$S
  doseNam <- lst$doseNam
  respNam <- lst$respNam
  doses <- unique(dd[[doseNam]])
  k <- length(doses)
  
  if (type == "normal") {
    dd[, doseNam] <- as.factor(dd[, doseNam])
    form <- paste(respNam, "~", doseNam, "+", "-1", sep = "")
    lm.fit <- lm(as.formula(form), data = dd)
    est <- coef(lm.fit)[1:k]
    vc <- vcov(lm.fit)[1:k, 1:k]
  }
  else {
    est <- dd[[respNam]]
    vc <- S
  }
  
  if (is.null(contMat)) {
    contMat <- optContr(models, doses, S=vc)$contMat
    rownames(contMat) <- doses
  }
  else {
    if (inherits(contMat, "optContr")) 
      contMat <- contMat$contMat
    if (nrow(contMat) != length(est)) 
      stop("contMat of incorrect dimensions")
  }
  
  ## calculate frequentist critical values if none supplied
  if(is.null(critV)){
    covMat <- t(contMat) %*% vc %*% contMat
    corMat <- cov2cor(covMat)
    critV <- critVal(corMat, alpha, df = Inf, alternative = "one.sided", mvtcontrol) ## using df = INF so values are derived using multivariate normal
    critV <- pnorm(critV)
    attr(critV, "Calc") <- TRUE
  }
  else{
    attr(critV, "Calc") <- FALSE
  }
  
  ## write complete multivariate normal specification for the prior
  n_comps <- unlist(lapply(prior, ncol))
  args <- lapply(1:k, function(x) 1:n_comps[x])
  comp_ind <- do.call("expand.grid", args) 
  
  n_comps_prior <- nrow(comp_ind)
  
  prior_weight <- matrix(sapply(1:k, function(x) sapply(1:n_comps_prior, function(y) prior[[x]][1, comp_ind[y,x]])), nrow = n_comps_prior)
  prior_weight <- apply(prior_weight, 1, prod)
  prior_mean <- matrix(sapply(1:k, function(x) sapply(1:n_comps_prior, function(y) prior[[x]][2, comp_ind[y,x]])), nrow = n_comps_prior)
  prior_sd <- matrix(sapply(1:k, function(x) sapply(1:n_comps_prior, function(y) prior[[x]][3, comp_ind[y,x]])), nrow = n_comps_prior)
  
  prior_weight <- as.list(prior_weight)
  prior_mean <- asplit(prior_mean, 1)
  prior_vc <- lapply(asplit(prior_sd^2, 1), diag)
  prior_mix <- list(prior_weight, prior_mean, prior_vc)
  
  ## Bayesian conjugate posterior
  post_res <- mvpostmix(prior_mix, est, vc)
  
  mu_mat <- do.call(rbind, lapply(post_res[[2]], as.numeric))
  ct <- t(contMat) %*% t(mu_mat) ## contrasts for each component (one candidate model per row)
  den <- lapply(post_res[[3]], function(x) t(contMat) %*% x %*% contMat)
  den <- sqrt(do.call(cbind, lapply(den, diag)))
  tStat <- ct/den
  
  dec_prob <- pnorm(tStat) %*% unlist(post_res[[1]])
  
  res <- list(contMat = contMat, tStat = tStat, 
              alpha = alpha, critVal = 1 - critV,
              posterior = post_res)
  attr(res$tStat, "pVal") <- dec_prob
  class(res) <- "bMCTtest"
  res
}

#' @export
print.bMCTtest <- function(x, digits = 3, eps = 1e-3, ...){
  cat("Bayesian MCP-Mod\n")
  cat("\n","Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Posterior Mixture Weights:","\n",sep="")
  w <- round(unlist(x$posterior[[1]]), digits = digits)
  names(w) <- paste("Comp.", 1:length(w))
  print(w)
  ord <- rev(order(attr(x$tStat, "pVal")))
  pval <- format.pval(attr(x$tStat, "pVal"),
                      digits = digits, eps = eps)
  dfrm <- data.frame(round(x$tStat, digits)[ord, , drop = FALSE],
                     pval[ord])
  names(dfrm) <- c(paste0("Comp. ", 1:ncol(x$tStat)), "posterior probability")
  cat("\n","Bayesian t-statistics:","\n",sep="")
  print(dfrm)
  if(!is.null(x$critVal)){
    cat("\n","Critical value (for maximum posterior probability): ", round(1- x$critVal, digits), sep="")
    if(attr(x$critVal, "Calc")){
      cat(" (alpha = ", x$alpha,", one-sided) \n", sep="")
    } else {
      cat("\n")
    }
  }
}



#' Prior to posterior updating for a multivariate normal mixture
#'
#' Calculate conjugate posterior mixture of multivariate normals with known covariance matrix
#'
#'
#' @param priormix Prior multivariate normal mixture given as a list of length 3. The first list entry contains the
#'   mixture weights, the second component the mean vectors and the third component of the list the covariance matrices.
#' @param mu_hat estimated mean response for each dose
#' @param S_hat estimated covariance matrix
#' @return Returns a posterior multivariate normal mixture as a list of length 3, containing mixture weights, mean
#'   vectors and covariance matrices.
#' @author Marius Thomas
#' @references Bernardo, J. M., and Smith, A. F. (1994). Bayesian theory. John Wiley & Sons.
#' @export
mvpostmix <- function(priormix, mu_hat, S_hat)
{

  logSumExp <- function(lx){
    lm <- max(lx)
    lm + log(sum(exp(lx - lm)))
  }
  
  dataPrec <- solve(S_hat)
  priorPrec <- lapply(priormix[[3]], solve)
  postPrec <- lapply(priorPrec, function(x) x + dataPrec)
  SigmaPred <- lapply(priormix[[3]], function(x) x + S_hat)
  
  lw <- numeric(length(priormix[[1]]))
  postmix <- vector("list", 3)
  names(postmix) <- c("weights", "mean", "covmat")
  for(i in 1:3)
    postmix[[i]] <- vector("list", length(lw))
  
  ## The posterior distribution is a mixture of multivariate normals with updated mixture weights.
  ## Posterior weights are updated based on the prior predictive (marginal) probabilities of the data under each
  ## component of the mixture. 
  ## In the case of a MVN likelihood with known covariance and MVN priors for the mean the 
  ## prior predictive distributions are MVN distribution with mean vectors equal to the prior components' mean vectors 
  ## and covariance matrices which are the sum of the prior components' covariance matrices and the "known" covariance 
  ## matrix of the data (for which S_hat is plugged in here)
  for(i in 1:length(lw)){
    lw[i] <- log(priormix[[1]][[i]]) + mvtnorm::dmvnorm(mu_hat, priormix[[2]][[i]], SigmaPred[[i]], log = TRUE)
    postmix[[2]][[i]] <- solve(priorPrec[[i]] + dataPrec) %*% (priorPrec[[i]] %*% priormix[[2]][[i]] + dataPrec %*% mu_hat)
    postmix[[3]][[i]] <- solve(priorPrec[[i]] + dataPrec)
  }
  postmix[[1]] <- as.list(exp(lw - logSumExp(lw)))
  
  for(i in 1:3)
    names(postmix[[i]]) <- paste0("Comp", 1:length(lw))

   postmix
}
