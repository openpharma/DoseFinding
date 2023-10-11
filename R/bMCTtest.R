#' bMCTtest - Performs Bayesian multiple contrast test for normal and general outcomes
#'
#' @param dose Either vectors of dose values, or names of dose variable in the data frame specified in data.
#' @param resp Either vectors of response values, or names of response variable in the data frame specified in data.
#' @param data Data frame containing the variables referenced in dose and resp if data is not specified it is assumed that dose and resp are variables referenced from data (and no vectors)
#' @param models An object of class Mods, see Mods for details
#' @param S The covariance matrix of resp when type = "general", see Description.
#' @param type Determines whether inference is based on an ANCOVA model under a homoscedastic normality assumption (when type = "normal"), or estimates at the doses and their covariance matrix are specified directly in resp, S. See also fitMod and Pinheiro et al. (2014).
#' @param prior list of length equal to the number of doses (including plc) giving priors for each arm. Each element needs to be of class "normMix"
#' @param critV Supply a critical value for the maximum posterior probability of the contrasts being greater (or less) than zero. If this argument is NULL, this will be calculated
#' based on frequentist critical values.
#' @param alpha Significance level for the frequentist multiple contrast test, that is used to derive critical values for Bayesian decision rule, if none supplied via critV.
#' @param na.action A function which indicates what should happen when the data contain NAs.
#' @param mvtcontrol A list specifying additional control parameters for the qmvt and pmvt calls in the code, which are used to obtain critical values from frequentist MCP-Mod see also mvtnorm.control for details.
#' @param contMat Contrast matrix to apply to the ANCOVA dose-response estimates. The contrasts need to be in the columns of the matrix (i.e. the column sums need to be 0). If no contrast matrix is supplied
#' optimal contrasts will be calculated as for MCTtest, ignoring any prior information.  
#'
#' @return An object of class "bMCTtest", a list providing the results of the tests.
#' @export
#'
#' @examples
bMCTtest <- function (dose, resp, data = NULL, models, S = NULL, type = c("normal", "general"), 
                      prior, alpha = 0.025, na.action = na.fail, mvtcontrol = mvtnorm.control(),
                      contMat = NULL, critV = NULL) 
{
  type <- match.arg(type)
  alternative <- match.arg(alternative)
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
  
  covMat <- t(contMat) %*% vc %*% contMat
  corMat <- cov2cor(covMat)
  
  ## calculate frequentist critical values if none supplied
  if(is.null(critV)){
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
  
  if (alternative == "greater") {
    dec_prob <- pnorm(tStat) %*% unlist(post_res[[1]])
  }
  else{
    dec_prob <- pnorm(-tStat) %*% unlist(post_res[[1]])
  }
  ## maxprob <- max(dec_prob)
  
  res <- list(contMat = contMat, corMat = corMat, tStat = tStat, 
              alpha = alpha, alternative = alternative[1],
              critVal = 1 -critV,
              posterior = post_res)
  attr(res$tStat, "pVal") <- dec_prob
  class(res) <- "bMCTtest"
  res
}

checkAnalyArgs_bMCP <- function (dose, resp, data, S, type, prior, na.action, cal) 
{
  
  if (!is.null(data)) {
    if (!is.data.frame(data)) 
      stop("data argument needs to be a data frame")
    nams <- c(cal[2], cal[3])
    ind <- match(nams, names(data))
    if (any(is.na(ind))) 
      stop("variable(s): ", paste(nams[is.na(ind)], collapse = ", "), 
           " not found in ", cal[4])
    dd <- na.action(data[, nams])
  }
  else {
    if (!(is.numeric(resp) && is.null(dim(resp)))) {
      warning(cal[3], " is not a numeric but a ", class(resp)[1], 
              ", converting with as.numeric()")
      resp <- as.numeric(resp)
    }
    if (length(dose) != length(resp)) 
      stop(cal[2], " and ", cal[3], " not of equal length")
    dd <- na.action(data.frame(dose, resp))
    cal[2:3] <- gsub("\\$", "", cal[2:3])
    cal[2:3] <- gsub("\\[|\\]", "", cal[2:3])
    colnames(dd) <- cal[2:3]
  }
  doseNam <- cal[2]
  respNam <- cal[3]
  if (any(dd[[doseNam]] < -.Machine$double.eps)) 
    stop("dose values need to be non-negative")
  if (!is.numeric(dd[[doseNam]])) 
    stop("dose variable needs to be numeric")
  if (!is.numeric(dd[[respNam]])) 
    stop("response variable needs to be numeric")
  if (type == "general" & is.null(S)) 
    stop("S argument missing")
  if (type == "normal" & !is.null(S)) 
    message("Message: S argument ignored for type == \"normal\"\n")
  if (!is.null(S)) {
    if (!is.matrix(S)) 
      stop("S needs to be of class matrix")
    nD <- length(dd[[doseNam]])
    if (nrow(S) != nD | ncol(S) != nD) 
      stop("S and dose have non-conforming size")
  }
  if (length(unique(dd[[doseNam]])) != length(prior)) 
    stop("Dose and prior have non-conforming size")
  if (!all(unlist(lapply(prior, function(x) "normMix" %in% class(x))))) 
    stop("Prior needs to be of class normMix")
  
  ord <- order(dd[[doseNam]])
  dd <- dd[ord, ]
  Sout <- NULL
  if (type == "general") 
    Sout <- S[ord, ord]
  return(list(dd = dd, type = type, S = Sout, ord = ord, doseNam = doseNam, 
              respNam = respNam))
}

print.bMCTtest <- function(x, digits = 3, eps = 1e-3, ...){
  cat("Bayesian MCP-Mod\n")
  cat("\n","Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Contrast Correlation:","\n", sep="")
  print(round(x$corMat, digits))
  cat("\n","Posterior Mixture Weights:","\n",sep="")
  w <- round(unlist(x$posterior[[1]]), digits = digits)
  names(w) <- paste("Comp.", 1:length(w))
  print(w)
  cat("\n","Multiple Contrasts:","\n",sep="")
  ord <- rev(order(attr(x$tStat, "pVal")))
  pval <- format.pval(attr(x$tStat, "pVal"),
                      digits = digits, eps = eps)
  dfrm <- data.frame(round(x$tStat, digits)[ord, ],
                     pval[ord])
  names(dfrm) <- c(paste0("t-Stat (Comp. ", 1:ncol(x$tStat), ")"), "posterior probability")

  print(dfrm)
  if(!is.null(x$critVal)){
    cat("\n","Critical value: ", round(1- x$critVal, digits), sep="")
    if(attr(x$critVal, "Calc")){
      cat(" (alpha = ", x$alpha,", one-sided) \n", sep="")
    } else {
      cat("\n")
    }
  }
}

#' mvpostmix - calculate conjugate posterior mixture of multivariate normals with known covariance matrix 
#' (DeGroot 1970, Bernardo and Smith 1994)
#'
#' @param priormix prior multivariate normal mixture given as a list of length 3, providings weights, mean vectors and covariance matrices
#' @param mu_hat estimated mean response for each dose
#' @param S_hat estimated covariance matrix
#'
#' @return returns posterior mixture distribution as a list with weights, mean vectors and covariance matrices
#' @export
#'
#' @examples
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
  for(i in 1:length(postmix))
    postmix[[i]] <- vector("list", length(lw))
  
  for(i in 1:length(lw)){
    lw[i] <- log(priormix[[1]][[i]]) + dmvnorm(mu_hat, priormix[[2]][[i]], SigmaPred[[i]], log = TRUE)
    postmix[[2]][[i]] <- solve(priorPrec[[i]] + dataPrec) %*% (priorPrec[[i]] %*% priormix[[2]][[i]] + dataPrec %*% mu_hat)
    postmix[[3]][[i]] <- solve(priorPrec[[i]] + dataPrec)
  }
  
  postmix[[1]] <- as.list(exp(lw - logSumExp(lw)))
  
  postmix
}
