## all design related functions for power calculations

mvtnorm.control <- function(maxpts = 30000, abseps = 0.001,
                            releps = 0, interval = NULL){
  res <- list(maxpts = maxpts, abseps = abseps,
              releps = releps, interval = interval)
  class(res) <- "GenzBretz"    
  res
}

powCalc <- function(alternative, critV, df, corMat, deltaMat, control){
  nC <- nrow(corMat) # number of contrasts
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
  ctrl$interval <- NULL      # not used with pmvt
  nScen <- ncol(deltaMat)
  res <- numeric(nScen)
  for(i in 1:nScen){
    pmvtCall <- c(list(lower, upper, df = df, corr = corMat, delta = deltaMat[,i],
                       algorithm = ctrl))
    res[i] <- as.vector(1 - do.call("pmvt", pmvtCall))
  }
  names(res) <- colnames(deltaMat)
  res
}

powMCT <- function(contMat, alpha = 0.025, altModels,
                   n, sigma, S, placAdj = FALSE,
                   alternative = c("one.sided", "two.sided"),
                   df, critV = TRUE,
                   control = mvtnorm.control()){
  alternative <- match.arg(alternative)
  if(inherits(contMat, "optContr")){
    if(attr(contMat, "placAdj") != placAdj){
      message("using \"placAdj\" specification from contMat object")
      placAdj <- attr(contMat, "placAdj")
    }
    contMat <- contMat$contMat
  }
  if(!is.matrix(contMat))
    stop("contMat needs to be a matrix")
  nD <- nrow(contMat) # nr of doses
  nC <- ncol(contMat) # nr of contrasts
  ## extract covariance matrix
  if(missing(S)){
    if(missing(n) | missing(sigma))
      stop("Either S or both n and sigma need to be specified")
    if(length(n) == 1)
      n <- rep(n, nD)
    if(length(n) != nD)
      stop("n needs to be of length nrow(contMat)")
    S <- sigma^2*diag(1/n)
    df <- sum(n) - nD
    if(df == 0)
      stop("cannot compute power: specified \"n\" and dose vector result in df = 0")
  } else {
    if(!missing(n)|!missing(sigma))
      stop("Need to specify either \"S\" or both \"n\" and \"sigma\"")
    if(nrow(S) != ncol(S))
      stop("S needs to be a square matrix")
    if(nrow(S) != nD)
      stop("S needs to have as many rows&cols as there are doses")
    if(missing(df))
      stop("need to specify degrees of freedom in \"df\", when specifying \"S\"")
  }
  ## extract means under the alternative
  if(missing(altModels))
    stop("altModels argument needs to be specified")
  muMat <- getResp(altModels)
  if(placAdj){
    muMat <- sweep(muMat, 2, muMat[1,], "-") # remove placebo column
    muMat <- muMat[-1, , drop=FALSE]
  }
  if(nrow(muMat) != nD)
    stop("Incompatible contMat and muMat")
  ## calculate non-centrality parameter
  deltaMat <- t(contMat) %*% muMat
  covMat <- t(contMat) %*% S %*% contMat
  den <- sqrt(diag(covMat))
  deltaMat <- deltaMat/den
  if(alternative == "two.sided"){
    deltaMat <- abs(deltaMat)
  }
  corMat <- cov2cor(covMat)
  
  if(!is.finite(df))
    df <- 0
  ## calculate critical value
  if(is.logical(critV) & critV == TRUE){
    critV <- critVal(corMat, alpha, df, alternative, control)
  } # else assume critV already contains critical value
  res <- powCalc(alternative, critV, df, corMat, deltaMat, control)
  ## class(res) <-  "powMCT"
  ## attr(res, "type") <- ifelse(missing(n), "S", "n&sigma")
  ## attr(res, "contMat") <- contMat
  ## attr(res, "muMat") <- muMat
  res
}

## print.powMCT <- function(x, ...){
##   attributes(x)[2:5] <- NULL
##   print(x)
## }

#' Function for power calculation for the multiple contrast test for
#' binary data (candidate models formulated on logit link scale) and
#' count data (negative binomial distribution with log link is
#' assumed; candidate models formulated on log mean scale)- Note this
#' function is somewhat limited (e.g. only one-sided testing allowed)
#' and not user-friendly (e.g. cannot hand over Mods objects). It is
#' not exported.
#'
#' @param n Vector of sample sizes per dose group
#' @param doses
#' @param candModList List containing the candidate models
#' @param respModList List containing the response models to assume
#' @param placEffu placebo effect on untransformed scale
#' @param maxEffu maximum treatment effect vs placebo on untransformed
#'   scale
#' @param type either "binary_logit" or "negative_binomial"
#' @param option whether the optimal contrast should be recalculated
#'   for each assumed true model (option "A") or not option "B"
#' @param alpha Significance level (one-sided testing is assumed)
#' @param theta Overdispersion parameter required for the negative
#'   binomial model. Parameterization of negative binomial: E(Y)=mu,
#'   Var(Y)=mu*(1+mu/theta)
#' @param control Control parameter for mvtnorm::qmvt, mvtnorm::pmvt
#'   (see ?mvtnorm.control)
#' @param contMat Contrast matrix (if non-model-based contrasts should
#'   be used). If a user-defined contrast matrix the argument option
#'   is ignored (automatically set to "B").
#' @param addArgs additional arguments for betaMod and linlog model
#'   (passed to Mods function)
#' 
#' @return Vector of calculated power values
#' @examples
#'mvt_control <- mvtnorm.control()
#'candModList <- list(emax = c(0.25, 1), sigEmax = rbind(c(1, 3), c(2.5, 4)), betaMod = c(1.1, 1.1))
#'powMCTBinCount(rep(20,5), doses = c(0, 0.5, 1.5, 2.5, 4),
#'               candModList=candModList, respModList=candModList,
#'               placEffu = 0.1, maxEffu = 0.25, 
#'               type = "binary_logit", option = "A",
#'               alpha = 0.1, theta, control = mvt_control,
#'               addArgs = list(scal = 4.8))
powMCTBinCount <- function(n, doses, candModList = NULL, respModList,
                           placEffu, maxEffu, 
                           type = c("binary_logit", "negative_binomial"),
                           option = c("A", "B"),
                           alpha, theta = NULL,
                           control, contMat = NULL, addArgs){
  type <- match.arg(type)
  stopifnot(length(doses) == length(n))
  if(is.null(candModList) & is.null(contMat))
    stop("either candModList or contMat need to be non-NULL")
  if(type == "binary_logit"){ 
    logit <- function(p)
      log(p/(1-p))
    trafo <- logit
  }
  if(type == "negative_binomial"){
    if(is.null(theta))
      stop("need argument theta for type = \"negative_binomial\"")
    trafo <- log
  }
  placEff_tr <- trafo(placEffu)
  maxEff_tr <- trafo(placEffu+maxEffu)-placEff_tr

  if(is.null(contMat)){
    mods <- do.call(Mods, append(candModList,
                                 list(placEff = placEff_tr, maxEff = maxEff_tr,
                                      doses = doses, addArgs=addArgs)))
    
    if(option == "B"){ # opt. contrasts not recalculated only calculated here
      cm <- optContr(mods, w=n)$contMat # assume diagonal cov matrix with entries 1/n_i
    }
  } else {
    option <- "B"
    cm <- contMat
  }

  ## calculate correlation matrix under null-hypothesis
  resp_null <- do.call(Mods, append(respModList,
                                    list(placEff = placEff_tr, maxEff = 0,
                                         doses = doses, addArgs=addArgs)))
  mu_null <- getResp(resp_null)[, 1, drop=FALSE] # under null all models are the same
  v_null <- getVarBinCount(mu_null, type, theta)
  S_null <- diag(as.vector(v_null)/n)
  
  resp_mods <- do.call(Mods, append(respModList,
                                    list(placEff = placEff_tr, maxEff = maxEff_tr,
                                         doses = doses, addArgs=addArgs)))
  resp <- getResp(resp_mods)
  nMod <- ncol(resp)
  pow <- numeric(nMod)
  for(i in 1:nMod){
    mu_vec <- resp[,i, drop=FALSE] # column i contains true response vector
    ## calculate covariance matrix
    v <- getVarBinCount(mu_vec, type, theta)
    S <- diag(as.vector(v)/n)
    if(option == "A"){ # (re)calculate optimal contrasts based on S
      contMat <- optContr(mods, S=S)
      cm <- contMat$contMat
    }
    
    ## calculate non-centrality parameter
    delta <- t(cm) %*% mu_vec
    covMat <- t(cm) %*% S %*% cm
    den <- sqrt(diag(covMat))
    delta <- delta/den
    corMat <- cov2cor(covMat)
    if(option == "A" | (option == "B" & i == 1)){ # for option B critV does not change
      if(option == "A"){ # ADDPLAN-DF appears to use corMat not corMat_null
        corMat_null <- corMat # for consistency do the same
      }
      if(option == "B"){
        covMat_null <- t(cm) %*% S_null %*% cm
        corMat_null <- cov2cor(covMat_null)
      }
      ## calculate critical value (df=0 corresponds to infinite degrees of freedom -> MVN distribution)
      qmvtCall <- c(list(1 - alpha, tail = "lower.tail", df = 0, corr = corMat_null,
                         algorithm = control, interval = control$interval))
      critV <- do.call("qmvt", qmvtCall)$quantile
    }
    ## calculate power
    lower <- rep(-Inf, ncol(corMat))
    upper <- rep(critV, ncol(corMat))
    control$interval <- NULL
    pmvtCall <- c(list(lower, upper, df = 0, corr = corMat,
                       delta = as.vector(delta), algorithm = control))
    pow[i] <- as.vector(1 - do.call("pmvt", pmvtCall))
  }
  names(pow) <- colnames(resp)
  pow
}

#' Function to calculate unit variances for mu_hat for binary (with
#' logit link) and negative-binomial (with log link). Resulting
#' variances need to be multiplied with factor 1/n.
#' @param muVec Vector of group means on transformed scale (logit
#'   scale for binary data, log scale for negative binomial)
#' @param type either "binary_logit" or "negative_binomial"
#' @param theta Overdispersion parameter required for the negative
#'   binomial model. Parameterization of negative binomial: E(Y)=mu,
#'   Var(Y)=mu*(1+mu/theta)
#' @return Vector of variances
getVarBinCount <- function(muVec, type, theta){
  if(type == "binary_logit"){
    inv_logit <- function(x)
      1/(1+exp(-x))
    p <- inv_logit(muVec) # mean responses on probability scale
    return(1 / (p * (1 - p)))
  }
  if(type == "negative_binomial")
    return((theta+exp(muVec))/(theta*exp(muVec)))
}
