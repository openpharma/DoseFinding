## here the multiple contrast test related functions

#' Performs multiple contrast test
#'
#' This function performs a multiple contrast test. The contrasts are either directly specified in \samp{contMat} or
#' optimal contrasts derived from the \samp{models} argument. The directionality of the data (i.e. whether an increase
#' or decrease in the response variable is beneficial is inferred from the \samp{models} object, see
#' [Mods()]).
#'
#' For \samp{type = "normal"} an ANCOVA model based on a homoscedastic normality assumption (with additive covariates
#' specified in \samp{addCovars}) is fitted.
#'
#' For \samp{type = "general"} it is assumed multivariate normally distributed estimates are specified in \samp{resp}
#' with covariance given by \samp{S}, and the contrast test statistic is calculated based on this assumption. Degrees of
#' freedom specified in \samp{df}.
#'
#' Integrals over the multivariate t and multivariate normal distribution are calculated using the \samp{mvtnorm}
#' package.
#'
#' @param dose,resp Either vectors of equal length specifying dose and response values, or names of variables in the
#'   data frame specified in \samp{data}.
#' @param data Data frame containing the variables referenced in dose and resp if \samp{data} is not specified it is
#'   assumed that \samp{dose} and \samp{resp} are variables referenced from data (and no vectors)
#' @param models An object of class \samp{Mods}, see [Mods()] for details
#' @param S The covariance matrix of \samp{resp} when \samp{type = "general"}, see Description.
#' @param type Determines whether inference is based on an ANCOVA model under a homoscedastic normality assumption (when
#'   \samp{type = "normal"}), or estimates at the doses and their covariance matrix and degrees of freedom are specified
#'   directly in \samp{resp}, \samp{S} and \samp{df}. See also [fitMod()] and \insertCite{pinheiro2014}{DoseFinding}.
#' @param addCovars Formula specifying additive linear covariates (for \samp{type = "normal"})
#' @param placAdj Logical, if true, it is assumed that placebo-adjusted
#'   estimates are specified in \samp{resp} (only possible for \samp{type =
#'   "general"}).
#' @param alpha Significance level for the multiple contrast test
#' @param df Specify the degrees of freedom to use in case \samp{type = "general"}.  If this argument is missing
#'   \samp{df = Inf} is used (which corresponds to the multivariate normal distribution).  For type = "normal" the
#'   degrees of freedom deduced from the AN(C)OVA fit are used and this argument is ignored.
#' @param critV Supply a pre-calculated critical value. If this argument is NULL, no critical value will be calculated
#'   and the test decision is based on the p-values. If \samp{critV = TRUE} the critical value will be calculated.
#' @param pVal Logical determining, whether p-values should be calculated.
#' @param alternative Character determining the alternative for the multiple contrast trend test.
#' @param na.action A function which indicates what should happen when the data contain NAs.
#' @param mvtcontrol A list specifying additional control parameters for the \samp{qmvt} and \samp{pmvt} calls in the
#'   code, see also [mvtnorm.control()] for details.
#' @param contMat Contrast matrix to apply to the ANCOVA dose-response estimates. The contrasts need to be in the
#'   columns of the matrix (i.e. the column sums need to be 0).
#' @return An object of class MCTtest, a list containing the output.
#' @author Bjoern Bornkamp
#' @seealso [powMCT()], [optContr()]
#' @references 
#' \insertRef{hothorn2008}{DoseFinding}
#'   
#' \insertRef{pinheiro2014}{DoseFinding}
#' 
#' @examples
#'
#' ## example without covariates
#' data(biom)
#' ## define shapes for which to calculate optimal contrasts
#' modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
#'                 linInt = c(0, 1, 1, 1), doses = c(0, 0.05, 0.2, 0.6, 1))
#' m1 <- MCTtest(dose, resp, biom, models=modlist)
#' ## now calculate critical value (but not p-values)
#' m2 <- MCTtest(dose, resp, biom, models=modlist, critV = TRUE, pVal = FALSE)
#' ## now hand over critical value
#' m3 <- MCTtest(dose, resp, biom, models=modlist, critV = 2.24)
#'
#' ## example with covariates
#' data(IBScovars)
#' modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
#'                 linInt = c(0, 1, 1, 1), doses = c(0, 1, 2, 3, 4))
#' MCTtest(dose, resp, IBScovars, models = modlist, addCovars = ~gender)
#'
#' ## example using general approach (fitted on placebo-adjusted scale)
#' ancMod <- lm(resp~factor(dose)+gender, data=IBScovars)
#' ## extract estimates and information to feed into MCTtest
#' drEst <- coef(ancMod)[2:5]
#' vc <- vcov(ancMod)[2:5, 2:5]
#' doses <- 1:4
#' MCTtest(doses, drEst, S = vc, models = modlist, placAdj = TRUE,
#'         type = "general", df = Inf)
#'
#' ## example with general alternatives handed over
#' data(biom)
#' ## calculate contrast matrix for the step-contrasts
#' ## represent them as linInt models
#' models <- Mods(linInt=rbind(c(1,1,1,1),
#'                             c(0,1,1,1),
#'                             c(0,0,1,1),
#'                             c(0,0,0,1)),
#'                 doses=c(0,0.05,0.2,0.6,1))
#' plot(models)
#' ## now calculate optimal contrasts for these means
#' ## use weights from actual sample sizes
#' weights <- as.numeric(table(biom$dose))
#' contMat <- optContr(models, w = weights)
#' ## plot contrasts
#' plot(contMat)
#' ## perform multiple contrast test
#' MCTtest(dose, resp, data=biom, contMat = contMat)
#'
#' ## example for using the Dunnett contrasts
#' ## Dunnett contrasts
#' doses <- sort(unique(biom$dose))
#' contMat <- rbind(-1, diag(4))
#' rownames(contMat) <- doses
#' colnames(contMat) <- paste("D", doses[-1], sep="")
#' MCTtest(dose, resp, data=biom, contMat = contMat)
#'
#' @export
MCTtest <- function(dose, resp, data = NULL, models, S = NULL,
                    type = c("normal", "general"),
                    addCovars = ~1, placAdj = FALSE, 
                    alpha = 0.025, df = NULL, critV = NULL, pVal = TRUE,
                    alternative = c("one.sided", "two.sided"),
                    na.action = na.fail, mvtcontrol = mvtnorm.control(),
                    contMat = NULL){
  ## perform multiple contrast test
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  ## check for valid arguments
  cal <- as.character(match.call())
  lst <- checkAnalyArgs(dose, resp, data, S, type,
                        addCovars, placAdj, na.action, cal)
  dd <- lst$dd;type <- lst$type;S <- lst$S
  doseNam <- lst$doseNam;respNam <- lst$respNam

  ## calculate optimal contrasts and test-statistics
  doses <- unique(dd[[doseNam]])
  k <- length(doses)
  if(type == "normal"){
    dd[, doseNam] <- as.factor(dd[, doseNam])
    form <- paste(respNam, "~", doseNam, "+", addCovars[2], "-1", sep="")
    lm.fit <- lm(as.formula(form), data = dd)
    est <- coef(lm.fit)[1:k]
    vc <- vcov(lm.fit)[1:k, 1:k]
    df <- lm.fit$df.residual
  } else {
    est <- dd[[respNam]]
    vc <- S
    if(is.null(df))
      df <- Inf
  }
  if(is.null(contMat)){ # calculate optimal contrasts
    contMat <- optContr(models, doses, S=vc, placAdj=placAdj)$contMat
    rownames(contMat) <- doses
  } else { # contrast matrix specified
    if(inherits(contMat, "optContr"))
      contMat <- contMat$contMat
    if(nrow(contMat) != length(est))
      stop("contMat of incorrect dimensions")
  }
  ct <- as.vector(est %*% contMat)
  covMat <- t(contMat) %*% vc %*% contMat
  den <- sqrt(diag(covMat))
  tStat <- ct/den
  
  if(alternative == "two.sided"){
    tStat <- abs(tStat)
  }
  corMat <- cov2cor(covMat)
  
  if(is.null(critV)){
    if(!pVal){
      stop("either p-values or critical value need to be calculated.")
    }
  } else if(is.logical(critV) & critV == TRUE){
    critV <- critVal(corMat, alpha, df, alternative, mvtcontrol)  
    attr(critV, "Calc") <- TRUE # determines whether cVal was calculated
  } else { 
    pVal <- FALSE # pvals are not calculated if critV is supplied
    attr(critV, "Calc") <- FALSE
  }
  if(pVal){
    pVals <- MCTpval(contMat, corMat, df, tStat,
                     alternative, mvtcontrol)
  }
  res <- list(contMat = contMat, corMat = corMat, tStat = tStat,
              alpha = alpha, alternative = alternative[1])
  if(pVal)
    attr(res$tStat, "pVal") <- pVals
  res$critVal <- critV
  class(res) <- "MCTtest"
  res
}

#' @export
print.MCTtest <- function(x, digits = 3, eps = 1e-3, ...){
  cat("Multiple Contrast Test\n")
  cat("\n","Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Contrast Correlation:","\n", sep="")
  print(round(x$corMat, digits))
  cat("\n","Multiple Contrast Test:","\n",sep="")
  ord <- rev(order(x$tStat))
  if(!any(is.null(attr(x$tStat, "pVal")))){
    pval <- format.pval(attr(x$tStat, "pVal"),
                        digits = digits, eps = eps)
    dfrm <- data.frame(round(x$tStat, digits)[ord],
                       pval[ord])
    names(dfrm) <- c("t-Stat", "adj-p")
  } else {
    dfrm <- data.frame(round(x$tStat, digits)[ord])
    names(dfrm) <- c("t-Stat")
  }
  print(dfrm)
  if(!is.null(x$critVal)){
    twoSide <- x$alternative == "two.sided"
    vec <- c(" one-sided)", " two-sided)")
    cat("\n","Critical value: ", round(x$critVal, digits), sep="")
    if(attr(x$critVal, "Calc")){
      cat(" (alpha = ", x$alpha,",", vec[twoSide+1], "\n", sep="")
    } else {
      cat("\n")
    }
  }
}



#' Calculate critical value for multiple contrast test
#'
#' Calculation of the critical value for a maximum contrast test. This is based on the equicoordinate quantile function
#' of the multivariate normal or t distribution as implemented in the `qmvt` function from the mvtnorm package.
#'
#' @inheritParams MCTtest
#' @param corMat Correlation matrix of contrasts
#' @param df Specify the degrees of freedom to use, if this argument is missing \samp{df = Inf} is used (which
#'   corresponds to the multivariate normal distribution).
#' @param control A list specifying additional control parameters for the \samp{qmvt} and \samp{pmvt} calls in the code,
#'   see also [mvtnorm.control()] for details.
#' @author Bjoern Bornkamp
#' @seealso [powMCT()], [optContr()], [MCTtest()]
#' @examples
#'
#' R <- matrix(c(1,0.5,0.5,1), nrow=2)
#' critVal(R, alpha = 0.05, df = 1)
#' critVal(R, alpha = 0.05, df = 20)
#' critVal(R, alpha = 0.05, df = Inf)
#'
#' @export
critVal <- function(corMat, alpha = 0.025, df = NULL,
                    alternative = c("one.sided", "two.sided"),
                    control = mvtnorm.control()){
  ## calculate critical value
  alternative <- match.arg(alternative)
  if(missing(corMat))
    stop("corMat needs to be specified")
  if(is.null(df))
    stop("degrees of freedom need to be specified")
  tail <- ifelse(alternative[1] == "two.sided",
                 "both.tails", "lower.tail")
  if (!missing(control)) {
    if(!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  if(!is.finite(df)) # normal case
    df <- 0
  qmvtCall <- c(list(1-alpha, tail = tail, df = df, corr = corMat,
                algorithm = ctrl, interval = ctrl$interval))
  do.call(mvtnorm::qmvt, qmvtCall)$quantile
}

#' Calculate multiplicity adjusted p-values for multiple contrast test
#'
#' Calculate multiplicity adjusted p-values for a maximum contrast test corresponding to a set of contrasts and given a
#' set of observed test statistics. This function is exported as it may be a useful building block and used in more
#' complex testing situations that are not covered by [MCTtest()]. Most users probably don't need to use this
#' function.
#'
#' @inheritParams critVal
#' @param contMat Contrast matrix to use. The individual contrasts should be saved in the columns of the matrix
#' @param df Degrees of freedom to use for calculation.
#' @param tStat Vector of contrast test statistics
#' @return Numeric containing the calculated p-values.
#' @author Bjoern Bornkamp
#' @seealso [MCTtest()], [optContr()]
#' @references \insertRef{pinheiro2006b}{DoseFinding}
#' @examples
#' data(biom)
#' ## define shapes for which to calculate optimal contrasts
#' modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
#'                 linInt = c(0, 1, 1, 1), doses = c(0, 0.05, 0.2, 0.6, 1))
#' contMat <- optContr(modlist, w=1)$contMat
#' ## calculate inputs needed for MCTpval
#' fit <- lm(resp~factor(dose)-1, data=biom)
#' est <- coef(fit)
#' vc <- vcov(fit)
#' ct <- as.vector(est %*% contMat)
#' covMat <- t(contMat) %*% vc %*% contMat
#' den <- sqrt(diag(covMat))
#' tStat <- ct/den
#' corMat <- cov2cor(t(contMat) %*% vc %*% contMat)
#' MCTpval(contMat, corMat, df=100-5, tStat)
#' ## compare to
#' test <- MCTtest(dose, resp, biom, models=modlist)
#' attr(test$tStat, "pVal")
#' @export
MCTpval <- function(contMat, corMat, df, tStat,
                    alternative = c("one.sided", "two.sided"),
                    control = mvtnorm.control()){
  ## function to calculate p-values
  nD <- nrow(contMat)
  nMod <- ncol(contMat)
  if(missing(corMat))
    stop("corMat needs to be specified")
  if(missing(df))
    stop("degrees of freedom need to be specified")
  if(length(tStat) != nMod)
    stop("tStat needs to have length equal to the number of models")
  alternative <- match.arg(alternative)
  ctrl <- mvtnorm.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if(!is.finite(df)) # normal case
    df <- 0
  lower <- switch(alternative[1],
                  one.sided = matrix(rep(-Inf, nMod^2), nrow = nMod),
                  two.sided = matrix(rep(-tStat, each = nMod), nrow = nMod))
  upper <- switch(alternative[1],
                  one.sided = matrix(rep(tStat, each = nMod), nrow = nMod),
                  two.sided = matrix(rep(tStat, each = nMod), nrow = nMod))
  pVals <- numeric(nMod)
  for(i in 1:nMod){
    tmp <- 1 - mvtnorm::pmvt(lower[,i], upper[,i], df = df,
                    corr = corMat, algorithm = ctrl)
    pVals[i] <- tmp
    if(attr(tmp,"msg") != "Normal Completion"){
      warning(sprintf("Warning from mvtnorm::pmvt: %s.", attr(tmp, "msg")))
      if(attr(tmp, "msg") == "Covariance matrix not positive semidefinite"){
        warning("Setting calculated p-value to NA")
        pVals[i] <- NA
      }
    }
  }
  pVals
}
