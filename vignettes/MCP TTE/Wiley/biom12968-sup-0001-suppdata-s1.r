DoseFinding:::critVal
critVal1 <- function (corMat, alpha = 0.025, df = NULL, alternative = c("one.sided", 
                                                                        "two.sided"), control = mvtnorm.control())  # same as critVal
{
  alternative <- match.arg(alternative)
  if (missing(corMat)) 
    stop("corMat needs to be specified")
  if (is.null(df)) 
    stop("degrees of freedom need to be specified")
  tail <- ifelse(alternative[1] == "two.sided", "both.tails", 
                 "lower.tail")
  if (!missing(control)) {
    if (!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  }
  else {
    ctrl <- control
  }
  if (!is.finite(df)) 
    df <- 0
  qmvtCall <- c(list(1 - alpha, tail = tail, df = df, corr = corMat, 
                     algorithm = ctrl, interval = ctrl$interval))
  do.call("qmvt", qmvtCall)$quantile # Quantiles of the Multivariate t Distribution
}


powCalc1 = function (alternative, critV, df, corMat, deltaMat, control) 
{
  nC <- nrow(corMat)
  if (alternative[1] == "two.sided") {
    lower <- rep(-critV, nC)
  }
  else {
    lower <- rep(-Inf, nC)
  }
  upper <- rep(critV, nC)
  if (!missing(control)) {
    if (!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  }
  else {
    ctrl <- control
  }
  ctrl$interval <- NULL
  nScen <- ncol(deltaMat)
  res <- numeric(nScen)
  for (i in 1:nScen) {
    pmvtCall <- c(list(lower, upper, df = df, corr = corMat, 
                       delta = deltaMat[, i], algorithm = ctrl))
    res[i] <- as.vector(1 - do.call("pmvt", pmvtCall))
  }
  names(res) <- colnames(deltaMat)
  res
}

### S0 is under null to calculate critical value, S is under alternative

# mcpmod.power1=powMCT1(contrast, altModels = model1, alpha = 0.05, S=S[,,1], S0=S0, placAdj=T, df=Inf)
contMat <- contrast
powMCT1=function (contMat, alpha = 0.025, altModels, n, sigma, S, S0, placAdj = FALSE, # S0 is the only difference compared to powMCT
                  alternative = c("one.sided", "two.sided"), df, critV = TRUE, 
                  control = mvtnorm.control()) 
{
  alternative <- match.arg(alternative)
  if (inherits(contMat, "optContr")) {
    if (attr(contMat, "placAdj") != placAdj) {
      message("using \"placAdj\" specification from contMat object")
      placAdj <- attr(contMat, "placAdj")
    }
    contMat <- contMat$contMat
  }
  if (!is.matrix(contMat)) 
    stop("contMat needs to be a matrix")
  nD <- nrow(contMat)
  nC <- ncol(contMat)
  if (missing(S)) {
    if (missing(n) | missing(sigma)) 
      stop("Either S or n and sigma need to be specified")
    if (length(n) == 1) 
      n <- rep(n, nD)
    if (length(n) != nD) 
      stop("n needs to be of length nrow(contMat)")
    S <- sigma^2 * diag(1/n)
    df <- sum(n) - nD
  }
  else {
    if (!missing(n) | !missing(sigma)) 
      stop("Need to specify exactly one of \"S\" or \"n\" and \"sigma\"")
    if (nrow(S) != ncol(S)) 
      stop("S needs to be a square matrix")
    if (nrow(S) != nD) 
      stop("S needs to have as many rows&cols as there are doses")
    if (missing(df)) 
      stop("need to specify degrees of freedom in \"df\", when specifying \"S\"")
  }
  if (missing(altModels)) 
    stop("altModels argument needs to be specified")
  muMat <- getResp(altModels) # START!!!!
  if (placAdj) {
    muMat <- sweep(muMat, 2, muMat[1, ], "-")
    muMat <- muMat[-1, , drop = FALSE]
  }
  if (nrow(muMat) != nD) 
    stop("Incompatible contMat and muMat")
  if (missing(S)) {
    if (missing(df)) 
      stop("degrees of freedom need to be specified in df")
    df <- sum(n) - nD
  }
  deltaMat <- t(contMat) %*% muMat
  #print(list(contMat=contMat, muMat=muMat, deltaMat=deltaMat))
  covMat <- t(contMat) %*% S %*% contMat
  den <- sqrt(diag(covMat))
  deltaMat <- deltaMat/den
  #print(list(covMat=covMat, den=den, deltaMat=deltaMat))
  if (alternative == "two.sided") {
    deltaMat <- abs(deltaMat)
  }
  corMat <- cov2cor(covMat)
  covMat0 <- t(contMat) %*% S0 %*% contMat      #!!!!!!!!! NEW !!!!!!!!!!
  corMat0 <- cov2cor(covMat0)      #!!!!!!!!! NEW !!!!!!!!!!

  if (!is.finite(df)) 
    df <- 0
  if (is.logical(critV) & critV == TRUE) {
    critV <- critVal1(corMat0, alpha, df, alternative, control)      # corMat0 instead of corMat & critVal1 instead of critVal !!!!!!!!
  }
  res <- powCalc1(alternative, critV, df, corMat, deltaMat, # critV is from S0, corMat & deltaMat are from S
                  control)
  #print(list(critV=critV, df=df, corMat=corMat, corMat0=corMat0, deltaMat=deltaMat, control=control))
  res}


 contMat <- contrast
powMCT2=function (contMat, alpha = 0.025, altModels, n, sigma, S, S0, placAdj = FALSE, # S0 is the only difference compared to powMCT
                  alternative = c("one.sided", "two.sided"), df, critV = TRUE, 
                  control = mvtnorm.control()) 
{
  alternative <- match.arg(alternative)
  if (inherits(contMat, "optContr")) {
    if (attr(contMat, "placAdj") != placAdj) {
      message("using \"placAdj\" specification from contMat object")
      placAdj <- attr(contMat, "placAdj")
    }
    contMat <- contMat$contMat
  }
  if (!is.matrix(contMat)) 
    stop("contMat needs to be a matrix")
  nD <- nrow(contMat) ########### HERE
  nC <- ncol(contMat) ########### HERE
  if (missing(S)) {
    if (missing(n) | missing(sigma)) 
      stop("Either S or n and sigma need to be specified")
    if (length(n) == 1) 
      n <- rep(n, nD)
    if (length(n) != nD) 
      stop("n needs to be of length nrow(contMat)")
    S <- sigma^2 * diag(1/n)
    df <- sum(n) - nD
  }
  else {
    if (!missing(n) | !missing(sigma)) 
      stop("Need to specify exactly one of \"S\" or \"n\" and \"sigma\"")
    if (nrow(S) != ncol(S)) 
      stop("S needs to be a square matrix")
    if (nrow(S) != nD) 
      stop("S needs to have as many rows&cols as there are doses")
    if (missing(df)) 
      stop("need to specify degrees of freedom in \"df\", when specifying \"S\"")
  }
  if (missing(altModels)) 
    stop("altModels argument needs to be specified")
  muMat <- getResp(altModels) ########### HERE
  if (placAdj) {
    muMat <- sweep(muMat, 2, muMat[1, ], "-") ########### HERE
    muMat <- muMat[-1, , drop = FALSE]        ########### HERE
  }
  if (nrow(muMat) != nD) 
    stop("Incompatible contMat and muMat")
  if (missing(S)) {
    if (missing(df)) 
      stop("degrees of freedom need to be specified in df")
    df <- sum(n) - nD
  }
  deltaMat <- t(contMat) %*% muMat ########### HERE
  #print(list(contMat=contMat, muMat=muMat, deltaMat=deltaMat))
  covMat <- t(contMat) %*% S %*% contMat
  den <- sqrt(diag(covMat)) ########### HERE
  deltaMat <- deltaMat/den  ########### HERE
  #print(list(covMat=covMat, den=den, deltaMat=deltaMat))
  if (alternative == "two.sided") {
    deltaMat <- abs(deltaMat)
  }
  corMat <- cov2cor(covMat) ########### HERE
  covMat0 <- t(contMat) %*% S0 %*% contMat      #!!!!!!!!! NEW !!!!!!!!!!
  corMat0 <- cov2cor(covMat0)      #!!!!!!!!! NEW !!!!!!!!!!
  
  
  if (!is.finite(df)) 
    df <- 0
  if (is.logical(critV) & critV == TRUE) {
    critV <- critVal(corMat0, alpha, df, alternative, control)      # corMat0 instead of corMat & critVal1 instead of critVal !!!!!!!!
  }
  res <- DoseFinding:::powCalc(alternative, critV, df, corMat, deltaMat, # critV is from S0, corMat & deltaMat are from S
                  control)
  print(list(res = res, critV=critV, df=df, corMat=corMat, corMat0=corMat0, covMat=covMat, covMat0=covMat0, deltaMat=deltaMat))
  res}

