## all design related functions for power calculations

#' Control options for pmvt and qmvt functions
#'
#' Returns a list (an object of class "GenzBretz") with control parameters for the \samp{pmvt} and \samp{qmvt} functions
#' from the \samp{mvtnorm} package. Note that the DoseFinding package always uses "GenzBretz" algorithm. See the mvtnorm
#' documentation for more information.
#'
#' @name mvtnorm-control
#' @param maxpts Maximum number of function values as integer.
#' @param abseps Absolute error tolerance as double.
#' @param releps Relative error tolerance as double.
#' @param interval Interval to be searched, when the quantile is calculated.
#' @export
mvtnorm.control <- function(maxpts = 30000, abseps = 0.001,
                            releps = 0, interval = NULL){
  res <- list(maxpts = maxpts, abseps = abseps,
              releps = releps, interval = interval)
  class(res) <- "GenzBretz"    
  res
}

#' Calculate power for multiple contrast test
#' 
#' Calculate power for a multiple contrast test for a set of specified
#' alternatives.
#' 
#' 
#' @param contMat Contrast matrix to use. The individual contrasts should be
#' saved in the columns of the matrix
#' @param alpha Significance level to use
#' @param altModels An object of class \samp{Mods}, defining the mean vectors
#' under which the power should be calculated
#' @param n,sigma,S Either a vector \samp{n} and \samp{sigma} or \samp{S} need
#' to be specified.  When \samp{n} and \samp{sigma} are specified it is assumed
#' computations are made for a normal homoscedastic ANOVA model with group
#' sample sizes given by \samp{n} and residual standard deviation \samp{sigma},
#' i.e. the covariance matrix used for the estimates is thus
#' `sigma^2*diag(1/n)` and the degrees of freedom are calculated as
#' `sum(n)-nrow(contMat)`. When a single number is specified for \samp{n}
#' it is assumed this is the sample size per group and balanced allocations are
#' used.\cr
#' 
#' When \samp{S} is specified this will be used as covariance matrix for the
#' estimates.
#' @param placAdj Logical, if true, it is assumed that the standard deviation
#' or variance matrix of the placebo-adjusted estimates are specified in
#' \samp{sigma} or \samp{S}, respectively. The contrast matrix has to be
#' produced on placebo-adjusted scale, see [optContr()], so that the
#' coefficients are no longer contrasts (i.e. do not sum to 0).
#' @param alternative Character determining the alternative for the multiple
#' contrast trend test.
#' @param df Degrees of freedom to assume in case \samp{S} (a general
#' covariance matrix) is specified. When \samp{n} and \samp{sigma} are
#' specified the ones from the corresponding ANOVA model are calculated.
#' @param critV Critical value, if equal to \samp{TRUE} the critical value will
#' be calculated. Otherwise one can directly specify the critical value here.
#' @param control A list specifying additional control parameters for the
#' \samp{qmvt} and \samp{pmvt} calls in the code, see also
#' \samp{mvtnorm.control} for details.
#' @return Numeric containing the calculated power values
#' @author Bjoern Bornkamp
#' @seealso [powN()], [sampSizeMCT()],
#' [MCTtest()], [optContr()], [Mods()]
#' @references Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and
#' analysis of dose finding studies combining multiple comparisons and modeling
#' procedures, *Journal of Biopharmaceutical Statistics*, **16**,
#' 639--656
#' @examples
#' 
#' ## look at power under some dose-response alternatives
#' ## first the candidate models used for the contrasts
#' doses <- c(0,10,25,50,100,150)
#' ## define models to use as alternative 
#' fmodels <- Mods(linear = NULL, emax = 25,
#'                 logistic = c(50, 10.88111), exponential= 85,
#'                 betaMod=rbind(c(0.33,2.31),c(1.39,1.39)),
#'                 doses = doses, addArgs=list(scal = 200),
#'                 placEff = 0, maxEff = 0.4)
#' ## plot alternatives
#' plot(fmodels)
#' ## power for to detect a trend
#' contMat <- optContr(fmodels, w = 1)
#' powMCT(contMat, altModels = fmodels, n = 50, alpha = 0.05, sigma = 1)
#' 
#' \dontrun{
#' ## power under the Dunnett test
#' ## contrast matrix for Dunnett test with informative names
#' contMatD <- rbind(-1, diag(5))
#' rownames(contMatD) <- doses
#' colnames(contMatD) <- paste("D", doses[-1], sep="")
#' powMCT(contMatD, altModels = fmodels, n = 50, alpha = 0.05, sigma = 1)
#' 
#' ## now investigate power of the contrasts in contMat under "general" alternatives
#' altFmods <- Mods(linInt = rbind(c(0, 1, 1, 1, 1),
#'                                   c(0.5, 1, 1, 1, 0.5)),
#'                  doses=doses, placEff=0, maxEff=0.5)
#' plot(altFmods)
#' powMCT(contMat, altModels = altFmods, n = 50, alpha = 0.05, sigma = 1)
#' 
#' ## now the first example but assume information only on the
#' ## placebo-adjusted scale
#' ## for balanced allocations and 50 patients with sigma = 1 one obtains
#' ## the following covariance matrix
#' S <- 1^2/50*diag(6)
#' ## now calculate variance of placebo adjusted estimates
#' CC <- cbind(-1,diag(5))
#' V <- (CC)%*%S%*%t(CC)
#' linMat <- optContr(fmodels, doses = c(10,25,50,100,150),
#'                    S = V, placAdj = TRUE)
#' powMCT(linMat, altModels = fmodels, placAdj=TRUE,
#'        alpha = 0.05, S = V, df=6*50-6) # match df with the df above
#' }
#' @export
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
