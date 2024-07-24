## function for sample size calculation and functions to evaluate
## performance metrics for different sample sizes



#' Sample size calculations
#'
#'
#' The \samp{sampSize} function implements a bisection search algorithm for sample size calculation. The user can hand
#' over a general target function (via \samp{targFunc}) that is then iterated so that a certain \samp{target} is
#' achieved. The \samp{sampSizeMCT} is a convenience wrapper of \samp{sampSize} for multiple contrast tests using the
#' power as target function.
#'
#' The \samp{targN} functions calculates a general target function for different given sample sizes. The \samp{powN}
#' function is a convenience wrapper of \samp{targN} for multiple contrast tests using the power as target function.
#'
#'
#' @aliases sampSize sampSizeMCT targN plot.targN powN
#' @param upperN,lowerN Upper and lower bound for the target sample size. \code{lowerN} defaults to
#'   \code{floor(upperN/2)}.
#' @param targFunc,target The target function needs to take as an input the vector of sample sizes in the different dose
#'   groups. For \samp{sampSize} it needs to return a univariate number. For function \samp{targN} it should return a
#'   numerical vector.\cr \cr Example: \samp{targFunc} could be a function that calculates the power of a test, and
#'   \samp{target} the desired target power value.  \cr For function \samp{sampSize} the bisection search iterates the
#'   sample size so that a specific target value is achieved (the implicit assumption is that targFunc is monotonically
#'   increasing in the sample size).\cr \cr Function \samp{targN} simply calculates \samp{targFunc} for a given set of
#'   sample sizes.
#' @param tol A positive numeric value specifying the tolerance level for the bisection search algorithm. Bisection is
#'   stopped if the \samp{targFunc} value is within \samp{tol} of \samp{target}.
#' @param alRatio Vector describing the relative patient allocations to the dose groups up to proportionality, e.g.
#'   \samp{rep(1, length(doses))} corresponds to balanced allocations.
#' @param Ntype One of "arm" or "total". Determines, whether the sample size in the smallest arm or the total sample
#'   size is iterated in bisection search algorithm.
#' @param verbose Logical value indicating if a trace of the iteration progress of the bisection search algorithm should
#'   be displayed.
#' @param ...  Arguments directly passed to the \code{\link{powMCT}} function in the \samp{sampSizeMCT} and \samp{powN}
#'   function.
#'
#'   The \samp{placAdj} argument needs to be \samp{FALSE} (which is the default value for this argument). If sample size
#'   calculations are desired for a placebo-adjusted formulation use \samp{sampSize} or \samp{targN} directly.
#'
#'   In case \code{S} is specified, the specified matrix needs to be proportional to the (hypothetical) covariance
#'   matrix of one single observation. The covariance matrix used for sample size calculation is 1/N*S, where N is the
#'   total sample size. Hence \samp{Ntype == "total"} needs to be used if
#' \code{S} is specified. When \code{S} is specified, automatically \samp{df =
#' Inf} is assumed in the underlying \samp{powMCT} calls.
#'
#'   For a homoscedastic normally distributed response variable only \samp{sigma} needs to be specified, as the sample
#'   size \samp{n} is iterated in the different \samp{powMCT} calls.
#'
#' @author Jose Pinheiro, Bjoern Bornkamp
#' @seealso \code{\link{powMCT}}
#' @references Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and analysis of dose finding studies
#'   combining multiple comparisons and modeling procedures, \emph{Journal of Biopharmaceutical Statistics}, \bold{16},
#'   639--656
#'
#'   Pinheiro, J.C., Bornkamp, B. (2017) Designing Phase II Dose-Finding Studies: Sample Size, Doses and Dose Allocation
#'   Weights, in O'Quigley, J., Iasonos, A. and Bornkamp, B. (eds) Handbook of methods for designing, monitoring, and
#'   analyzing dose-finding trials, CRC press
#' @examples
#'
#' ## sampSize examples
#'
#' ## first define the target function
#' ## first calculate the power to detect all of the models in the candidate set
#' fmodels <- Mods(linear = NULL, emax = c(25),
#'                 logistic = c(50, 10.88111), exponential=c(85),
#'                 betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=TRUE, nrow=2),
#'                 doses = c(0,10,25,50,100,150), placEff=0, maxEff=0.4,
#'                 addArgs = list(scal=200))
#' ## contrast matrix to use
#' contMat <- optContr(fmodels, w=1)
#' ## this function calculates the power under each model and then returns
#' ## the average power under all models
#' tFunc <- function(n){
#'   powVals <- powMCT(contMat, altModels=fmodels, n=n, sigma = 1,
#'                     alpha=0.05)
#'   mean(powVals)
#' }
#'
#' ## assume we want to achieve 80% average power over the selected shapes
#' ## and want to use a balanced allocations
#' \dontrun{
#' sSize <- sampSize(upperN = 80, targFunc = tFunc, target=0.8,
#'                   alRatio = rep(1,6), verbose = TRUE)
#' sSize
#'
#'
#' ## Now the same using the convenience sampSizeMCT function
#' sampSizeMCT(upperN=80, contMat = contMat, sigma = 1, altModels=fmodels,
#'             power = 0.8, alRatio = rep(1, 6), alpha = 0.05)
#' ## Alternatively one can also specify an S matrix
#' ## covariance matrix in one observation (6 total observation result in a
#' ## variance of 1 in each group)
#' S <- 6*diag(6)
#' ## this uses df = Inf, hence a slightly smaller sample size results
#' sampSizeMCT(upperN=500, contMat = contMat, S=S, altModels=fmodels,
#'             power = 0.8, alRatio = rep(1, 6), alpha = 0.05, Ntype = "total")
#'
#'
#' ## targN examples
#' ## first calculate the power to detect all of the models in the candidate set
#' fmodels <- Mods(linear = NULL, emax = c(25),
#'                 logistic = c(50, 10.88111), exponential=c(85),
#'                 betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=TRUE, nrow=2),
#'                 doses = c(0,10,25,50,100,150), placEff=0, maxEff=0.4,
#'                 addArgs = list(scal=200))
#' ## corresponding contrast matrix
#' contMat <- optContr(fmodels, w=1)
#' ## define target function
#' tFunc <- function(n){
#'   powMCT(contMat, altModels=fmodels, n=n, sigma = 1, alpha=0.05)
#' }
#' powVsN <- targN(upperN = 100, lowerN = 10, step = 10, tFunc,
#'                 alRatio = rep(1, 6))
#' plot(powVsN)
#'
#' ## the same can be achieved using the convenience powN function
#' ## without the need to specify a target function
#' powN(upperN = 100, lowerN=10, step = 10, contMat = contMat,
#'      sigma = 1, altModels = fmodels, alpha = 0.05, alRatio = rep(1, 6))
#' }
#' @export
sampSize <- function (upperN, lowerN = floor(upperN/2),
                      targFunc, target, tol = 0.001, alRatio,
                      Ntype = c("arm", "total"), verbose = FALSE){
  ## target function to iterate
  func <- function(n){
    targFunc(n) - target
  }

  Ntype <- match.arg(Ntype)
  if (!missing(alRatio)) {
    if (any(alRatio <= 0)) {
      stop("all entries of alRatio need to be positive")
    } else {
      alRatio <- alRatio/sum(alRatio)
    }
    if(Ntype == "arm") {
      alRatio <- alRatio/min(alRatio)
    } 
  } else { ## by default assume
    stop("allocation ratios need to be specified")
  }
  
  ## first call
  upper <- func(round(upperN*alRatio))
  if(length(upper) > 1)
    stop("targFunc(n) needs to evaluate to a vector of length 1.")
  if(!is.numeric(upper))
    stop("targFunc(n) needs to evaluate to a numeric.")

  ## bracket solution
  if (upper < 0)
    message("upper limit for sample size is raised")

  while (upper < 0) {
    upperN <- 2 * upperN
    upper <- func(round(upperN*alRatio))
  }
  
  lower <- func(round(lowerN*alRatio))
  
  if (lower > 0) 
    message("lower limit for sample size is decreased")

  while (lower > 0) {
    lowerN <- round(lowerN/2)
    if (lowerN == 0) 
      stop("cannot find lower limit on n")
    lower <- func(round(lowerN*alRatio))
  }

  ## now start bisection
  if (verbose) {
    cat("Upper N:", upperN, "Upper value", round(upper+target, 4), "\n")
    cat("Lower N:", lowerN, "Lower value", round(lower+target, 4), "\n\n")
  }
  
  current <- tol+1
  niter <- 0
  ## bisect sample size until tolerance is achieved
  while (abs(current) > tol & (upperN > lowerN + 1)) {
    currN <- round((upperN + lowerN)/2)
    current <- func(round(currN * alRatio))
    if (current > 0) {
      upperN <- currN
    } else {
      lowerN <- currN
    }
    niter <- niter + 1
    if (verbose) {
      cat("Iter: ", niter, ", N = ", currN, ", current value = ",
          round(current+target, 4), "\n", sep = "")
    }
  }
  ## increase sample size so that the obtained value is larger than the target
  while (current < 0) {
    currN <- currN + 1
    current <- func(round(currN * alRatio))
  }

  res <- list(samp.size = round(currN * alRatio),
              target = round(current+target, 4))
  attr(res, "alRatio") <- round(alRatio/min(alRatio), 4)
  attr(res, "target") <- target
  attr(res, "Ntype") <- Ntype
  class(res) <- "sampSize"
  res
}

#' @export
print.sampSize <- function(x, ...){
  cat("Sample size calculation\n\n")
  cat("alRatio:", attr(x, "alRatio"), "\n")
  cat("Total sample size:", sum(x$samp.size), "\n")
  cat("Sample size per arm:", x$samp.size, "\n")
  cat("targFunc:", x$target,"\n")
}

#' Sample size calculations for multiple contrast tests
#'
#' @inheritParams sampSize
#' @param ...  Arguments directly passed to the \code{\link{powMCT}} function in the \samp{sampSizeMCT} and \samp{powN}
#'   function.
#' @param power,sumFct power is a numeric defining the desired summary power to achieve (in \samp{sampSizeMCT}). sumFct
#'   needs to be a function that combines the power values under the different alternatives into one value (in
#'   \samp{sampSizeMCT}).
#' @rdname sampSize
#' @export
sampSizeMCT <- function(upperN, lowerN = floor(upperN/2),
                        ...,
                        power, sumFct = mean,
                        tol = 0.001, alRatio, Ntype = c("arm", "total"), verbose = FALSE){
  ## function to calculate sample size for multiple contrast test
  ## if S is specified this needs to be the (hypothetical) covariance matrix
  ## for a total sample size of 1 patient
  Ntype <- match.arg(Ntype)
  args <- list(...)
  namargs <- names(args)
  if(is.element("placAdj", namargs)){
    if(args$placAdj)
      stop("placAdj needs to be FALSE for sampSizeMCT.
  Use sampSize directly in placebo-adjusted case.")
  }
  if(is.element("S", namargs)){
    S <- args[["S"]]
    if(Ntype == "arm"){
      Ntype <- "total"
      message("Only Ntype == \"total\" possible if S is specified")
    }
    if(is.element("df", namargs)){
      if(is.finite(args$df))
        message("df argument set to Inf, if S is specified.
Use sampSize directly in case exact df are required.")
    }
    args$df <- Inf
    tFunc <- function(n){
      N <- sum(n)
      Sn <- 1/N*S
      args$S <- Sn
      powVals <- do.call("powMCT", args)
      sumFct(powVals)
    }
  } else {
    if(is.element("n", namargs))
      stop("n is not allowed to be specified for sample size calculation")
    if(!is.element("sigma", namargs))
      stop("need sigma if S is not specified")
    tFunc <- function(n){
      powVals <- powMCT(n=n, ...)
      sumFct(powVals)
    }
  }
  sampSize(upperN, lowerN, targFunc = tFunc, target = power,
           alRatio = alRatio, Ntype = Ntype, verbose = verbose)
}

#' Calculate target function for given sample size
#'
#' @inheritParams sampSize
#' @param step Only needed for functions \samp{targN} and \samp{powN}. Stepsize for the sample size at which the target
#'   function is calculated. The steps are calculated via \code{seq(lowerN,upperN,by=step)}.
#' @param power,sumFct power is a numeric defining the desired summary power to achieve (in \samp{sampSizeMCT}).
#' @rdname sampSize
#' @export
targN <- function(upperN, lowerN, step, targFunc,
                  alRatio, Ntype = c("arm", "total"), sumFct = c("min", "mean", "max")){

  if(!is.character(sumFct))
    stop("sumFct needs to be a character vector")
  Ntype <- match.arg(Ntype)
  if (!missing(alRatio)) {
    if (any(alRatio <= 0)) {
      stop("all entries of alRatio need to be positive")
    } else {
      alRatio <- alRatio/sum(alRatio)
    }
    if(Ntype == "arm") {
      alRatio <- alRatio/min(alRatio)
    } 
  } else { ## by default assume 
    stop("allocation ratios need to be specified")
  }
  
  nseq <- seq(lowerN, upperN, by=step)
  out <-t(sapply(nseq, function(x){
    targFunc(round(x * alRatio))
  }))
  if(nrow(out) == 1 & length(nseq) > 1){
    out <- t(out)
    colnames(out) <- ""
  }
  out2 <- out
  for(i in 1:length(sumFct)){
    out2 <- cbind(out2, apply(out, 1, sumFct[i]))
  }
  dimnames(out2) <- list(nseq, c(colnames(out), sumFct))
  attr(out2, "alRatio") <- alRatio
  attr(out2, "sumFct") <- sumFct
  attr(out2, "Ntype") <- Ntype
  class(out2) <- "targN"
  out2
}

#' Calculate power for given sample size
#'
#' @inheritParams targN
#'
#' @rdname sampSize
#' @export
powN <- function(upperN, lowerN, step,
                 ...,
                 alRatio, Ntype = c("arm", "total"), sumFct = c("min", "mean", "max")){
  args <- list(...)
  namargs <- names(args)
  if(is.element("placAdj", namargs)){
    if(args$placAdj)
      stop("placAdj needs to be FALSE for powN.
  Use targN directly in placebo-adjusted case.")
  }
  Ntype <- match.arg(Ntype)
  if(is.element("S", namargs)){
    S <- args[["S"]]
    if(Ntype == "arm"){
      Ntype <- "total"
      message("Only Ntype == \"total\" possible if S is specified")
    }
    if(is.element("df", namargs)){
      if(is.finite(args$df))
        message("df argument set to Inf, if S is specified.
Use sampSize directly in case exact df are required.")
    }
    args$df <- Inf
    tFunc <- function(n){
      N <- sum(n)
      Sn <- 1/N*S
      args$S <- Sn
      do.call("powMCT", args)
    }
  } else {
    if(is.element("n", namargs))
      stop("n is not allowed to be specified for sample size calculation")
    if(!is.element("sigma", namargs))
      stop("need sigma if S is not specified")
    tFunc <- function(n)
      powMCT(n=n, ...)
  }
  targN(upperN=upperN, lowerN=lowerN, step=step, targFunc=tFunc,
        alRatio=alRatio, Ntype = Ntype, sumFct = sumFct)
}

#' Produce Trellis plot of targN object
#'
#' @param x,superpose,line.at,xlab,ylab arguments for the plot method of \samp{targN} and \samp{powN}, additional
#'   arguments are passed down to the low-level lattice plotting routines.
#'
#' @rdname sampSize
#' @method plot targN
#' @export
plot.targN <- function(x, superpose = TRUE, line.at = NULL, 
                       xlab = NULL, ylab = NULL, ...){
  nSeq <- as.integer(dimnames(x)[[1]])
  alRatio <- attr(x, "alRatio")
  unbN <- (length(unique(alRatio)) > 1)
  if (is.null(xlab)) {
    if(attr(x, "Ntype") == "total" | unbN){
      xlab <- "Overall sample size"
      nSeq <- sapply(nSeq, function(x){
        sum(round(x*alRatio))
      })
    } else {
      xlab <- "Sample size per dose (balanced)"
    }
  }
  nams <- dimnames(x)[[2]]
  ## separating model data from summary data
  x <- as.data.frame(unclass(x))
  nams <- names(x)
  nC <- ncol(x)
  pMatTr <- data.frame(targ = as.vector(unlist(x)), n = rep(nSeq, nC),
                       type = factor(rep(nams, each = length(nSeq)), levels = nams))
  if(superpose){
    panelFunc1 <- function(x, y, subscripts, groups, lineAt, ...) {
      lattice::panel.grid(h = -1, v = -1, col = "lightgrey", lty = 2)
      if(!is.null(line.at))
        lattice::panel.abline(h = lineAt, lty = 3, ..., col = "red")
      lattice::panel.superpose(x, y, subscripts, groups, ...)
    }
    trLn <- lattice::trellis.par.get("superpose.line")[c("col", "lwd", "lty")]
    for(i in seq(along = trLn)) {
      if(length(trLn[[i]]) > nC) trLn[[i]] <- trLn[[i]][1:nC]
    }
    ltplot <- lattice::xyplot(targ ~ n, pMatTr, groups = pMatTr$type, subscripts = TRUE,
                              panel = panelFunc1, type = "l", lineAt = line.at,
                              xlab = xlab, ylab = ylab,
                              key = list(lines = trLn, text = list(lab = nams), transparent = TRUE, 
                                         columns = ifelse(nC < 5, nC, min(4,ceiling(nC/min(ceiling(nC/4),3))))), ...)
  } else {                              # models in different panels
    panelFunc2 <- function(x, y, lineAt, ...) {
      lattice::panel.grid(h = -1, v = -1, col = "lightgrey", lty = 2)
      if(!is.null(line.at))
        lattice::panel.abline(h = lineAt, lty = 3, ..., col = "red") ## used 2 for consistency with above
      lattice::panel.xyplot(x, y, ...)
    }
    ltplot <- lattice::xyplot(targ ~ n | type, pMatTr, panel = panelFunc2,
                              type = "l", lineAt = line.at,
                              xlab = xlab, ylab = ylab, 
                              strip = function(...) lattice::strip.default(..., style = 1), ...)
  }
  print(ltplot)
}
