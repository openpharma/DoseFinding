

#' Fit a dose-response model using Bayesian or bootstrap methods.
#'
#' For \samp{type = "Bayes"}, MCMC sampling from the posterior distribution of the dose-response model is done. The
#' function assumes a multivariate normal distribution for `resp` with covariance matrix `S`, and this is
#' taken as likelihood function and combined with the prior distributions specified in prior to form the posterior
#' distribution.
#'
#' For \samp{type = "bootstrap"}, a multivariate normal distribution for `resp` with covariance matrix `S` is
#' assumed, and a large number of samples is drawn from this distribution. For each draw the fitMod function with
#' \samp{type = "general"} is used to fit the draws from the multivariate normal distribution.
#'
#' Componentwise univariate slice samplers are implemented (see Neal, 2003) to sample from the posterior distribution.
#'
#' @aliases bFitMod coef.bFitMod predict.bFitMod plot.bFitMod
#' @param dose Numeric specifying the dose variable.
#' @param resp Numeric specifying the response estimate corresponding to the doses in `dose`
#' @param S Covariance matrix associated with the dose-response estimate specified via `resp`
#' @param model Dose-response model to fit, possible models are "linlog", "linear", "quadratic", "emax", "exponential",
#'   "sigEmax", "betaMod" and "logistic", see [drmodels()].
#' @param placAdj Whether or not estimates in "placAdj" are placebo-adjusted (note that the linear in log and the
#'   logistic model cannot be fitted for placebo-adjusted data)
#' @param type Character with allowed values "Bayes" and "bootstrap", Determining whether samples are drawn from the
#'   posterior, or the bootstrap distribution.
#' @param start Optional starting values for the dose-response parameters in the MCMC algorithm.
#' @param prior List containing the information regarding the prior distributions for \samp{type = "Bayes"}.  The list
#'   needs to have as many entries as there are model parameters. The ordering of the list entries should be the same as
#'   in the arguments list of the model see (see [drmodels()]).  For example for the Emax model the first
#'   entry determines the prior for e0, the second to eMax and the third to ed50.
#'
#'   For each list entry the user has the choice to choose from 4 possible
#' distributions: \itemize{ \item `norm`: Vector of length 2 giving mean
#' and standard deviation.  \item `t`: Vector of length 3 giving median,
#' scale and degrees of freedom of the t-distribution.  \item `lnorm`:
#' Vector of length 2 giving mean and standard deviation on log scale.  \item
#' `beta`: Vector of length 4 giving lower and upper bound of the beta
#' prior as well as the alpha and beta parameters of the beta distribution }
#' @param nSim Desired number of samples to produce with the algorithm
#' @param MCMCcontrol List of control parameters for the MCMC algorithm
#' \itemize{ \item `thin` Thinning rate. Must be a positive integer.
#' \item `w` Numeric of same length as number of parameters in the model,
#' specifies the width parameters of the slice sampler.  \item `adapt`
#' Logical whether to adapt the `w` (width) parameter of the slice sampler
#' in a short trial run. The widths are chosen as IQR/1.3 of the trial run.  }
#' @param control Same as the control argument in [fitMod()].
#' @param bnds Bounds for non-linear parameters, in case \samp{type = "bootstrap"}. If missing the the default bounds
#'   from [defBnds()] is used.
#' @param addArgs List containing two entries named "scal" and "off" for the "betaMod" and "linlog" model. When addArgs
#'   is NULL the following defaults are used \samp{list(scal = 1.2*max(doses), off = 0.01*max(doses))}
#' @param x,object A bFitMod object
#' @param ...  Additional arguments are ignored.
#' @return An object of class bFitMod, which is a list containing the matrix of posterior simulations plus some
#'   additional information on the fitted model.
#' @author Bjoern Bornkamp
#' @seealso [fitMod()]
#' @references Neal, R. M. (2003), Slice sampling, Annals of Statistics, 31, 705-767
#' @examples
#' data(biom)
#  ## produce first stage fit (using dose as factor)
#' anMod <- lm(resp~factor(dose)-1, data=biom)
#' drFit <- coef(anMod)
#' S <- vcov(anMod)
#' dose <- sort(unique(biom$dose))
#' ## define prior list
#' ## normal prior for E0 (mean=0 and sdev=10)
#' ## normal prior for Emax (mean=0 and sdev=100)
#' ## beta prior for ED50: bounds: [0,1.5] parameters shape1=0.45, shape2=1.7
#' prior <- list(norm = c(0, 10), norm = c(0,100), beta=c(0,1.5,0.45,1.7))
#' ## now fit an emax model
#' gsample <- bFitMod(dose, drFit, S, model = "emax",
#'                    start = c(0, 1, 0.1), nSim = 1000, prior = prior)
#' ## summary information
#' gsample
#' ## samples are stored in
#' head(gsample$samples)
#' ## predict 0.025, 0.25, 0.5, 0.75, 0.975 Quantile at 0, 0.5 and 1
#' predict(gsample, doseSeq = c(0, 0.5, 1))
#' ## simple plot function
#' plot(gsample)
#'
#' ## now look at bootstrap distribution
#' gsample <- bFitMod(dose, drFit, S, model = "emax", type = "bootstrap",
#'                    nSim = 100, bnds = defBnds(1)$emax)
#' plot(gsample)
#'
#' ## now fit linear interpolation
#' prior <- list(norm = c(0,1000), norm = c(0,1000),
#'               norm = c(0,1000), norm = c(0,1000), norm = c(0,100))
#' gsample <- bFitMod(dose, drFit, S, model = "linInt",
#'                    start = rep(1,5), nSim = 1000, prior = prior)
#' gsample <- bFitMod(dose, drFit, S, model = "linInt", type = "bootstrap",
#'                    nSim = 100)
#'
#' ## data fitted on placebo adjusted scale
#' data(IBScovars)
#' anovaMod <- lm(resp~factor(dose)+gender, data=IBScovars)
#' drFit <- coef(anovaMod)[2:5] # placebo adjusted estimates at doses
#' vCov <- vcov(anovaMod)[2:5,2:5]
#' dose <- sort(unique(IBScovars$dose))[-1]
#' prior <- list(norm = c(0,100), beta=c(0,6,0.45,1.7))
#' ## Bayes fit
#' gsample <- bFitMod(dose, drFit, vCov, model = "emax", placAdj=TRUE,
#'                    start = c(1, 0.1), nSim = 1000, prior = prior)
#' ## bootstrap fit
#' gsample <- bFitMod(dose, drFit, vCov, model = "emax", placAdj=TRUE,
#'                    type = "bootstrap", start = c(1, 0.1),
#'                    nSim = 100, prior = prior, bnds = c(0.01,6))
#' ## calculate target dose estimate
#' TD(gsample, Delta = 0.2)
#' ## now fit linear interpolation
#' prior <- list(norm = c(0,1000), norm = c(0,1000), norm = c(0,1000), norm = c(0,100))
#' gsample <- bFitMod(dose, drFit, vCov, model = "linInt", placAdj=TRUE,
#'                    start = rep(1,4), nSim = 1000, prior = prior)
#' gsample <- bFitMod(dose, drFit, vCov, model = "linInt", type = "bootstrap",
#'                    placAdj = TRUE, nSim = 100)
#' @export                    
bFitMod <- function(dose, resp, model, S, placAdj = FALSE,
                    type = c("Bayes", "bootstrap"),
                    start = NULL, prior = NULL, nSim = 1000,
                    MCMCcontrol = list(), control = NULL, bnds, 
                    addArgs = NULL){
  if(placAdj & model %in% c("linlog", "logistic"))
    stop("logistic and linlog models can only be fitted with placAdj")
  nD <- length(dose)
  if (length(resp) != nD) 
    stop("dose and resp need to be of the same size")
  dose <- as.numeric(dose)
  if (any(dose < -.Machine$double.eps)) 
    stop("dose values need to be non-negative")
  if (!is.numeric(dose)) 
    stop("dose variable needs to be numeric")
  resp <- as.numeric(resp)
  ## order dose and resp increasingly
  ord <- order(dose)
  dose <- dose[ord]
  resp <- resp[ord]
  if (nrow(S) != nD | ncol(S) != nD) 
    stop("S and dose have non-conforming size")
  if (missing(model)) 
    stop("need to specify the model that should be fitted")
  scal <- off <- nodes <- NULL
  if(model %in% c("linlog", "betaMod")){
    lst <- getAddArgs(addArgs, dose)
    if(model == "betaMod")
      scal <- lst$scal
    if(model == "linlog")
      off <- lst$off
  }
  if(model == "linInt")
    nodes <- dose
  
  ## model number
  builtIn <- c("linear", "linlog", "quadratic", "linInt", "emax",
               "logistic", "exponential", "sigEmax", "betaMod")
  modNr <- match(model, builtIn)
  if(is.na(modNr))
    stop("invalid model selected")
  ## number of parameters
  nPar <- as.integer(c(2, 2, 3, length(dose), 3, 4, 3, 4, 4)[modNr])

  type <- match.arg(type)
  if(type == "Bayes"){
    res <- bFitMod.Bayes(dose, resp, S, model, placAdj,
                         start, prior, nSim, MCMCcontrol,
                         off, scal, nPar, modNr)
    res <- matrix(res, nrow = nSim, ncol = nPar)
    if(placAdj & model != "linInt")
      res <- res[,-1, drop = FALSE]
  } else { ## bootstrap
    res <- bFitMod.bootstrap(dose, resp, S, model, placAdj,
                             nSim, control, bnds, off, scal,
                             nodes)
  }
  
  out <- list(samples = res)
  if(model != "linInt"){
    nams <- names(formals(model))[-1]
  } else {
    nams <- paste("d", dose, sep="")
  }
  if(modNr %in% c(2,9))
    nams <- nams[-length(nams)]
  if(placAdj & model != "linInt")
    nams <- nams[-1]
  colnames(out$samples) <- nams
  attr(out, "model") <- model
  lst <- list(dose, resp, S)
  doseNam <- as.list(match.call())$dose
  respNam <- as.list(match.call())$resp
  attr(out, "doseRespNam") <- as.character(c(doseNam, respNam))
  names(lst) <- c(doseNam, respNam, "S")
  attr(out, "data") <- lst
  attr(out, "type") <- type
  attr(out, "call") <- match.call()
  attr(out, "placAdj") <- placAdj
  attr(out, "prior") <- prior
  attr(out, "scal") <- scal
  attr(out, "off") <- off
  attr(out, "nodes") <- nodes
  class(out) <- "bFitMod"
  out
}

#' @export
print.bFitMod <- function(x, digits = 3, ...){
  ## print brief summary of MCMC samples
  doseNam <- attr(x, "doseRespNam")[1]
  respNam <- attr(x, "doseRespNam")[2]
  resp <- attr(x, "data")[[respNam]]
  names(resp) <- attr(x, "data")[[doseNam]]
  cat("Dose Response Model\n\n")
  cat(paste("Model:", attr(x, "model")), "\n\n")
  if(attr(x, "type") == "Bayes"){
    cat("Summary of posterior draws\n")
    func <- function(x){
      c(mean=mean(x), sdev=sd(x),
        quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)),
        n.eff=ess.mcmc(x))
    }
    print(t(apply(x$samples, 2, func)), digits=digits)
  } else {
    cat("Summary of bootstrap draws\n")
    func <- function(x){
      c(mean=mean(x), sdev=sd(x),
        quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)))
    }
    print(t(apply(x$samples, 2, func)), digits=digits)
  }
  cat("\nFitted to:\n")
  print(signif(resp, digits+2))
}

#' Make predictions from fitted dose-response model
#'
#' @param predType,summaryFct,doseSeq,lenSeq Arguments for the predict method.
#'
#'   \samp{predType}: predType determines whether predictions are returned for the dose-response curve or the effect
#'   curve (difference to placebo).
#'
#'   \samp{summaryFct}: If equal to NULL predictions are calculated for each sampled parameter value. Otherwise a
#'   summary function is applied to the dose-response predictions for each parameter value.  The default is to calculate
#'   0.025, 0.25, 0.5, 0.75, 0.975 quantiles of the predictions for each dose.
#'
#'   \samp{doseSeq}: Where to calculate predictions. If not specified predictions are calculated on a grid of length
#'   \samp{lenSeq} between minimum and maximum dose.
#'
#'   \samp{lenSeq}: Length of the default grid where to calculate predictions.
#'
#' @rdname bFitMod
#' @method predict bFitMod
#' @export
predict.bFitMod <- function(object, predType = c("full-model", "effect-curve"),
                            summaryFct = function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                            doseSeq = NULL, lenSeq = 101, ...){
  predType <- match.arg(predType)
  doseNam <- attr(object, "doseRespNam")[1]
  if (is.null(doseSeq)) {
    doseSeq <- seq(0, max(attr(object, "data")[[doseNam]]), length = lenSeq)
  }
  model <- attr(object, "model")
  scal <- attr(object, "scal")
  off <- attr(object, "off")
  placAdj <- attr(object, "placAdj")
  if(placAdj){
    nodes <- c(0,attr(object, "data")[[doseNam]])
  } else {
    nodes <- attr(object, "data")[[doseNam]]
  }
  out <- predSamples(samples = object$samples, doseSeq = doseSeq,
                     placAdj = placAdj, model = model, scal = scal,
                     off = off, nodes = nodes)

  if(predType == "effect-curve"){
    out <- out - out[,1]
  }
  if(!is.null(summaryFct)){
    out0 <- apply(out, 2, summaryFct)
    out <- matrix(out0, ncol = ncol(out))
  }
  colnames(out) <- doseSeq
  out
}


#' Plot fitted dose-response model
#'
#' @param plotType,quant,plotData,level,lenDose Arguments for plot method.
#'
#'   \samp{plotType}: Determining whether the dose-response curve or the effect curve should be plotted.
#'
#'   \samp{quant}: Vector of quantiles to display in plot
#'
#'   \samp{plotData}: Determines how the original data are plotted: Either as means or as means with CI or not. The
#'   level of the CI is determined by the argument \samp{level}.
#'
#'   \samp{level}: Level for CI, when plotData is equal to \samp{meansCI}.
#'
#'   \samp{lenDose}: Number of grid values to use for display.
#'
#' @rdname bFitMod
#' @method plot bFitMod
#' @export
plot.bFitMod <- function (x, plotType = c("dr-curve", "effect-curve"),
                          quant = c(0.025, 0.5, 0.975), 
                          plotData = c("means", "meansCI", "none"),
                          level = 0.95, lenDose = 201, ...){
  addArgs <- list(...)
  plotType <- match.arg(plotType)
  doseNam <- attr(x, "doseRespNam")[1]
  respNam <- attr(x, "doseRespNam")[2]
  dose <- attr(x, "data")[[doseNam]]
  resp <- attr(x, "data")[[respNam]]
  doseSeq <- seq(0, max(dose), length = lenDose)
  plotData <- match.arg(plotData)
  placAdj <- attr(x, "placAdj")
  sumFct <- function(x){
    quantile(x, probs = quant)
  }
  if (placAdj) 
    plotType <- "effect-curve"
  if (plotType == "effect-curve") {
    pred <- predict(x, predType = plotType, summaryFct = sumFct,
                    doseSeq = doseSeq)
    main <- "Effect Curve"
    if (placAdj) {
      if (plotData == "meansCI") {
        sdev <- sqrt(diag(attr(x, "data")$S))
        q <- qnorm(1 - (1 - level)/2)
        LBm <- UBm <- numeric(length(dose))
        for (i in 1:length(dose)) {
          LBm[i] <- resp[i] - q * sdev[i]
          UBm[i] <- resp[i] + q * sdev[i]
        }
      }
      else {
        LBm <- UBm <- NULL
      }
    }
    else {
      LBm <- UBm <- NULL
    }
  }
  if (plotType == "dr-curve") {
    pred <- predict(x, predType = "full-model", summaryFct = sumFct,
                    doseSeq = doseSeq)
    main <- "Dose-Response Curve\n"
    if (plotData == "meansCI") {
      sdev <- sqrt(diag(attr(x, "data")$S))
      q <- qnorm(1 - (1 - level)/2)
      LBm <- UBm <- numeric(length(dose))
      for (i in 1:length(dose)) {
        LBm[i] <- resp[i] - q * sdev[i]
        UBm[i] <- resp[i] + q * sdev[i]
      }
    }
    else {
      LBm <- UBm <- NULL
    }
  }
  rng <- range(c(pred, resp, LBm, UBm))
  dff <- diff(rng)
  ylim <- c(rng[1] - 0.02 * dff, rng[2] + 0.02 * dff)
  callList <- list(doseSeq, t(pred), type = "l", xlab = doseNam, 
                   ylim = ylim, ylab = respNam, main = main,
                   lty=1, col=1)
  callList[names(addArgs)] <- addArgs
  do.call("matplot", callList)
  
  if (plotType == "dr-curve" | placAdj) {
    if (plotData == "means") 
      points(dose, resp, pch = 19, cex = 0.75)
    if (plotData == "meansCI") {
      points(dose, resp, pch = 19, cex = 0.75)
      for (i in 1:length(dose)) {
        lines(c(dose[i], dose[i]), c(LBm[i], UBm[i]), 
              lty = 2)
      }
    }
  }
  res <- list(doseSeq = doseSeq)
  attr(res, "level") <- level
  attr(res, "ylim") <- ylim
  res$mean <- pred
  invisible(res)
}

#' Extract dose-response model coefficients
#'
#' @param x,object A bFitMod object
#' @param ... Additional arguments are ignored.
#'
#' @rdname bFitMod
#' @method coef bFitMod
#' @export
coef.bFitMod <- function (object, ...){
  object$samples
}
