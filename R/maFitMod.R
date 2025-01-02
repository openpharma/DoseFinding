#' This function fits dose-response models in a bootstrap model
#' averaging approach motivated from the bagging procedure (Breiman
#' XYZ). Given summary estimates for the outcome at doses, the
#' function samples summary data from the multivariate normal
#' distribution. For each sample dose-response models are fit to these
#' summary estimates and within each bootstrap the best model
#' according to the gAIC is selected.
#'
#' @title maFitMod - Fit dose-response models via bootstrap model
#'   averaging (bagging)
#' @aliases predict.maFit plot.maFit print.maFit
#' @param dose Numeric specifying the dose variable.
#' @param resp Numeric specifying the response estimate corresponding
#'   to the doses in \code{dose}
#' @param S Covariance matrix associated with the dose-response
#'   estimate specified via \code{resp}
#' @param models dose-response models to fit
#' @param nSim Number of bootstrap simulations
#' @param control Same as the control argument in
#'   \code{\link{fitMod}}.
#' @param bnds Bounds for non-linear parameters. This needs to be a
#'   list with list entries corresponding to the selected bounds. The
#'   names of the list entries need to correspond to the model
#'   names. The \code{\link{defBnds}} function provides the default
#'   selection.
#' @param addArgs List containing two entries named "scal" and "off"
#'   for the "betaMod" and "linlog" model. When addArgs is NULL the
#'   following defaults are used \samp{list(scal = 1.2*max(doses), off
#'   = 0.01*max(doses))}
#' @return An object of class \samp{maFit}, which contains the fitted
#'   dose-response models \samp{DRMod} objects, as well as which model
#'   was selected in each bootstrap and basic input parameters.
#' @author XYZ
#' @seealso \code{\link{fitMod}}, \code{\link{bFitMod}}, \code{\link{drmodels}}
#' @examples
#' data(biom)
#' ## produce first stage fit (using dose as factor)
#' anMod <- lm(resp~factor(dose)-1, data=biom)
#' drFit <- coef(anMod)
#' S <- vcov(anMod)
#' dose <- sort(unique(biom$dose))
#' ## fit an emax and sigEmax model (increase nSim for real use)
#' mFit <- maFitMod(dose, drFit, S, model = c("emax", "sigEmax"), nSim = 10)
#' mFit
#' plot(mFit, plotData = "meansCI")
#' ED(mFit, direction = "increasing", p = 0.9)
#' @export
maFitMod <- function(dose, resp, S, models, 
                     nSim = 1000,
                     control, bnds, addArgs = NULL){

  builtIn <- c("linlog", "linear", "quadratic", "linInt", "emax",
               "exponential", "logistic", "betaMod", "sigEmax")
  if(missing(models))
    stop("Need to specify the models that should be fitted")
  modelNum <- match(models, builtIn)
  if(any(is.na(modelNum))){
    stop_str <- sprintf("Invalid dose-response model specified: %s",
                        paste(models[is.na(modelNum)], collapse = ", "))
    stop(stop_str)
  }

  if(!missing(bnds)){
    if(!is.list(bnds))
      stop("bnds needs to be a list")
  }

  if(any(modelNum > 4)){ # non-linear model -> needs bounds
    if(missing(bnds)){
      message("Message: Need bounds in \"bnds\" for nonlinear models, using default bounds from \"defBnds\".")
      bnds <- defBnds(max(dose))
    } else{
      nonlin_models <- builtIn[modelNum[modelNum > 4]]
      bnds <- sapply(nonlin_models, function(x) if(x %in% names(bnds)){bnds[[x]]} else{
        message("Message: Need bounds in \"bnds\" for nonlinear models, using default bounds from \"defBnds\".");
        defBnds(max(dose))[[x]]}, simplify = FALSE)
    }
  }

  ## parametric bootstrap
  sims <- mvtnorm::rmvnorm(nSim, resp, S)
  fits <- vector("list", nSim)
  selModel <- character(nSim)
  for(i in 1:nSim){
    mod_fits <- lapply(models, function(mod){
      fitMod(dose, sims[i,], model = mod, S = S,
             type = "general", bnds = bnds[[mod]],
             addArgs = addArgs)
    })
    index <- which.min(sapply(mod_fits, gAIC))
    fits[[i]] <- mod_fits[[index]]
    selModel[i] <- models[index]
  }
  out <- list(fits = fits,
              selModels = selModel,
              args = list(dose = dose, resp = resp, S=S,
                          models = models))
  class(out) <- "maFit"
  out
}

#' @export
#' @param object Object of class maFit
#' @param summaryFct If equal to NULL predictions are calculated for
#'   each sampled parameter value. Otherwise a summary function is
#'   applied to the dose-response predictions for each parameter
#'   value.  The default is to calculate 0.025, 0.25, 0.5, 0.75, 0.975
#'   quantiles of the predictions for each dose.
#' @param doseSeq Where to calculate predictions. 
#' @param ... Further arguments (currently ignored)
predict.maFit <- function(object,
                          summaryFct = function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
                          doseSeq = NULL,
                          ...){
  if(is.null(doseSeq))
    stop("Need to provide doseSeq argument")
  nSim <- length(object$selModel)
  pred <- matrix(nrow = nSim, ncol = length(doseSeq))
  colnames(pred) <- doseSeq
  rownames(pred) <- 1:nSim
  for(i in 1:nSim){
    pred[i,] <- predict(object$fits[[i]], doseSeq = doseSeq, predType = "ls-means")
  }
  if(!is.null(summaryFct)){
    out0 <- apply(pred, 2, summaryFct)
    out <- matrix(out0, ncol = length(doseSeq))
  } else {
    out <- pred
  }
  colnames(out) <- doseSeq
  out
}

#' @export
print.maFit <- function(x, digits = 3, ...){
  cat("Bootstrap model averaging fits\n")
  
  cat("\nSpecified summary data:\n")
  dose_str <- paste(x$args$dose, collapse = ", ")
  cat(sprintf("doses: %s\n", dose_str))
  mn_str <- paste(round(x$args$resp, digits), collapse = ", ")
  cat(sprintf("mean: %s\n", mn_str))
  cat("Covariance Matrix:\n")
  S2 <- x$args$S
  rownames(S2) <- colnames(S2) <- x$args$dose
  print(round(S2, digits))

  cat(sprintf("\nModels fitted: %s\n", paste(x$args$models, collapse = ", ")))
  nSim <- length(x$selModels)
  cat(sprintf("\nModels selected by gAIC on bootstrap samples (nSim = %s)\n", nSim))
  tab0 <- table(x$selModels)
  out <- as.numeric(tab0)
  names(out) <- names(tab0)
  print(out)
}

#' @param x object of class maFit
#' @param plotData Determines how the original data are plotted:
#'   Either as means or as means with CI or not at all. The level of the CI
#'   is determined by the argument \samp{level}.
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param level Level for CI, when plotData is equal to
#'   \samp{meansCI}.
#' @param trafo Plot the fitted models on a transformed scale
#'   (e.g. probability scale if models have been fitted on log-odds
#'   scale). The default for \samp{trafo} is the identity function.
#' @param lenDose Number of grid values to use for display.
#' @param ... Additional parametes (unused)
#' @export
plot.maFit <- function(x, 
                       plotData = c("means", "meansCI", "none"),
                       xlab = "Dose", ylab = "Response",
                       title = NULL,
                       level = 0.95, trafo = function(x) x,
                       lenDose = 201, ...){
  if(!inherits(trafo, "function"))
    stop("trafo needs to be a function")
  plotData <- match.arg(plotData)
  dsq <- seq(0, max(x$args$dose), length = lenDose)
  preds <- predict(x, doseSeq = dsq, summaryFct = NULL)
  tail_prob <- (1-level)/2
  pdat <- data.frame(
    dose = dsq,
    median = trafo(apply(preds, 2, function(x) quantile(x, 0.5))),
    UB = trafo(apply(preds, 2, function(x) quantile(x, 1-tail_prob))),
    LB = trafo(apply(preds, 2, function(x) quantile(x, tail_prob)))
  )
  if (plotData %in% c("meansCI", "means")) {
    pmdat <- data.frame(dose = x$args$dose,
                        median = trafo(x$args$resp))
    sdev <- sqrt(diag(x$args$S))
    crit <- qnorm(1 - tail_prob)
    LBm <- UBm <- numeric(length(x$args$dose))
    for (i in 1:length(x$args$dose)) {
      LBm[i] <- trafo(x$args$resp[i] - crit * sdev[i])
      UBm[i] <- trafo(x$args$resp[i] + crit * sdev[i])
    }
    pmdat$LBm <- LBm
    pmdat$UBm <- UBm
  }
  pp <- ggplot2::ggplot(pdat, ggplot2::aes(x=dose, y=median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = LB, ymax = UB), alpha=0.2)+
    ggplot2::geom_line()+
    ggplot2::xlab(xlab)+
    ggplot2::ylab(ylab)+
    ggplot2::theme_bw()+
    ggplot2::scale_x_continuous(breaks=x$args$dose)+
    ggplot2::scale_y_continuous(breaks=pretty(c(pdat$UB, pdat$LB), 8))
  if(plotData %in% c("means", "meansCI")){
    pp <- pp +
      ggplot2::geom_point(ggplot2::aes(x=dose, y=median), data=pmdat)
    if(plotData == "meansCI")
      pp <- pp +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=LBm, ymax=UBm), data=pmdat, width = 0)
  }
  if(!is.null(title)){
    if(!is.character(title))
      stop("title needs to be a character")
    pp <- pp + ggplot2::ggtitle(title)
  }
  pp
}
