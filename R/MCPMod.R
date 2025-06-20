## wrapper function for MCTtest and fitMod calls


#' MCPMod - Multiple Comparisons and Modeling
#'
#' Tests for a dose-response effect using a model-based multiple contrast test (see [MCTtest()]), selects one
#' (or several) model(s) from the significant shapes, fits them using [fitMod()].  For details on the method
#' see Bretz et al. (2005).
#'
#'
#' @aliases MCPMod predict.MCPMod plot.MCPMod
#' @inheritParams MCTtest
#' @param selModel Optional character vector specifying the model selection
#'   criterion for dose estimation.  Possible values are \itemize{ \item
#'   `AIC`: Selects model with smallest AIC (this is the default) \item
#'   `maxT`: Selects the model corresponding to the largest t-statistic.
#'   \item `aveAIC`: Uses a weighted average of the models corresponding to
#'   the significant contrasts.  The model weights are chosen by the formula:
#'   \eqn{w_i = \exp(-0.5AIC_i)/\sum_i(\exp(-0.5AIC_i))}{w_i =
#'   exp(-0.5AIC_i)/sum(exp(-0.5AIC_i))} See Buckland et al. (1997) for details.
#'   } For \samp{type = "general"} the "gAIC" is used.
#' @param df Specify the degrees of freedom to use in case \samp{type = "general"}, for the call to
#'   [MCTtest()] and [fitMod()]. Infinite degrees of (\samp{df=Inf}) correspond to the multivariate
#'   normal distribution.  For type = "normal" the degrees of freedom deduced from the AN(C)OVA fit are used and this
#'   argument is ignored.
#' @param doseType,Delta,p \samp{doseType} determines the dose to estimate, ED or TD (see also [Mods()]), and
#'   \samp{Delta} and \samp{p} need to be specified depending on whether TD or ED is to be estimated.  See
#'   [TD()] and [ED()] for details.
#' @param bnds Bounds for non-linear parameters. This needs to be a list with list entries corresponding to the selected
#'   bounds. The names of the list entries need to correspond to the model names. The [defBnds()] function
#'   provides the default selection.
#' @param control Control list for the optimization.\cr A list with entries: "nlminbcontrol", "optimizetol" and
#'   "gridSize".
#'
#'   The entry nlminbcontrol needs to be a list and is passed directly to control argument in the nlminb function, that
#'   is used internally for models with 2 nonlinear parameters (e.g. sigmoid Emax or beta model).
#'
#'   The entry optimizetol is passed directly to the tol argument of the optimize function, which is used for models
#'   with 1 nonlinear parameters (e.g. Emax or exponential model).
#'
#'   The entry gridSize needs to be a list with entries dim1 and dim2 giving the size of the grid for the gridsearch in
#'   1d or 2d models.
#' @return An object of class \samp{MCPMod}, which contains the fitted \samp{MCTtest} object as well as the \samp{DRMod}
#'   objects and additional information (model selection criteria, dose estimates, selected models).
#' @author Bjoern Bornkamp
#' @seealso [MCTtest()], [fitMod()], [drmodels()]
#' @references Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining multiple comparisons and modeling
#'   techniques in dose-response studies, *Biometrics*, **61**, 738--748
#'
#'   Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and analysis of dose finding studies combining multiple
#'   comparisons and modeling procedures, *Journal of Biopharmaceutical Statistics*, **16**, 639--656
#'
#'   Pinheiro, J. C., Bretz, F., and Branson, M. (2006). Analysis of dose-response studies - modeling approaches,
#'   *in* N. Ting (ed.). *Dose Finding in Drug Development*, Springer, New York, pp. 146--171
#'
#'   Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014) Model-based dose finding under model uncertainty
#'   using general parametric models, *Statistics in Medicine*, **33**, 1646--1661
#'
#'   Schorning, K., Bornkamp, B., Bretz, F., & Dette, H. (2016). Model selection
#' versus model averaging in dose finding studies. *Statistics in
#' Medicine*, **35**, 4021--4040
#'
#'   Xun, X. and Bretz, F. (2017) The MCP-Mod methodology: Practical Considerations and The DoseFinding R package, in
#'   O'Quigley, J., Iasonos, A. and Bornkamp, B. (eds) Handbook of methods for designing, monitoring, and analyzing
#'   dose-finding trials, CRC press
#'
#'   Buckland, S. T., Burnham, K. P. and Augustin, N. H. (1997). Model selection an integral part of inference,
#'   *Biometrics*, **53**, 603--618
#'
#'   Seber, G.A.F. and Wild, C.J. (2003). Nonlinear Regression, Wiley.
#' @examples
#'
#' data(biom)
#' ## first define candidate model set (only need "standardized" models)
#' models <- Mods(linear = NULL, emax=c(0.05,0.2), linInt=c(1, 1, 1, 1),
#'                doses=c(0,0.05,0.2,0.6,1))
#' plot(models)
#' ## perform MCPMod procedure
#' MM <- MCPMod(dose, resp, biom, models, Delta=0.5)
#' ## a number of things can be done with an MCPMod object
#' MM # print method provides basic information
#' summary(MM) # more information
#' ## predict all significant dose-response models
#' predict(MM, se.fit=TRUE, doseSeq=c(0,0.2,0.4, 0.9, 1),
#'         predType="ls-means")
#' ## display all model functions
#' plot(MM, plotData="meansCI", CI=TRUE)
#'
#' ## now perform model-averaging
#' MM2 <- MCPMod(dose, resp, biom, models, Delta=0.5, selModel = "aveAIC")
#' sq <- seq(0,1,length=11)
#' pred <- predict(MM, doseSeq=sq, predType="ls-means")
#' modWeights <- MM2$selMod
#' ## model averaged predictions
#' pred <- do.call("cbind", pred)%*%modWeights
#' ## model averaged dose-estimate
#' TDEst <- MM2$doseEst%*%modWeights
#'
#' ## now an example using a general fit and fitting based on placebo
#' ## adjusted first-stage estimates
#' data(IBScovars)
#' ## ANCOVA fit model including covariates
#' anovaMod <- lm(resp~factor(dose)+gender, data=IBScovars)
#' drFit <- coef(anovaMod)[2:5] # placebo adjusted estimates at doses
#' vCov <- vcov(anovaMod)[2:5,2:5]
#' dose <- sort(unique(IBScovars$dose))[-1] # no estimate for placebo
#' ## candidate models
#' models <- Mods(emax = c(0.5, 1), betaMod=c(1,1), doses=c(0,4))
#' plot(models)
#' ## hand over placebo-adjusted estimates drFit to MCPMod
#' MM3 <- MCPMod(dose, drFit, S=vCov, models = models, type = "general",
#'               placAdj = TRUE, Delta=0.2)
#' plot(MM3, plotData="meansCI")
#'
#' ## The first example, but with critical value handed over
#' ## this is useful, e.g. in simulation studies
#' MM4 <- MCPMod(dose, resp, biom, models, Delta=0.5, critV = 2.31)
#' @export
MCPMod <- function(dose, resp, data = NULL, models = NULL, S=NULL,
                   type = c("normal", "general"), 
                   addCovars = ~1, placAdj = FALSE,
                   selModel = c("AIC", "maxT", "aveAIC"),
                   alpha = 0.025, df = NULL, critV = NULL,
                   doseType = c("TD", "ED"), Delta, p,
                   pVal = TRUE, alternative = c("one.sided", "two.sided"), 
                   na.action = na.fail, mvtcontrol = mvtnorm.control(),
                   bnds, control = NULL){

  direction <- attr(models, "direction")
  ## first perform multiple contrast test
  if(!is.null(data)){
    callMCT <- list(deparse(substitute(dose)), deparse(substitute(resp)), data,
                    models, S, type, addCovars, placAdj, alpha, df,
                    critV, pVal, alternative, na.action, mvtcontrol)
    test <- do.call(MCTtest, callMCT)
  } else {
    test <- MCTtest(dose, resp, data, models, S, type, addCovars, placAdj, alpha, df,
                    critV, pVal, alternative, na.action, mvtcontrol)
  }

  
  ## now pre-select models based on contrasts
  tstat <- test$tStat
  pvals <- attr(tstat, "pVal")
  if(!is.null(pvals)){
    tstat <- tstat[pvals < alpha]
  } else {
    tstat <- tstat[tstat > test$critVal]
  }
  if(length(tstat) == 0) ## stop if no model significant
    return(list(MCTtest = test, mods = NULL, modcrit = NULL, selMod = NULL, TD = NULL))

  ## fit models and calculate model selection criteria
  addArgs <- list(off=attr(models, "off"), scal=attr(models, "scal"))
  selModel <- match.arg(selModel)
  builtIn <- c("linlog", "linear", "quadratic", "linInt", "emax",
               "exponential", "logistic", "betaMod", "sigEmax")
  nams <- gsub("[0-9]", "", names(tstat)) ## remove numbers from model-names
  namsU <- unique(nams)
  
  mods <- vector("list", length(namsU));z <- 1
  if(missing(bnds)){
    if(!is.null(data)){
      cal <- as.character(match.call())
      doseVec <- data[, cal[2]]
    } else {
      doseVec <- dose
    }
    bnds <- defBnds(max(doseVec))
  } else {
    if(!is.list(bnds))
      stop("bnds needs to be a list")
  }
  if(selModel %in% c("AIC", "aveAIC")){
    if(type[1] == "normal"){
      modcrit <- AIC
    } else {
      modcrit <- gAIC
    }
  } else {
    modcrit <- function(x)
      max(tstat[attr(x, "model") == nams])
  }
  for(i in 1:length(namsU)){
    if(!is.null(data)){
      callMod <- list(deparse(substitute(dose)), deparse(substitute(resp)), data,
                      namsU[i], S, type, addCovars, placAdj, bnds[[namsU[i]]],
                      df, NULL, na.action, control, addArgs)
      mods[[i]] <- do.call(fitMod, callMod)
    } else {
      mods[[i]] <- fitMod(dose, resp, data, namsU[i], S, type, addCovars,
                          placAdj, bnds[[namsU[i]]], df, NULL,
                          na.action, control, addArgs)
    }
  }
  crit <- sapply(mods, modcrit)
  names(crit) <- names(mods) <- namsU
  attr(crit, "crit") <- selModel

  if(selModel %in% c("maxT", "AIC")){
    if(selModel == "AIC"){
      ind <- which.min(crit)
    }
    if(selModel == "maxT"){
      nam <- names(tstat)[which.max(tstat)]
      ind <- which(gsub("[0-9]", "", nam) == names(mods))
    }
    selMod <- namsU[ind] # name of selected model
  } else {
    aic <- crit-mean(crit)
    selMod <- exp(-0.5*aic)/sum(exp(-0.5*aic)) # model weights
    names(selMod) <- namsU
  }

  ## calculate target dose estimate
  tds <- NULL
  doseType <- match.arg(doseType)
  if(doseType == "TD"){
    if(missing(Delta))
      stop("\"Delta\" needs to be specified for TD estimation")
    tds <- sapply(mods, TD, Delta=Delta, direction = direction)
    attr(tds, "addPar") <- Delta
  }
  if(doseType == "ED"){
    if(missing(p))
      stop("\"p\" needs to be specified for TD estimation")
    tds <- sapply(mods, ED, p=p)
    attr(tds, "addPar") <- p
  }

  out <- list(MCTtest = test, mods = mods, modcrit=crit,
              selMod=selMod, doseEst=tds, doseType = doseType)
  class(out) <- "MCPMod"
  out
}

#' Predict from the fitted dose-response model
#'
#' @param object,x MCPMod object
#' @param predType,newdata,doseSeq,se.fit,...  predType determines whether predictions are returned for the full model
#'   (including potential covariates), the ls-means (SAS type) or the effect curve (difference to placebo).
#'
#'   newdata gives the covariates to use in producing the predictions (for \samp{predType = "full-model"}), if missing
#'   the covariates used for fitting are used.
#'
#'   doseSeq dose-sequence on where to produce predictions (for \samp{predType =
#'   "effect-curve"} and \samp{predType = "ls-means"}). If missing the doses used
#'   for fitting are used.
#'
#'   se.fit: logical determining, whether the standard error should be calculated.
#'
#'   \ldots: Additional arguments, for plot.MCPMod these are passed to plot.DRMod.
#'
#' @rdname MCPMod
#' @method predict MCPMod
#' @export
predict.MCPMod <- function(object,
                           predType = c("full-model", "ls-means", "effect-curve"),
                           newdata = NULL, doseSeq = NULL, se.fit = FALSE, ...){
  lapply(object$mods, function(x) predict(x, predType, newdata, doseSeq, se.fit))
}

#' @export
print.MCPMod <- function(x, digits=3, eps=1e-03, ...){
  cat("MCPMod\n")

  xx <- x$MCTtest
  cat("\nMultiple Contrast Test:\n")
  ord <- rev(order(xx$tStat))
  if (!any(is.null(attr(xx$tStat, "pVal")))) {
    pval <- format.pval(attr(xx$tStat, "pVal"), digits = digits, 
                        eps = eps)
    dfrm <- data.frame(round(xx$tStat, digits)[ord], pval[ord])
    names(dfrm) <- c("t-Stat", "adj-p")
  }
  else {
    dfrm <- data.frame(round(xx$tStat, digits)[ord])
    names(dfrm) <- c("t-Stat")
  }
  print(dfrm)
  if (!is.null(xx$critVal)) {
    twoSide <- xx$alternative == "two.sided"
    vec <- c(" one-sided)", " two-sided)")
    cat("\n", "Critical value: ", round(xx$critVal, digits), 
        sep = "")
    if (attr(xx$critVal, "Calc")) {
      cat(" (alpha = ", xx$alpha, ",", vec[twoSide + 1], 
          sep = "")
    }
    cat("\n")
  }
  cat("\n")

  cat("Estimated Dose Response Models:")
  for(i in 1:length(x$mods)){
    cat("\n")
    cat(names(x$mods)[i], "model\n")
    cofList <- coef(x$mods[[i]], sep = TRUE)
    cof <- do.call("c", cofList)
    namcof <- c(names(cofList$DRpars), names(cofList$covarPars))
    namcof <- gsub(" ", "", namcof)  # remove white spaces for GUI
    names(cof) <- gsub("doseM", "dose", namcof) # use more obvious names
    print(round(cof, digits))
  }
  if(attr(x$modcrit, "crit") != "aveAIC"){
    cat("\nSelected model (",attr(x$modcrit, "crit"),"): ", x$selMod, "\n", sep="")
  } else {
    cat("\nModel weights (AIC):\n")
    attr(x$selMod, "crit") <- NULL
    print(round(x$selMod, 4))
  }
  
  if(is.null(length(x$doseEst)))
    return()
  if(x$doseType == "TD")
    strn <- ", Delta="
  if(x$doseType == "ED")
    strn <- ", p="
  cat("\nEstimated ",x$doseType,strn,attr(x$doseEst, "addPar"),"\n", sep="")
  attr(x$doseEst, "addPar") <- NULL
  print(round(x$doseEst, 4))
}

#' @export
summary.MCPMod <- function(object, ...){
  class(object) <- "summary.MCPMod"
  print(object, digits = 3)
}

#' @export
print.summary.MCPMod <- function(x, ...){
  cat("MCPMod\n\n")

  cat(rep("*", 39), "\n", sep="")
  cat("MCP part \n")
  cat(rep("*", 39), "\n", sep="")
  print(x$MCTtest)
  cat("\n")

  if(length(x$mods) == 0)
    return()
  cat(rep("*", 39), "\n", sep="")
  cat("Mod part \n")
  cat(rep("*", 39), "\n", sep="")
  for(i in 1:length(x$mods)){
    if(i > 1)
      cat("\n")
    if(length(x$mods) > 1)
      cat("** Fitted model", i,"\n")
    summary(x$mods[[i]])
  }
  cat("\n")
  cat(rep("*", 39), "\n", sep="")
  cat("Model selection criteria (",attr(x$modcrit, "crit"),"):\n", sep="")
  cat(rep("*", 39), "\n", sep="")
  crit <- attr(x$modcrit, "crit")
  attr(x$modcrit, "crit") <- NULL
  print(x$modcrit)
  if(crit != "aveAIC"){
    cat("\nSelected model:", x$selMod, "\n")
  } else {
    cat("\nModel weights (AIC):\n")
    attr(x$selMod, "crit") <- NULL
    print(round(x$selMod, 4))
  }
  
  if(is.null(length(x$doseEst)))
    return()
  cat("\n")
  cat(rep("*", 39), "\n", sep="")
  if(x$doseType == "TD")
    strn <- ", Delta="
  if(x$doseType == "ED")
    strn <- ", p="
  cat("Estimated ",x$doseType,strn,attr(x$doseEst, "addPar"),"\n", sep="")
  cat(rep("*", 39), "\n", sep="")
  attr(x$doseEst, "addPar") <- NULL
  print(round(x$doseEst, 4))
}

#' Plot fitted dose-response model
#'
#' @inheritParams predict.MCPMod
#' @param CI,level,plotData,plotGrid,colMn,colFit Arguments for plot method: \samp{CI} determines whether confidence
#'   intervals should be plotted. \samp{level} determines the level of the confidence intervals. \samp{plotData}
#'   determines how the data are plotted: Either as means or as means with CI, raw data or none. In case of \samp{type =
#'   "normal"} and covariates the ls-means are displayed, when \samp{type = "general"} the option "raw" is not
#'   available.  \samp{colMn} and \samp{colFit} determine the colors of fitted model and the raw means.
#'
#' @rdname MCPMod
#' @method plot MCPMod
#' @export
plot.MCPMod <- function(x, CI = FALSE, level = 0.95,
                        plotData = c("means", "meansCI", "raw", "none"),
                        plotGrid = TRUE, colMn = 1, colFit = 1, ...){
  if(is.null(x$mods))
    stop("No models significant, nothing to plot")
  plotFunc(x, CI, level, plotData, plotGrid, colMn, colFit, ...)
}

