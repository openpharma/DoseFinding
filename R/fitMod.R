## functions related to fitting dose-response models using ML or generalized approach


#' Calculates default bounds for non-linear parameters in dose-response models
#'
#' Calculates reasonable bounds for non-linear parameters for the built-in non-linear regression model based on the dose
#' range under investigation.
#'
#' For the logistic model the first row corresponds to the ED50 parameter and the second row to the delta parameter. For
#' the sigmoid Emax model the first row corresponds to the ED50 parameter and the second row to the h parameter, while
#' for the beta model first and second row correspond to the delta1 and delta2 parameters. See \code{\link{logistic}},
#' \code{\link{sigEmax}} and \code{\link{betaMod}} for details.
#'
#'
#' @param mD Maximum dose in the study.
#' @param emax,exponential,logistic,sigEmax,betaMod values for the nonlinear parameters for these model-functions
#' @return List containing bounds for the model parameters.
#' @author Bjoern Bornkamp
#' @seealso \code{\link{fitMod}}
#' @examples
#'
#'   defBnds(mD = 1)
#'   defBnds(mD = 200)
#' @export
defBnds <- function(mD, emax = c(0.001, 1.5)*mD,
                    exponential = c(0.1, 2)*mD, 
                    logistic = matrix(c(0.001, 0.01, 1.5, 1/2)*mD, 2),
                    sigEmax = matrix(c(0.001*mD, 0.5, 1.5*mD, 10), 2),
                    betaMod = matrix(c(0.05,0.05,4,4), 2)){
  list(emax = emax, logistic = logistic, sigEmax = sigEmax,
       exponential = exponential, betaMod = betaMod)
}

#' Fit non-linear dose-response model
#'
#' Fits a dose-response model. Built-in dose-response models are "linlog", "linear", "quadratic", "emax", "exponential",
#' "sigEmax", "betaMod" and "logistic" (see \code{\link{drmodels}}).
#'
#' When \samp{type = "normal"} ordinary least squares is used and additional additive covariates can be specified in
#' \samp{addCovars}. The underlying assumption is hence normally distributed data and homoscedastic variance.
#'
#' For \samp{type = "general"} a generalized least squares criterion is used
#' \deqn{}{(f(dose,theta)-resp)'S^{-1}(f(dose,theta)-resp)}\deqn{
#' (f(dose,\theta)-resp)'S^{-1}(f(dose,\theta)-resp)}{(f(dose,theta)-resp)'S^{-1}(f(dose,theta)-resp)}
#' and an inverse weighting matrix is specified in \samp{S}, \samp{type =
#' "general"} is primarily of interest, when fitting a model to AN(C)OVA type
#' estimates obtained in a first stage fit, then \samp{resp} contains the estimates and \samp{S} is the estimated
#' covariance matrix for the estimates in \samp{resp}. Statistical inference (e.g. confidence intervals) rely on
#' asymptotic normality of the first stage estimates, which makes this method of interest only for sufficiently large
#' sample size for the first stage fit. A modified model-selection criterion can be applied to these model fits (see
#' also Pinheiro et al. 2014 for details).
#'
#' For details on the implemented numerical optimizer see the Details section below.
#'
#' Details on numerical optimizer for model-fitting:\cr For linear models fitting is done using numerical linear algebra
#' based on the QR decomposition.  For nonlinear models numerical optimization is performed only in the nonlinear
#' parameters in the model and optimizing over the linear parameters in each iteration (similar as the Golub-Pereyra
#' implemented in \code{\link{nls}}). For models with 1 nonlinear parameter the \code{\link{optimize}} function is used
#' for 2 nonlinear parameters the \code{\link{nlminb}} function is used. The starting value is generated using a
#' grid-search (with the grid size specified via \samp{control$gridSize}), or can directly be handed over via
#' \samp{start}.
#'
#' For details on the asymptotic approximation used for \samp{type = "normal"}, see Seber and Wild (2003, chapter 5).
#' For details on the asymptotic approximation used for \samp{type = "general"}, and the gAIC, see Pinheiro et al.
#' (2014).
#'
#' @aliases fitMod coef.DRMod vcov.DRMod predict.DRMod plot.DRMod logLik.DRMod AIC.DRMod gAIC gAIC.DRMod
#' @param dose,resp Either vectors of equal length specifying dose and response values, or names of variables in the
#'   data frame specified in \samp{data}.
#' @param data Data frame containing the variables referenced in dose and resp if \samp{data} is not specified it is
#'   assumed that \samp{dose} and \samp{resp} are variables referenced from data (and no vectors)
#' @param model The dose-response model to be used for fitting the data. Built-in models are "linlog", "linear",
#'   "quadratic", "emax", "exponential", "sigEmax", "betaMod" and "logistic" (see \link{drmodels}).
#' @param S The inverse weighting matrix used in case, when \samp{type = "general"}, see Description. For later
#'   inference statements (vcov or predict methods) it is assumed this is the estimated covariance of the estimates in
#'   the first stage fit.
#' @param type Determines whether inference is based on an ANCOVA model under a homoscedastic normality assumption (when
#'   \samp{type = "normal"}), or estimates at the doses and their covariance matrix and degrees of freedom are specified
#'   directly in \samp{resp}, \samp{S} and \samp{df}. See also the Description above and Pinheiro et al. (2014).
#' @param addCovars Formula specifying additional additive linear covariates (only for \samp{type = "normal"})
#' @param placAdj Logical, if true, it is assumed that placebo-adjusted
#'   estimates are specified in \samp{resp} (only possible for \samp{type =
#'   "general"}).
#' @param bnds Bounds for non-linear parameters. If missing the the default bounds from \code{\link{defBnds}} is used.
#'
#'   When the dose-response model has only one non-linear parameter (for example Emax or exponential model), \samp{bnds}
#'   needs to be a vector containing upper and lower bound. For models with two non-linear parameters \samp{bnds} needs
#'   to be a matrix containing the bounds in the rows, see the Description section of \code{\link{defBnds}} for details
#'   on the formatting of the bounds for the individual models.
#' @param df Degrees of freedom to use in case of \samp{type = "general"}. If this argument is missing \samp{df = Inf}
#'   is used. For \samp{type = "normal"} this argument is ignored as the exact degrees of freedom can be deduced from
#'   the model.
#' @param start Vector of starting values for the nonlinear parameters (ignored for linear models). When equal to NULL,
#'   a grid optimization is performed and the best value is used as starting value for the local optimizer.
#' @param na.action A function which indicates what should happen when the data contain NAs.
#' @param control A list with entries: "nlminbcontrol", "optimizetol" and "gridSize".
#'
#'   The entry nlminbcontrol needs to be a list and it is passed directly to control argument in the nlminb function,
#'   that is used internally for models with 2 nonlinear parameters.
#'
#'   The entry optimizetol is passed directly to the tol argument of the optimize function, which is used for models
#'   with 1 nonlinear parameters.
#'
#'   The entry gridSize needs to be a list with entries dim1 and dim2 giving the size of the grid for the gridsearch in
#'   1d or 2d models.
#' @param addArgs List containing two entries named "scal" and "off" for the "betaMod" and "linlog" model. When addArgs
#'   is NULL the following defaults is used \samp{list(scal = 1.2*max(doses), off = 0.01*max(doses))}.
#' @return An object of class DRMod. Essentially a list containing information about the fitted model coefficients, the
#'   residual sum of squares (or generalized residual sum of squares),
#' @author Bjoern Bornkamp
#' @seealso \code{\link{defBnds}}, \code{\link{drmodels}}
#' @references Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014) Model-based dose finding under model
#'   uncertainty using general parametric models, \emph{Statistics in Medicine}, \bold{33}, 1646--1661
#'
#'   Seber, G.A.F. and Wild, C.J. (2003). Nonlinear Regression, Wiley.
#' @examples
#'
#' ## Fit the emax model to the IBScovars data set
#' data(IBScovars)
#' fitemax <- fitMod(dose, resp, data=IBScovars, model="emax",
#'                   bnds = c(0.01, 4))
#'
#' ## methods for DRMod objects
#' summary(fitemax)
#' ## extracting coefficients
#' coef(fitemax)
#' ## (asymptotic) covariance matrix of estimates
#' vcov(fitemax)
#' ## predicting
#' newdat <- data.frame(dose = c(0,0.5,1), gender=factor(1))
#' predict(fitemax, newdata=newdat, predType = "full-model", se.fit = TRUE)
#' ## plotting
#' plot(fitemax, plotData = "meansCI", CI=TRUE)
#'
#' ## now include (additive) covariate gender
#' fitemax2 <- fitMod(dose, resp, data=IBScovars, model="emax",
#'                    addCovars = ~gender, bnds = c(0.01, 4))
#' vcov(fitemax2)
#' plot(fitemax2)
#' ## fitted log-likelihood
#' logLik(fitemax2)
#' ## extracting AIC (or BIC)
#' AIC(fitemax2)
#'
#' ## Illustrating the "general" approach for a binary regression
#' ## produce first stage fit (using dose as factor)
#' data(migraine)
#' PFrate <- migraine$painfree/migraine$ntrt
#' doseVec <- migraine$dose
#' doseVecFac <- as.factor(migraine$dose)
#' ## fit logistic regression with dose as factor
#' fitBin <- glm(PFrate~doseVecFac-1, family = binomial,
#'               weights = migraine$ntrt)
#' drEst <- coef(fitBin)
#' vCov <- vcov(fitBin)
#' ## now fit an Emax model (on logit scale)
#' gfit <- fitMod(doseVec, drEst, S=vCov, model = "emax", bnds = c(0,100),
#'                 type = "general")
#' ## model fit on logit scale
#' plot(gfit, plotData = "meansCI", CI = TRUE)
#' ## model on probability scale
#' logitPred <- predict(gfit, predType ="ls-means", doseSeq = 0:200,
#'                      se.fit=TRUE)
#' plot(0:200, 1/(1+exp(-logitPred$fit)), type = "l", ylim = c(0, 0.5),
#'      ylab = "Probability of being painfree", xlab = "Dose")
#' LB <- logitPred$fit-qnorm(0.975)*logitPred$se.fit
#' UB <- logitPred$fit+qnorm(0.975)*logitPred$se.fit
#' lines(0:200, 1/(1+exp(-LB)))
#' lines(0:200, 1/(1+exp(-UB)))
#'
#'
#' ## now illustrate "general" approach for placebo-adjusted data (on
#' ## IBScovars) note that the estimates are identical to fitemax2 above)
#' anovaMod <- lm(resp~factor(dose)+gender, data=IBScovars)
#' drFit <- coef(anovaMod)[2:5] # placebo adjusted estimates at doses
#' vCov <- vcov(anovaMod)[2:5,2:5]
#' dose <- sort(unique(IBScovars$dose))[-1]
#' ## now fit an emax model to these estimates
#' gfit2 <- fitMod(dose, drFit, S=vCov, model = "emax", type = "general",
#'                placAdj = TRUE, bnds = c(0.01, 2))
#' ## some outputs
#' summary(gfit2)
#' coef(gfit2)
#' vcov(gfit2)
#' predict(gfit2, se.fit = TRUE, doseSeq = c(1,2,3,4), predType = "effect-curve")
#' plot(gfit2, CI=TRUE, plotData = "meansCI")
#' gAIC(gfit2)
#'
#' @export
fitMod <- function(dose, resp, data = NULL, model = NULL, S = NULL,
                   type = c("normal", "general"),
                   addCovars = ~1, placAdj = FALSE, bnds, df = NULL,
                   start = NULL, na.action = na.fail, control = NULL,
                   addArgs = NULL){
  ## check for valid dose, resp and data
  cal <- as.character(match.call())
  type <- match.arg(type)
  lst <- checkAnalyArgs(dose, resp, data, S, type,
                        addCovars, placAdj, na.action, cal)
  doseNam <- lst$doseNam;respNam <- lst$respNam
  dose <- lst$dd[[doseNam]];type <- lst$type
  resp <- lst$dd[[respNam]];data <- lst$dd;S <- lst$S
  covarsUsed <- addCovars != ~1
  
  ## check type related arguments
  if(type == "general"){
    if(placAdj & model %in% c("linlog", "logistic")) # stop as fitting algorithm assumes f^0(0) = 0
      stop("logistic and linlog models cannot be fitted to placebo adjusted data") 
    if(covarsUsed)
      stop("addCovars argument ignored for type == \"general\"")
    if(is.null(df))
      df <- Inf
  }
  ## check whether model has been specified correctly
  builtIn <- c("linlog", "linear", "quadratic", "linInt", "emax",
               "exponential", "logistic", "betaMod", "sigEmax")
  if(missing(model))
    stop("Need to specify the model that should be fitted")
  modelNum <- match(model, builtIn)
  if(is.na(modelNum))
    stop("Invalid dose-response model specified")
  ## check for start argument
  if(modelNum < 5 & !is.null(start))
    message("Message: Starting values in \"start\" ignored for linear models")
  ## check for valid bnds
  if(modelNum > 4){
    if(missing(bnds)){
      message("Message: Need bounds in \"bnds\" for nonlinear models, using default bounds from \"defBnds\".")
      bnds <- defBnds(max(dose))[[model]]
    } else {
      if(is.null(bnds)){
        message("Message: Need bounds in \"bnds\" for nonlinear models, using default bounds from \"defBnds\".")
        bnds <- defBnds(max(dose))[[model]]
      }
    }
  }
  ## addArgs argument
  scal <- off <- nodes <- NULL
  if(model %in% c("linlog", "betaMod")){
    aPar <- getAddArgs(addArgs, sort(unique(dose)))
    if(model == "betaMod")
      scal <- aPar$scal
    if(model == "linlog")
      off <- aPar$off
  }
  if(model == "linInt"){ ## not allowed to use nodes different from used doses
    nodes <- sort(unique(dose))
  }

  ## call fit-model raw!
  out <- fitMod.raw(dose, resp, data, model, S, type,
                    addCovars, placAdj, bnds, df, start,
                    na.action, control, doseNam=doseNam,
                    respNam=respNam, off = off, scal = scal,
                    nodes=nodes, covarsUsed)
  ## attach data to object
  reord <- order(lst$ord)
  if(type == "normal"){
    if(covarsUsed){
      attr(out, "data") <- data[reord,]
    } else {
      dat <- data.frame(dose=dose, resp=resp)
      colnames(dat) <- c(doseNam, respNam)
      attr(out, "data") <- dat[reord,]
    }
  } else {
    lst <- list(dose=dose[reord], resp=resp[reord], S=S[reord,reord]) 
    names(lst) <- c(doseNam, respNam, "S")
    attr(out, "data") <- lst
  }
  out
}

#' @export
print.DRMod <- function(x, digits = 4, ...){
  if (length(x) == 1) {
    cat("NA\n")
    return()
  }
  cat("Dose Response Model\n\n")
  cat(paste("Model:", attr(x, "model")), "\n")
  cat(paste("Fit-type:", attr(x, "type")), "\n\n")
  Coefs <- sepCoef(x)
  cat("Coefficients dose-response model\n")
  print(signif(Coefs$DRpars, digits))
  if(attr(x, "type") == "normal"){
    if(x$addCovars != ~1){
      cat("Coefficients additional covariates\n")
      print(signif(Coefs$covarPars, digits))
    }
    cat("\nDegrees of freedom:", x$df, "\n")
    cat("Residual standard error:",
        signif(sqrt(x$RSS/x$df), digits),"\n")
  }
  if(attr(x, "type") == "general"){
    cat("\nFitted to:\n")
    doseRespNam <- attr(x, "doseRespNam")
    resp <- attr(x, "data")[[doseRespNam[2]]]
    names(resp) <- attr(x, "data")[[doseRespNam[1]]]
    print(signif(resp, digits))
    cat("\nGeneralized residual sum of squares:",
        signif(x$gRSS, digits),"\n")
  }
}

#' @export
summary.DRMod <- function(object, digits = 3, ...){
  class(object) <- "summary.DRMod"
  print(object, digits = digits)
}

#' @export
print.summary.DRMod <- function(x, digits = 3, data, ...){
  if(length(x) == 1){
    cat("NA\n")
    return()
  }
  data <- attr(x, "data")
  cat("Dose Response Model\n\n")
  cat(paste("Model:", attr(x, "model")), "\n")
  type <- attr(x, "type")
  cat(paste("Fit-type:", type), "\n")
  if(type == "normal"){
    ## residual information
    cat("\nResiduals:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    respNam <- attr(x, "doseRespNam")[2]
    resid <- predict.DRMod(x, predType = "full-model")-data[[respNam]]
    rq <- structure(quantile(resid), names = nam)
    print(rq, digits = digits, ...)
  }
  cat("\nCoefficients with approx. stand. error:\n")
  coefs <- x$coef
  sdv <- sqrt(diag(vcov.DRMod(x)))
  datf <- matrix(nrow = length(coefs), ncol = 2)
  datf[,1] <- coefs
  datf[,2] <- sdv
  colnam <- c("Estimate", "Std. Error")
  dimnames(datf) <- list(names(coefs), colnam)
  print(datf, digits = digits)
  if(type == "normal"){
    cat("\nResidual standard error:",
        signif(sqrt(x$RSS/x$df), digits), "\n")
    cat("Degrees of freedom:", x$df, "\n")
  }
  if(type == "general"){
    doseRespNam <- attr(x, "doseRespNam")
    dose <- attr(x, "data")[[doseRespNam[1]]]
    drEst <- attr(x, "data")[[doseRespNam[2]]]
    names(drEst) <- dose
    S <- attr(x, "data")$S
    dimnames(S) <- list(dose, dose)
    cat("\nFitted to:\n")
    print(signif(drEst, digits))
    cat("\nwith Covariance Matrix:\n")
    print(signif(S, digits))
  }
}

#' Extract dose-response model coefficients
#'
#' @param object,x DRMod object
#' @param sep Logical determining whether all coefficients should be returned in one numeric or separated in a list.
#' @param ... Additional arguments for plotting for the plot method. For all other cases additional arguments are
#'   ignored.
#'
#' @rdname fitMod
#' @method coef DRMod
#' @export
#' 
coef.DRMod <- function(object, sep = FALSE, ...){
  if(length(object) == 1){ # object does not contain a converged fit
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }
  if(sep){
    return(sepCoef(object))
  }
  object$coefs
}

#' Extract dose-response vcov matrix
#'
#' @rdname fitMod
#' @method vcov DRMod
#' @export
vcov.DRMod <- function(object, ...){
  ## object - DRMod object
  ## uGrad - function returning gradient for userModel
  if(length(object) == 1){ # object does not contain a converged fit
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }
  type <- attr(object, "type")
  model <- attr(object, "model")
  cf <- sepCoef(object)$DRpars
  nams <- names(coef(object))
  scal <- attr(object, "scal")
  off <- attr(object, "off")
  nodes <- attr(object, "nodes")  
  doseNam <- attr(object, "doseRespNam")[1]
  if(type == "normal"){
    addCovars <- object$addCovars
    xlev <- attr(object, "xlev")
    RSS <- object$RSS
    df <- object$df
    data <- attr(object, "data")
    dose <- attr(object, "data")[[doseNam]]
    m <- model.matrix(addCovars, data, xlev = xlev)
  }
  if(type == "general"){
    placAdj <- attr(object, "placAdj")
    if(placAdj) # no intercept
      cf <- c(0, cf)
    dose <- attr(object, "data")[[doseNam]]
    inS <- solve(attr(object, "data")$S)
  }
  grd <- gradCalc(model, cf, dose, off, scal, nodes)
  if(type == "normal"){
    J <- cbind(grd, m[,-1])
    JtJ <- crossprod(J)
    covMat <- try(solve(JtJ)*RSS/df, silent=TRUE)
    if(!inherits(covMat, "matrix")){
      covMat <- try(chol2inv(qr.R(qr(J)))*RSS/df, silent=TRUE) # more stable (a little slower)
      if(!inherits(covMat, "matrix")){
        warning("cannot calculate covariance matrix. singular matrix in calculation of covariance matrix.")
        nrw <- length(grd[1,])
        covMat <- matrix(NA, nrow=nrw, ncol=nrw)
      }
      dimnames(covMat) <- dimnames(JtJ)
    }
  }
  if(type == "general"){
    if(placAdj){
      if(model != "linInt")
        grd <- grd[,-1]
    }
    covMat <- try(solve(t(grd)%*%inS%*%grd), silent = TRUE)
    if(!inherits(covMat, "matrix")) {
      warning("cannot calculate covariance matrix. singular matrix in calculation of covariance matrix.")
      nrw <- length(grd[1,])
      covMat <- matrix(NA, nrow=nrw, ncol=nrw)
    }
    
  }
  dimnames(covMat) <- list(nams, nams)
  covMat
}

#'Make predictions from dose-response model
#'
#'@param predType,newdata,doseSeq,se.fit predType determines whether predictions are returned for the full model
#'  (including potential covariates), the ls-means (SAS type) or the effect curve (difference to placebo).
#'
#'  newdata gives the covariates to use in producing the predictions (for predType = "full-model"), if missing the
#'  covariates used for fitting are used.
#'
#'  doseSeq dose-sequence on where to produce predictions (for predType = "effect-curve" and predType = "ls-means"). If
#'  missing the doses used for fitting are used.
#'
#'  se.fit: logical determining, whether the standard error should be calculated.
#'
#'@rdname fitMod
#'@method predict DRMod
#'@export
predict.DRMod <- function(object, predType = c("full-model", "ls-means", "effect-curve"),
                          newdata = NULL, doseSeq = NULL, se.fit = FALSE, ...){
  ## Extract relevant information from object
  scal <- attr(object, "scal")
  off <- attr(object, "off")
  nodes <- attr(object, "nodes")
  model <- attr(object, "model")
  addCovars <- attr(object, "addCovars")
  xlev <- attr(object, "xlev")
  doseNam <- attr(object, "doseRespNam")[1]
  data <- attr(object, "data")
  type <- attr(object, "type")

  if(missing(predType))
    stop("need to specify the type of prediction in \"predType\"")
  predType <- match.arg(predType)
  ## if model fitted on plac-adj. data can only produce predictions for effect-curve
  if(attr(object, "placAdj") & predType != "effect-curve"){ 
    message("Message: Setting predType to \"effect-curve\" for placebo-adjusted data")
    predType <- "effect-curve"
  }
  if(type == "general" & predType == "full-model"){ ## there are no covariates
    message("Message: Setting predType to \"ls-means\" for \"type = general\"")
    predType <- "ls-means"
  }
  
  if(predType %in% c("ls-means", "full-model")){
    ## create design-matrix according to the SAS predType ls-means
    if(predType == "ls-means"){
      if(!is.null(newdata))
        stop("newdata is ignored for \"predType = \"ls-means\"")
      if(is.null(doseSeq)){ ## use doses used for fitting
        if(type == "normal")
          doseSeq <- data[, doseNam]
        if(type == "general")
          doseSeq <- data[[doseNam]]
      }
      covarsUsed <- addCovars != ~1
      if(covarsUsed){
        nams <- all.vars(addCovars)
        out <- list()
        z <- 1
        for(covar in nams){
          varb <- data[,covar]
          if(is.numeric(varb)){
            out[[z]] <- mean(varb)
          } else if(is.factor(varb)){
            k <- nlevels(varb)
            out[[z]] <- rep(1/k, k-1)
          }
          z <- z+1
        }
        out <- do.call("c", out)
        m <- matrix(rep(out, length(doseSeq)), byrow=TRUE, nrow = length(doseSeq))
      }
    }
    ## create design-matrix either from newdata or data used for fitting
    if(predType == "full-model"){
      if(!is.null(doseSeq) & predType == "full-model")
        stop("doseSeq should only be used when predType = \"effect-curve\" or \"ls-means\"")
      if(is.null(newdata)){
        ## if not provided use covariates in observed data
        if(type == "normal"){
          m <- model.matrix(addCovars, data)
          doseSeq <- data[, doseNam]
        } else {
          doseSeq <- data[[doseNam]]
        }
      } else {
        tms <- c(doseNam, attr(terms(addCovars), "term.labels"))
        missind <- !is.element(tms, names(newdata))
        if(any(missind)){
          chct <- paste("No values specified in newdata for", paste(tms[missind], collapse=", "))
          stop(chct)
        } else {
          m <- model.matrix(addCovars, newdata, xlev = xlev)
          doseSeq <- newdata[, doseNam]
          if(nrow(m) != length(doseSeq))
            stop("incompatible model matrix and doseSeq created from newdata")
        } 
      }
      m <- m[,-1, drop=FALSE] # remove intercept column (is necessary)
    }
    coeflist <- sepCoef(object) # separate coefs of DR model and additional covars
    DRpars <- coeflist$DRpars   
    covarPars <- coeflist$covarPars
    ## predictions
    if(model != "linInt"){
      call <- c(list(doseSeq), as.list(c(DRpars, scal, off)))
    } else {
      call <- c(list(doseSeq), as.list(list(DRpars, nodes)))
    }
    mn <- do.call(model, call)
    if(addCovars != ~1)
      mn <- mn + as.numeric(m%*%covarPars)
    if(!se.fit){
      return(as.numeric(mn))
    } else { ## calculate standard error of predictions
      covMat <- vcov(object)
      if(any(is.na(covMat))){
        seFit <- (rep(NA, length(doseSeq)))
      } else {
        grd <- gradCalc(model, DRpars, doseSeq, off, scal, nodes)
        if(addCovars != ~1)
          grd <- cbind(grd, m)
        cholcovMat <- try(chol(covMat), silent = TRUE)
        if (!inherits(cholcovMat, "matrix")) {
          warning("Cannot cannot calculate standard deviation for ", 
                  model, " model.\n")
          seFit <- rep(NA, length(doseSeq))
        } else {
          seFit <- sqrt(rowSums((grd%*%t(cholcovMat))^2)) # t(grd)%*%covMat%*%grd
        }
      }
      return(list(fit = mn, se.fit = as.vector(seFit)))
    }
  }
  if(predType == "effect-curve") {  ## predict effect curve
    if(!is.null(newdata))
      stop("newdata is ignored for \"predType = \"effect-curve\"")
    if(is.null(doseSeq)){
      if(type == "normal")
        doseSeq <- data[, doseNam]
      if(type == "general")
        doseSeq <- data[[doseNam]]
    }
    coeflist <- sepCoef(object) 
    DRpars <- coeflist$DRpars   
    if(attr(object, "placAdj")){
      DRpars <- c(0, DRpars)
      if(model == "linInt")
        nodes <- c(0, nodes)
    } else {
      if(model != "linInt"){
        DRpars[1] <- 0
      } else {
        DRpars <- DRpars - DRpars[1]
      }
    }
    ## predictions
    if(model != "linInt"){
      call <- c(list(doseSeq), as.list(c(DRpars, scal, off)))
    } else {
      call <- c(list(doseSeq), as.list(list(DRpars, nodes)))
    }
    mn <- do.call(model, call)
    if(is.element(model,c("logistic", "linlog"))){ # if standardized model not 0 at placebo
      call <- c(0, as.list(c(DRpars, scal, off)))      
      predbase <- do.call(model, call)
      mn <- mn-predbase
    }
    if(!se.fit){
      return(as.numeric(mn))
    } else { ## calculate st. error (no need to calculate full covMat here)
      covMat <- vcov(object)
      if(addCovars != ~1) ## remove columns corresponding to covariates
        covMat <- covMat[1:length(DRpars), 1:length(DRpars)]
      if(!attr(object, "placAdj")){ ## remove intercept from cov-matrix
        if(model != "linInt"){
          covMat <- covMat[-1,-1]
        } else {
          diffMat <- cbind(-1,diag(length(DRpars)-1))
          covMat <- diffMat%*%covMat%*%t(diffMat)
        }
      }
      if(any(is.na(covMat))){
        seFit <- (rep(NA, length(doseSeq)))
      } else {
        grd <- gradCalc(model, DRpars, doseSeq, off, scal, nodes)
        if(!is.matrix(grd)){ # can happen if length(doseSeq) == 1
          grd <- matrix(grd, nrow = 1)
        }
        if(model == "linInt"){
          grd <- grd[,-1, drop = FALSE]
        } else {
          grd0 <- gradCalc(model, DRpars, 0, off, scal, nodes)
          grd <- grd[, -1, drop=FALSE]
          grd0 <- grd0[,-1]
          grd <- sweep(grd, 2, grd0, "-")
        }
        cholcovMat <- try(chol(covMat), silent = TRUE)
        if (!inherits(cholcovMat, "matrix")) {
          warning("Cannot cannot calculate standard deviation for ", 
                  model, " model.\n")
          seFit <- rep(NA, length(doseSeq))
        } else {
          seFit <- sqrt(rowSums((grd%*%t(cholcovMat))^2)) # t(grd)%*%covMat%*%grd
        }
      }
      res <- list(fit = mn, se.fit = as.vector(seFit))
      return(res)
    }    
  }
}

## plot.DRMod <- function(x, CI = FALSE, level = 0.95,
##                        plotData = c("means", "meansCI", "none"),
##                        lenDose = 201, ...){
##   ## arguments passed to plot
##   pArgs <- list(...)
##   ## Extract relevant information from object
##   scal <- attr(x, "addArgs")$scal
##   off <- attr(x, "addArgs")$off
##   model <- attr(x, "model")
##   addCovars <- attr(x, "addCovars")
##   covarsUsed <- addCovars != ~1
##   xlev <- attr(x, "xlev")
##   doseNam <- attr(x, "doseRespNam")[1]
##   respNam <- attr(x, "doseRespNam")[2]
##   data <- attr(x, "data")
##   type <- attr(x, "type")
##   placAdj <- attr(x, "placAdj")
##   doseSeq <- seq(0, max(data[[doseNam]]), length=lenDose)

##   plotData <- match.arg(plotData)
##   if(type == "normal"){
##     ## first produce estimates for ANOVA type model
##     if(plotData %in% c("means", "meansCI")){
##       data$doseFac <- as.factor(data[[doseNam]])
##       form <- as.formula(paste(respNam, "~ doseFac +", addCovars[2]))
##       fit <- lm(form, data=data)
##       ## build design matrix for prediction
##       dose <- sort(unique(data[[doseNam]]))
##       preddat <- data.frame(doseFac=factor(dose))
##       m <- model.matrix(~doseFac, data=preddat)
##       if(covarsUsed){
##         ## get sas type ls-means
##         nams <- all.vars(addCovars)
##         out <- list()
##         z <- 1
##         for(covar in nams){
##           varb <- data[,covar]
##           if(is.numeric(varb)){
##             out[[z]] <- mean(varb)
##           } else if(is.factor(varb)){
##             k <- nlevels(varb)
##             out[[z]] <- rep(1/k, k-1)
##           }
##           z <- z+1
##         }
##         out <- do.call("c", out)
##         m0 <- matrix(rep(out, length(dose)), byrow=TRUE, nrow = length(dose))
##         m <- cbind(m, m0)
##       }
##       mns <- as.numeric(m%*%coef(fit))
##       lbndm <- ubndm <- rep(NA, length(mns))
##       if(plotData == "meansCI"){
##         sdv <- sqrt(diag(m%*%vcov(fit)%*%t(m)))
##         quant <- qt(1 - (1 - level)/2, df=x$df)
##         lbndm <- mns-quant*sdv
##         ubndm <- mns+quant*sdv
##       }
##     }
##   }
##   if(type == "general"){
##     ## extract ANOVA estimates
##     if(plotData %in% c("means", "meansCI")){
##       dose <- data[[doseNam]]
##       mns <- data[[respNam]]
##       sdv <- sqrt(diag(data$S))
##       lbndm <- ubndm <- rep(NA, length(dose))
##       if(plotData == "meansCI"){
##         quant <- qnorm(1 - (1 - level)/2)
##         lbndm <- mns-quant*sdv
##         ubndm <- mns+quant*sdv
##       }
##     }
##   }
##   ## curve produced (use "ls-means" apart when data are fitted on placAdj scale)
##   predtype <- ifelse(placAdj, "effect-curve", "ls-means")
##   predmn <- predict(x, doseSeq = doseSeq, predType = predtype, se.fit = CI)
##   lbnd <- ubnd <- rep(NA, length(doseSeq))
##   if(CI){
##     quant <- qt(1 - (1 - level)/2, df=x$df)
##     lbnd <- predmn$fit-quant*predmn$se.fit
##     ubnd <- predmn$fit+quant*predmn$se.fit
##     predmn <- predmn$fit
##   }
##   ## determine plotting range
##   if(plotData %in% c("means", "meansCI")){
##     rng <- range(lbndm, ubndm, mns, predmn, ubnd, lbnd, na.rm=TRUE)
##   } else {
##     rng <- range(predmn, ubnd, lbnd, na.rm=TRUE)    
##   }
##   dff <- diff(rng)
##   ylim <- c(rng[1] - 0.02 * dff, rng[2] + 0.02 * dff)
##   ## default title
##   main <- "Dose-Response Curve"
##   main2 <- ifelse(placAdj, "(placebo-adjusted)", "(ls-means)")
##   main <- paste(main, main2)
##   ## plot
##   callList <- list(doseSeq, predmn, type = "l", col = "white",
##                    xlab = doseNam, ylim = ylim,
##                    ylab = respNam, main = main)
##   callList[names(pArgs)] <- pArgs
##   do.call("plot", callList)
##   grid()
##   if(plotData %in% c("means", "meansCI")){
##     points(dose, mns, pch = 19, cex = 0.75)
##     if(plotData == "meansCI"){
##       for(i in 1:length(dose)){
##         lines(c(dose[i],dose[i]), c(lbndm[i], ubndm[i]), lty=2)
##       }
##     }
##   }
##   lines(doseSeq, predmn, lwd=1.5)
##   lines(doseSeq, ubnd, lwd=1.5)
##   lines(doseSeq, lbnd, lwd=1.5)
## }


#'Plot fitted dose-response model
#'
#'@param CI,level,plotData,plotGrid,colMn,colFit Arguments for plot method: \samp{CI} determines whether confidence
#'  intervals should be plotted. \samp{level} determines the level of the confidence intervals. \samp{plotData}
#'  determines how the data are plotted: Either as means or as means with CI, raw data or none. In case of \samp{type =
#'  "normal"} and covariates the ls-means are displayed, when \samp{type = "general"} the option "raw" is not available.
#'  \samp{colMn} and \samp{colFit} determine the colors of fitted model and the raw means.
#'
#'@rdname fitMod
#'@method plot DRMod
#'@export
plot.DRMod <- function(x, CI = FALSE, level = 0.95,
                       plotData = c("means", "meansCI", "raw", "none"),
                       plotGrid = TRUE, colMn = 1, colFit = 1, ...){
  plotFunc(x, CI, level, plotData, plotGrid, colMn, colFit, ...)
}


#' Extract log-likelihood of dose-response model
#'
#' @rdname fitMod
#' @method logLik DRMod
#' @export
logLik.DRMod <- function(object, ...){

  type <- attr(object, "type")
  data <- attr(object, "data")
  if(type == "normal"){
    RSS <- object$RSS
    n <- nrow(data)
    sig2 <- RSS/n
    val <- -n/2*(log(2*pi) + 1 + log(sig2))
    attr(val, "df") <- length(object$coefs)+1 # +1 because of sigma parameter
    class(val) <- "logLik"
    return(val)
  }
  if(type == "general")
    stop("method glogLik only available for type == \"normal\"")
}

#' Extract AIC of dose-response model
#'
#' @param k Penalty to use for model-selection criterion (AIC uses 2, BIC uses log(n)).
#'
#' @rdname fitMod
#' @method AIC DRMod
#' @export
AIC.DRMod <- function(object, ..., k = 2){
  type <- attr(object, "type")
  if(type == "general")
    stop("use method gAIC for type == \"general\"")
  logL <- logLik(object)
  -2*as.vector(logL) + k*(attr(logL, "df")) 
}


#' Extract gAIC of dose-response model
#'
#' @rdname fitMod
#' @method gAIC DRMod
#' @export
gAIC.DRMod <- function(object, ..., k = 2){
  type <- attr(object, "type")
  if(type == "normal")
    stop("use method AIC for type == \"normal\"")
  object$gRSS+k*length(object$coefs)
}




