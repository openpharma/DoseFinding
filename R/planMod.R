## various functions for assessing the operating characteristics of a design
## for model-based estimation of dose-response functions

#' Evaluate performance metrics for fitting dose-response models
#'
#' This function evaluates, the performance metrics for fitting dose-response models (using asymptotic approximations or
#' simulations). Note that some metrics are available via the print method and others only via the summary
#' method applied to planMod objects. The implemented metrics are \itemize{
#' \item Root of the mean-squared error to estimate the placebo-adjusted
#' dose-response averaged over the used dose-levels, i.e. a rather discrete set
#' (`dRMSE`). Available via the print method of planMod objects.  \item
#' Root of the mean-squared error to estimate the placebo-adjusted
#' dose-response (`cRMSE`) averaged over fine (almost continuous) grid at
#' 101 equally spaced values between placebo and the maximum dose. NOTE:
#' Available via the summary method applied to planMod objects.  \item Ratio of
#' the placebo-adjusted mean-squared error (at the observed doses) of
#' model-based vs ANOVA approach (`Eff-vs-ANOVA`). This can be interpreted
#' on the sample size scale. NOTE: Available via the summary method applied to
#' planMod objects.  \item Power that the (unadjusted) one-sided \samp{1-alpha}
#' confidence interval comparing the dose with maximum effect vs placebo is
#' larger than \samp{tau}. By default \samp{alpha = 0.025} and \samp{tau = 0}
#' (`Pow(maxDose)`). Available via the print method of planMod objects.
#' \item Probability that the EDp estimate is within the true \[EDpLB, EDpUB\]
#' (by default \samp{p=0.5}, \samp{pLB=0.25} and \samp{pUB=0.75}). This metric
#' gives an idea on the ability to characterize the increasing part of the
#' dose-response curve (`P(EDp)`). Available via the print method of
#' planMod objects.  \item Length of the quantile range for a target dose (TD
#' or EDp). This is calculated by taking the difference of the dUB and dLB
#' quantile of the empirical distribution of the dose estimates.
#' (`lengthTDCI` and `lengthEDpCI`). It is NOT calculated by
#' calculating confidence interval lengths in each simulated data-set and
#' taking the mean. NOTE: Available via the summary method of planMod objects.
#' }
#'
#' A plot method exists to summarize dose-response and dose estimations graphically.
#'
#'
#' @aliases planMod plot.planMod summary.planMod
#' @param model Character vector determining the dose-response model(s) to be used for fitting the data.  When more than
#'   one dose-response model is provided the best fitting model is chosen using the AIC. Built-in models are "linlog",
#'   "linear", "quadratic", "emax", "exponential", "sigEmax", "betaMod" and "logistic" (see [drmodels]).
#' @param altModels An object of class \samp{Mods}, defining the true mean vectors under which operating characteristics
#'   should be calculated.
#' @param n,sigma,S Either a vector \samp{n} and \samp{sigma} or \samp{S} need to be specified.  When \samp{n} and
#'   \samp{sigma} are specified it is assumed computations are made for a normal homoscedastic ANOVA model with group
#'   sample sizes given by \samp{n} and residual standard deviation \samp{sigma}, i.e. the covariance matrix used for
#'   the estimates is thus `sigma^2*diag(1/n)` and the degrees of freedom are calculated as
#'   `sum(n)-nrow(contMat)`. When a single number is specified for \samp{n} it is assumed this is the sample size
#'   per group and balanced allocations are used.\cr
#'
#'   When \samp{S} is specified this will be used as covariance matrix for the estimates.
#' @param doses Doses to use
#' @param asyApprox,simulation Logicals determining, whether asymptotic approximations or simulations should be
#'   calculated. If multiple models are specified in \samp{model} asymptotic approximations are not available.
#' @param alpha,tau Significance level for the one-sided confidence interval for model-based contrast of best dose vs
#'   placebo. Tau is the threshold to compare the confidence interval limit to. CI(MaxDCont) gives the percentage that
#'   the bound of the confidence interval was larger than tau.
#' @param p,pLB,pUB p determines the type of EDp to estimate. pLB and pUB define the bounds for the EDp estimate. The
#'   performance metric Pr(Id-ED) gives the percentage that the estimated EDp was within the true EDpLB and EDpUB.
#' @param nSim Number of simulations
#' @param cores Number of cores to use for simulations. By default 1 cores is used, note that cores > 1 will have no
#'   effect Windows, as the mclapply function is used internally.
#' @param showSimProgress In case of simulations show the progress using a progress-bar.
#' @param bnds Bounds for non-linear parameters. This needs to be a list with list entries corresponding to the selected
#'   bounds. The names of the list entries need to correspond to the model names. The [defBnds()] function
#'   provides the default selection.
#' @param addArgs See the corresponding argument in function [fitMod()]. This argument is directly passed to
#'   fitMod.
#' @author Bjoern Bornkamp
#' @seealso [fitMod()]
#' @references Pinheiro, J.C. and Bornkamp, B. (2017). Designing phase II dose-finding studies: sample size, 
#' doses, and dose allocation weights.  Handbook of Methods for Designing, Monitoring, and Analyzing Dose-Finding 
#' Trials, Chapter 13, Chapman and Hall/CRC, 229-246.
#'
#' @examples
#'
#' \dontrun{
#' doses <- c(0,10,25,50,100,150)
#' fmodels <- Mods(linear = NULL, emax = 25,
#'                 logistic = c(50, 10.88111), exponential= 85,
#'                 betaMod=rbind(c(0.33,2.31),c(1.39,1.39)),
#'                 doses = doses, addArgs=list(scal = 200),
#'                 placEff = 0, maxEff = 0.4)
#' sigma <- 1
#' n <- rep(62, 6)*2
#'
#' model <- "quadratic"
#' pObj <- planMod(model, fmodels, n, sigma, doses=doses,
#'                simulation = TRUE,
#'                alpha = 0.025, nSim = 200,
#'                p = 0.5, pLB = 0.25, pUB = 0.75)
#' print(pObj)
#' ## to get additional metrics (e.g. Eff-vs-ANOVA, cRMSE, lengthTDCI, ...)
#' summary(pObj, p = 0.5, Delta = 0.3)
#' plot(pObj)
#' plot(pObj, type = "TD", Delta=0.3)
#' plot(pObj, type = "ED", p = 0.5)
#' }
#' @export
planMod <- function(model, altModels, n, sigma, S, doses,
                    asyApprox = TRUE, simulation = FALSE,
                    alpha = 0.025, tau = 0,
                    p = 0.5, pLB = 0.25, pUB = 0.75,
                    nSim = 100, cores = 1,  showSimProgress = TRUE,
                    bnds, addArgs = NULL){
  if(any(is.element(model, "linInt")))
    stop("planMod works for all built-in models but not linInt")
  if(length(model) > 1 & asyApprox){
    stop("\"asyApprox\" needs to be FALSE for multiple models")
  } 
  ## off and scal
  off <- scal <- NULL
  if(any(is.element(model, c("linlog", "betaMod")))) {
    lst <- getAddArgs(addArgs, sort(unique(doses)))
    if ("betaMod" %in% model) 
      scal <- lst$scal
    if ("linlog" %in% model) 
      off <- lst$off
  }
  if(missing(doses))
    doses <- attr(altModels, "doses")
  ## calculate mean response at doses
  muMat <- getResp(altModels, doses)

  nD <- length(doses)
  if(missing(S)){
    if(missing(n) | missing(sigma))
      stop("either S or n and sigma need to be specified")
    if (length(n) == 1) 
      n <- rep(n, nD)
    if (length(n) != nD) 
      stop("\"n\" and \"doses\" need to be of same length")
    S <- sigma^2 * diag(1/n)
  }
  
  ## calculate parameters, gradients and results for the asymptotic approximation
  if(missing(bnds)) {
    if(any(!is.element(model, c("linear", "linlog", "quadratic")))){
      message("Message: Need bounds in \"bnds\" for nonlinear models, using default bounds from \"defBnds\".")
      bnds <- defBnds(max(doses))
    }
  }
  nams <- colnames(muMat)
  covMat <- list()
  approx <- matrix(nrow = ncol(muMat), ncol = 3)
  maxdose <- apply(abs(muMat-muMat[1,]), 2, function(x) doses[which.max(x)])
  EDs <- ED(altModels, p)
  EDsUB <- ED(altModels, pUB)
  EDsLB <- ED(altModels, pLB)

  if(!asyApprox & !simulation)
    stop("Need to select either \"asyApprox = TRUE\" or \"simulation = TRUE\"")
  
  if(asyApprox){
    npar <- switch(model,
                   linInt = length(doses),
                   nPars(model))
    bestPar <- matrix(nrow = ncol(muMat), ncol = npar) ## best fit by model to models in altModels
    for(i in 1:ncol(muMat)){
      ## if other model-class approximate best fit
      nam <- gsub("[0-9]", "", nams[i]) # model name (number removed)
      if(nam == model){
        pars <- attr(muMat, "parList")[[i]]
        if(is.element(model, c("betaMod", "linlog")))
          bestPar[i,] <- pars[-length(pars)]
        else
          bestPar[i,] <- pars
        bias <- 0
      } else { ## find the best fit 
        fit <- fitMod(doses, muMat[,i], model=model, S=S,
                      bnds = bnds[[model]], type="general")
        bias <- predict(fit, predType = "effect-curve" , doseSeq = doses[-1])-(muMat[-1,i]-muMat[1,i])
        bestPar[i,] <- coef(fit)
      }
      ## now calculate approximate covariance matrix
      covMat[[i]] <- aprCov(doses, model, bestPar[i,], S, off, scal)
      if(!is.matrix(covMat[[i]])){
        approx[i,] <- NA
      } else {
        ## root-mse
        paVar <- getPredVar(model, bestPar[i,], covMat[[i]], 
                            pDose=doses[-1], scal=scal, off=off)
        approx[i,1] <- sqrt(mean(paVar+bias^2))
        ## Pr(eff_maxdose > 0)
        ind <- which(doses[-1] == maxdose[i])
        paVar <- paVar[ind]
        call <- c(list(c(0,maxdose[i])), as.list(c(bestPar[i,], scal, off)))
        pa <- abs(diff(do.call(model, call)))
        LBmn <- qnorm(alpha, pa, sqrt(paVar))
        approx[i,2] <- pnorm(tau, LBmn, sqrt(paVar), lower.tail = FALSE)
        ## Pr(eff_ED50)
        edvar <- getEDVar(model, bestPar[i,], covMat[[i]], "unrestricted", p,
                          maxdose[i], off=off, scal=scal)
        ed <- calcED(model, bestPar[i,], p, maxdose[i], "continuous",
                                   off=off, scal=scal)
        edsd <- sqrt(edvar)
        approx[i,3] <- pnorm(EDsUB[i], ed, edsd) - pnorm(EDsLB[i], ed, edsd)
      }
    }
    colnames(approx) <- c("dRMSE", "Pow(maxDose)", "P(EDp)")
    rownames(approx) <-   rownames(bestPar) <- nams
    colnames(bestPar) <- rownames(covMat[[1]])
    attr(approx, "bestPar") <- bestPar
    attr(approx, "covMat") <- covMat
  }
  
  if(simulation){
    cat("Running simulations\n")
    requireNamespace("parallel", quietly = TRUE)
    sim <- parallel::mclapply(1:ncol(muMat), function(i){
      if(showSimProgress){
        if(cores == 1){
          cat(sprintf("Scenario %d/%d\n", i, ncol(muMat)))
          pb <- txtProgressBar(style=3, char="*")
        } else {
          cat(sprintf("Scenario %d/%d started\n", i, ncol(muMat)))
        }
      }
      dat <- mvtnorm::rmvnorm(nSim, mean = muMat[,i], sigma = S)
      sims <- numeric(3)
      mse <- LBmn <- edpred <- resp <- numeric(nSim)
      coefs <- vector("list", length = nSim)
      modelSel <- character(nSim)
      for(j in 1:nSim){
        if(showSimProgress & cores == 1)
          setTxtProgressBar(pb, j/nSim)
        fit <- vector("list", length = length(model))
        k <- 1
        for(namMod in model){
          fit[[k]] <- fitMod(dose=doses, dat[j,], model=namMod,
                        S=S, type="general", bnds=bnds[[namMod]])
          k <- k+1
          ## ## this would be faster
          ## fit <- fitMod.raw(doses, dat[j,], model=model,
          ##                                 off=off, scal=scal, nodes=NULL,
          ##                                 S=S, type="general", bnds=bnds,
          ##                                 covarsUsed = FALSE, df = Inf,
          ##                                 control = NULL, 
          ##                                 doseNam = "dose", respNam = "resp")
        }
        aics <- sapply(fit, gAIC)
        fit <- fit[[which.min(aics)]]
        coefs[[j]] <- coef(fit)
        modelSel[j] <- attr(fit, "model")
        
        ## root-MSE of plac-adj dr at doses
        respDoses <- predict(fit, predType = "effect-curve", doseSeq = doses[-1])
        call <- c(list(doses), as.list(c(coef(fit), scal, off)))
        trm <- muMat[-1,i] - muMat[1,i]
        mse[j] <- mean((respDoses-trm)^2)
        ## Pr(LB_maxdose > tau) > 1-alpha
        respMaxD <- predict(fit, predType = "effect-curve", doseSeq = maxdose[i], se.fit=TRUE)
        if(is.na(respMaxD$se.fit)){
          LBmn[j] <- NA
        } else {
          LBmn[j] <- qnorm(alpha, abs(respMaxD$fit), respMaxD$se.fit)
        }
        resp[j] <- respMaxD$fit
        ## ED estimation
        edpred[j] <- ED(fit, p=p)
      }
      ind <- is.na(LBmn)
      NAind <- sum(ind)
      LBmn[ind] <- qnorm(alpha, abs(resp[ind]), sd(resp, na.rm=TRUE))
      sims[1] <- sqrt(mean(mse))
      sims[2] <- mean(LBmn > tau)
      sims[3] <- mean(edpred > EDsLB[i] & edpred < EDsUB[i])
      attr(sims, "NAind") <- NAind
      attr(sims, "coefs") <- coefs
      attr(sims, "model") <- modelSel 
      if(showSimProgress){
        if(cores == 1){
          close(pb)
        } else {
          cat(sprintf("Scenario %d/%d finished\n", i, ncol(muMat)))
        }
      }
      sims
    }, mc.cores=cores)
    NAind <- sapply(sim, function(x) attr(x, "NAind"))
    coefs <- lapply(sim, function(x) attr(x, "coefs"))
    modelSel <- sapply(sim, function(x) attr(x, "model"))
    names(NAind) <- colnames(modelSel) <- names(coefs) <- nams
    rownames(modelSel) <- 1:nSim
    sim <- do.call("rbind", sim)
    colnames(sim) <- c("dRMSE", "Pow(maxDose)", "P(EDp)")
    rownames(sim) <- nams
    attr(sim, "NAind") <- NAind
    attr(sim, "coefs") <- coefs
    attr(sim, "modelSel") <- modelSel
  }

  out <- list(approx = NULL, sim = NULL)
  if(asyApprox)
    out$approx <- approx
  if(simulation){
    out$sim <- sim
    attr(out$sim, "nSim") <- nSim
  }
  attr(out, "model") <- model
  attr(out, "altModels") <- altModels
  attr(out, "doses") <- doses
  attr(out, "off") <- off
  attr(out, "scal") <- scal
  attr(out, "S") <- S
  class(out) <- "planMod"
  out
}

#' @export
print.planMod <- function(x, digits = 3,...){
  model <- attr(x, "model")
  multiMod <- length(model) > 1
  str <- ifelse(multiMod, "s", "")
  cat(sprintf("Fitted Model%s: %s\n\n", str, paste(model, collapse=" ")))
  if(!is.null(x$approx)){
    attr(x$approx, "bestPar") <- NULL
    attr(x$approx, "NAind") <- NULL
    attr(x$approx, "covMat") <- NULL
    cat("Asymptotic Approximations\n")
    print(signif(x$approx, digits))
    cat("\n")
  }
  if(!is.null(x$sim)){
    pltsim <- x$sim
    attr(pltsim, "NAind") <- NULL
    attr(pltsim, "coefs") <- NULL
    attr(pltsim, "modelSel") <- NULL
    attr(pltsim, "nSim") <- NULL
    cat(sprintf("Simulation Results (nSim = %i)\n", attr(x$sim, "nSim")))
    print(signif(pltsim, digits))
    if(multiMod){
      cat("\nSelected models\n")
      res <- apply(attr(x$sim, "modelSel"), 2, tableMatch, match = model)
      print(signif(t(res)/colSums(res), digits))
    }
  }
}

#' Summarize performance metrics for dose-response models
#'
#' @param object,digits object: A planMod object. digits: Digits in summary output
#' @param len Number of equally spaced points to determine the mean-squared error on a grid (cRMSE).
#' @param Delta Additional arguments determining what dose estimate to plot, when \samp{type = "ED"} or \samp{type =
#'   "TD"}
#' @param dLB,dUB Which quantiles to use for calculation of `lengthTDCI` and `lengthEDpCI`. By default dLB =
#'   0.05 and dUB = 0.95, so that this corresponds to a 90% interval.
#' @param ...  Additional arguments (currently ignored)
#'
#' @rdname planMod
#' @method summary planMod
#' @export
summary.planMod <- function(object, digits = 3, len = 101,
                            Delta=NULL, 
                            p=NULL, dLB = 0.05, dUB = 0.95, ...){
  class(object) <- "summary.planMod"
  print(object, digits, len, Delta, p, dLB, dUB, ...)
}

#' @export
print.summary.planMod <- function(x, digits = 3, len = 101,
                                  Delta=NULL, 
                                  p=NULL, dLB = 0.05, dUB = 0.95, ...){
  ## provide more information than print method
  modelSel <- attr(x$sim, "modelSel") 
  model <- attr(x, "model")
  coefs <- attr(x$sim, "coefs")
  altModels <- attr(x, "altModels")
  direction <- attr(altModels, "direction")
  doses <- attr(x, "doses")
  S <- attr(x, "S")
  off <- attr(x, "off")
  scal <- attr(x, "scal")
  ## calculate mean response at doses
  doseSeq <- seq(min(doses), max(doses), length=len)
  muMat <- getResp(altModels, doseSeq)

  if(is.null(x$sim))
    stop("Additional metrics only available if simulations were performed")
  ## calculate average mse of placebo-adjusted dose-response for ANOVA
  CM <- cbind(-1, diag(length(doses)-1))
  mseANOVA <- mean(diag(CM%*%S%*%t(CM)))
  ## calculate predictions
  predList <- getSimEst(x, "dose-response", doseSeq=doseSeq)
  
  out <- matrix(ncol = 5, nrow = ncol(muMat))
  colnames(out) <- c("Eff-vs-ANOVA", "cRMSE", "lengthTDCI", "P(no TD)", "lengthEDCI")
  rownames(out) <- colnames(muMat)
  if(!is.null(Delta)){
    tds <- getSimEst(x, "TD", Delta=Delta, direction=direction)
  }
  if(!is.null(p)){
    eds <- getSimEst(x, "ED", p=p)
  }
  for(i in 1:ncol(muMat)){
    out[i,1] <- mseANOVA/x$sim[i,1]^2
    ## calculate mse of estimating the plac-adj dose-response at fine grid
    ## first calculate placebo-adjusted predictions
    pred <- predList[[i]]
    pred <- (pred-pred[,1])[,-1]
    ## placebo-adjusted response
    mn <- (muMat[-1,i]-muMat[1,i])
    clmn <- colMeans(sweep(pred, 2, mn)^2)
    out[i,2] <- sqrt(mean(clmn))
    ## calculate length of CI for TD
    if(!is.null(Delta)){
      out[i,3] <- diff(quantile(tds[[i]], c(dLB, dUB), na.rm = TRUE))
      out[i,4] <- mean(is.na(tds[[i]]))
    } else {
      out[i,3] <- out[i,4] <- NA
    }
    ## calculate length of CI for ED
    if(!is.null(p)){
      out[i,5] <- diff(quantile(eds[[i]], c(dLB, dUB)))
    } else {
      out[i,5] <- NA
    }
  }
  cat(sprintf("Additional simulation metrics (nSim=%i)\n",
              attr(x$sim, "nSim")))
  print(signif(out, digits=digits))
}


#' Plot to summarize dose-response and dose estimations
#'
#' @param x An object of class planMod
#' @param type Type of plot to produce
#' @param placAdj When \samp{type = "dose-response"}, this determines whether dose-response estimates are shown on
#'   placebo-adjusted or original scale
#' @param xlab,ylab Labels for the plot (ylab only applies for \samp{type = "dose-response"})
#'
#' @rdname planMod
#' @method plot planMod
#' @export
plot.planMod <- function(x, type = c("dose-response", "ED", "TD"),
                         p, Delta, placAdj = FALSE,
                         xlab = "Dose", ylab = "", ...){
  type <- match.arg(type)
  if(type == "dose-response"){
    plotDRSims(x, placAdj = placAdj, xlab=xlab, ylab = ylab)
  } else {
    plotDoseSims(x, type=type, p=p, Delta=Delta, xlab = xlab)
  }
}
