## functions related to creating, plotting candidate model sets

#' Define dose-response models
#'
#' The Mods functions allows to define a set of dose-response models.  The function is used as input object for a number
#' of other different functions.
#'
#' The dose-response models used in this package (see [drmodels()] for details) are of form
#'
#' \deqn{f(d) = \theta_0+\theta_1 f^0(d,\theta_2)}{f(d) = theta0+theta1
#' f0(d,theta2)}
#'
#' where the parameter \eqn{\theta_2}{theta2} is the only non-linear parameter and can be one- or two-dimensional,
#' depending on the used model.
#'
#' One needs to hand over the effect at placebo and the maximum effect in the dose range, from which
#' \eqn{\theta_0,\theta_1}{theta0,theta1} are then back-calculated, the output object is of class \samp{"Mods"}. This
#' object can form the input for other functions to extract the mean response (\samp{getResp}) or target doses
#' ([TD()] and [ED()]) corresponding to the models. It is also needed as input to the functions
#' [powMCT()], [optDesign()]
#'
#' Some models, for example the beta model (\samp{scal}) and the linlog model (\samp{off}) have parameters that are not
#' estimated from the data, they need to be specified via the \samp{addArgs} argument.
#'
#' The default plot method for \samp{Mods} objects is based on a plot using the \samp{lattice} package for backward
#' compatibility. The function \samp{plotMods} function implements a plot using the \samp{ggplot2} package.
#'
#' NOTE: If a decreasing effect is beneficial for the considered response
#' variable it needs to specified here, either by using \samp{direction =
#' "decreasing"} or by specifying a negative "maxEff" argument.
#'
#'
#' @aliases Mods getResp plot.Mods plotMods
#' @param ...  In function Mods:\cr Dose-response model names with parameter values specifying the guesstimates for the
#'   \eqn{\theta_2}{theta2} parameters. See [drmodels()] for a complete list of dose-response models
#'   implemented. See below for an example specification.\cr \cr In function plot.Mods:\cr Additional arguments to the
#'   \samp{xyplot} call.
#' @param doses Dose levels to be used, this needs to include placebo.
#' @param addArgs List containing two entries named "scal" and "off" for the "betaMod" and "linlog". When addArgs is
#'   NULL the following defaults are used \samp{list(scal = 1.2*max(doses), off = 0.01*max(doses), nodes = doses)}.
#' @param fullMod Logical determining, whether the model parameters specified in the Mods function (via the ...
#'   argument) should be interpreted as standardized or the full model parameters.
#' @param placEff,maxEff Specify used placebo effect and the maximum effect over placebo.  Either a numeric vector of
#'   the same size as the number of candidate models or of length one.\cr When these parameters are not specified
#'   \samp{placEff = 0} is assumed, for \samp{maxEff = 1} is assumed, if \samp{direction = "increasing"} and
#'   \samp{maxEff = -1} is assumed, for \samp{direction = "decreasing"}.
#' @param direction Character determining whether the beneficial direction is \samp{increasing} or \samp{decreasing}
#'   with increasing dose levels. This argument is ignored if \samp{maxEff} is specified.
#' @return Returns an object of class \samp{"Mods"}. The object contains the specified model parameter values and the
#'   derived linear parameters (based on \samp{"placEff"} and \samp{"maxEff"}) in a list.
#' @author Bjoern Bornkamp
#' @seealso [Mods()], [drmodels()], [optDesign()], [powMCT()]
#' @references Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and analysis of dose finding studies
#'   combining multiple comparisons and modeling procedures, *Journal of Biopharmaceutical Statistics*, **16**,
#'   639--656
#' @examples
#'
#' ## Example on how to specify candidate models
#'
#' ## Suppose one would like to use the following models with the specified
#' ## guesstimates for theta2, in a situation where the doses to be used are
#' ## 0, 0.05, 0.2, 0.6, 1
#'
#' ## Model            guesstimate(s) for theta2 parameter(s) (name)
#' ## linear           -
#' ## linear in log    -
#' ## Emax             0.05 (ED50)
#' ## Emax             0.3 (ED50)
#' ## exponential      0.7 (delta)
#' ## quadratic       -0.85 (delta)
#' ## logistic         0.4  0.09 (ED50, delta)
#' ## logistic         0.3  0.1 (ED50, delta)
#' ## betaMod          0.3  1.3 (delta1, delta2)
#' ## sigmoid Emax     0.5  2 (ED50, h)
#' ## linInt           0.5 0.75 1 1 (perc of max-effect at doses)
#' ## linInt           0.5 1 0.7 0.5 (perc of max-effect at doses)
#'
#' ## for the linInt model one specifies the effect over placebo for
#' ## each active dose.
#' ## The fixed "scal" parameter of the betaMod is set to 1.2
#' ## The fixed "off"  parameter of the linlog is set to 0.1
#' ## These (standardized) candidate models can be specified as follows
#'
#' models <- Mods(linear = NULL, linlog = NULL, emax = c(0.05, 0.3),
#'                exponential = 0.7, quadratic = -0.85,
#'                logistic = rbind(c(0.4, 0.09), c(0.3, 0.1)),
#'                betaMod = c(0.3, 1.3), sigEmax = c(0.5, 2),
#'                linInt = rbind(c(0.5, 0.75, 1, 1), c(0.5, 1, 0.7, 0.5)),
#'                doses = c(0, 0.05, 0.2, 0.6, 1),
#'                addArgs = list(scal=1.2, off=0.1))
#' ## "models" now contains the candidate model set, as placEff, maxEff and
#' ## direction were not specified a placebo effect of 0 and an effect of 1
#' ## is assumed
#'
#' ## display of specified candidate set using default plot (based on lattice)
#' plot(models)
#' ## display using ggplot2
#' plotMods(models)
#'
#' ## example for creating a candidate set with decreasing response
#' doses <- c(0, 10, 25, 50, 100, 150)
#' fmodels <- Mods(linear = NULL, emax = 25,
#'                    logistic = c(50, 10.88111), exponential = 85,
#'                    betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)),
#'                    linInt = rbind(c(0, 1, 1, 1, 1),
#'                                   c(0, 0, 1, 1, 0.8)),
#'                    doses=doses, placEff = 0.5, maxEff = -0.4,
#'                    addArgs=list(scal=200))
#' plot(fmodels)
#' plotMods(fmodels)
#' ## some customizations (different model names, symbols, line-width)
#' plot(fmodels, lwd = 3, pch = 3, cex=1.2, col="red",
#'      modNams = paste("mod", 1:8, sep="-"))
#'
#' ## for a full-model object one can calculate the responses
#' ## in a matrix
#' getResp(fmodels, doses=c(0, 20, 100, 150))
#'
#' ## calculate doses giving an improvement of 0.3 over placebo
#' TD(fmodels, Delta=0.3, direction = "decreasing")
#' ## discrete version
#' TD(fmodels, Delta=0.3, TDtype = "discrete", doses=doses, direction = "decreasing")
#' ## doses giving 50% of the maximum effect
#' ED(fmodels, p=0.5)
#' ED(fmodels, p=0.5, EDtype = "discrete", doses=doses)
#'
#' plot(fmodels, plotTD = TRUE, Delta = 0.3)
#'
#' ## example for specifying all model parameters (fullMod=TRUE)
#' fmods <- Mods(emax = c(0, 1, 0.1), linear = cbind(c(-0.4,0), c(0.2,0.1)),
#'               sigEmax = c(0, 1.1, 0.5, 3),
#'               doses = 0:4, fullMod = TRUE)
#' getResp(fmods, doses=seq(0,4,length=11))
#' ## calculate doses giving an improvement of 0.3 over placebo
#' TD(fmods, Delta=0.3)
#' ## discrete version
#' TD(fmods, Delta=0.3, TDtype = "discrete", doses=0:4)
#' ## doses giving 50% of the maximum effect
#' ED(fmods, p=0.5)
#' ED(fmods, p=0.5, EDtype = "discrete", doses=0:4)
#' plot(fmods)
#' @export
Mods <- function(..., doses, placEff = 0, maxEff, direction = c("increasing", "decreasing"),
                 addArgs = NULL, fullMod = FALSE){
  if(missing(doses))
    stop("Need to specify dose levels")
  doses <- sort(doses)
  if(doses[1] < -.Machine$double.eps ^ 0.5)
    stop("Only dose-levels >= 0 allowed")
  if(abs(doses[1]) > .Machine$double.eps ^ 0.5)
    stop("Need to include placebo dose")
  ## check for adequate addArgs
  lst <- getAddArgs(addArgs, doses)
  if(lst$scal < max(doses))
    stop("\"scal\" parameter needs to be >= max(doses)")
  if(lst$scal < 0)
    stop("\"scal\" parameter needs to be positive")    
  if(lst$off < 0)
    stop("\"off\" parameter needs to be positive")    
  ## obtain model list
  modL <- list(...)
  nams <- names(modL)
  ## perform some simple check for a valid standModel list
  if(length(nams) != length(unique(nams)))
    stop("only one list entry allowed for each model class")
  checkEntries(modL, doses, fullMod)
  if(!fullMod){ ## assume standardized models
    direction <- match.arg(direction)
    if (missing(maxEff)) 
      maxEff <- ifelse(direction == "increasing", 1, -1)
    modL <- fullMod(modL, doses, placEff, maxEff, lst$scal, lst$off)
  } else {
    ## calculate placEff and maxEff from model pars. For unimodal
    ## models maxEff determination might fail if the dose with maximum
    ## efficacy is not among those used!
    resp <- calcResp(modL, doses, lst$off, lst$scal, lst$nodes)
    placEff <- resp[1,]
    maxEff <- apply(resp, 2, function(x){
      difs <- x-x[1]
      indMax <- which.max(difs)
      indMin <- which.min(difs)
      if(difs[indMax] > 0)
        return(difs[indMax])
      if(difs[indMin] < 0)
        return(difs[indMin])
    })
  }
  attr(modL, "placEff") <- placEff
  attr(modL, "maxEff") <- maxEff
  direc <- unique(ifelse(maxEff > 0, "increasing", "decreasing"))
  if(length(direc) > 1)
    stop("Inconsistent direction of effect specified in maxEff")
  attr(modL, "direction") <- direc
  class(modL) <- "Mods"
  attr(modL, "doses") <- doses
  attr(modL, "scal") <- lst$scal
  attr(modL, "off") <- lst$off
  return(modL)
}


#' Extract mean response from set of dose-response models
#'
#' @param fmodels An object of class Mods
#'
#' @rdname Mods
#' @export
getResp <- function(fmodels, doses){
  ## convenience function for getting the mean responses of
  ## the models in a Mods object (output in matrix)
  if(!inherits(fmodels, "Mods"))
    stop("\"fmodels\" needs to be of class Mods")
  if(missing(doses))
    doses <- attr(fmodels, "doses")
  off <- attr(fmodels, "off")
  scal <- attr(fmodels, "scal")
  nodes <- attr(fmodels, "doses")
  calcResp(fmodels, doses, off=off, scal=scal, nodes=nodes)
}

#' Plot dose-response models using ggplot
#'
#' @param ModsObj For function \samp{plotMods} the \samp{ModsObj} should contain an object of class \samp{Mods}.
#' @param trafo For function \samp{plotMods} there is the option to plot the candidate model set on a transformed scale
#'   (e.g. probability scale if the candidate models are formulated on log-odds scale). The default for \samp{trafo} is
#'   the identity function.
#' @param nPoints Number of points for plotting
#' @param superpose Logical determining, whether model plots should be superposed
#' @param xlab,ylab Label for y-axis and x-axis.
#' @param modNams When \samp{modNams == NULL}, the names for the panels are determined by the underlying model
#'   functions, otherwise the contents of \samp{modNams} are used.
#'
#' @rdname Mods
#' @export
plotMods <- function(ModsObj, nPoints = 200, superpose = FALSE,
                     xlab = "Dose", ylab = "Model means",
                     plotED = FALSE, p = NULL, pLB = NULL, pUB = NULL,
                     modNams = NULL, trafo = function(x) x){
  ## candidate model plot using ggplot2
  ## check for class Mods
  if(!inherits(ModsObj, "Mods"))
    stop("\"ModsObj\" needs to be of class Mods")
  doses <- nodes <- attr(ModsObj, "doses")
  placEff <- attr(ModsObj, "placEff")
  maxEff <- attr(ModsObj, "maxEff")
  off <- attr(ModsObj, "off")
  scal <- attr(ModsObj, "scal")
  nM <- modCount(ModsObj, fullMod = TRUE)
  if(nM > 50)
    stop("too many models in Mods object to plot (> 50 models).")
  
  doseSeq <- sort(union(seq(min(doses), max(doses), length = nPoints), 
                        doses))
  resp <- calcResp(ModsObj, doseSeq, off, scal, nodes)
  resp <- trafo(resp)
  
  if(is.null(modNams)){ # use default model names
    parList <- attr(resp, "parList")
    mod_nams <- getModNams(parList)
  } else { # use specified model names
    if(length(modNams) != nM)
      stop("specified model-names in \"modNams\" of invalid length")
    mod_nams <- modNams
  }
  
  modelfact <- factor(rep(mod_nams, each = length(doseSeq)),
                      levels = mod_nams)
  respdata <- data.frame(response = c(resp), 
                         dose = rep(doseSeq, ncol(resp)),
                         model = modelfact)
  if(superpose){
    pp <- ggplot2::ggplot(respdata, ggplot2::aes(x=.data$dose, y=.data$response, col=.data$model))+
      ggplot2::geom_line(linewidth=1.2)+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
  } else {
    pp <- ggplot2::ggplot(respdata, ggplot2::aes(x=.data$dose, y=.data$response))+
      ggplot2::geom_line(linewidth=1.2)+
      ggplot2::theme_bw()+
      ggplot2::facet_wrap(~model, labeller = ggplot2::label_wrap_gen())
  }
  resp2 <- calcResp(ModsObj, doses, off, scal, nodes)
  resp2 <- trafo(resp2)
  modelfact2 <- factor(rep(mod_nams, each = length(doses)),
                       levels = mod_nams)
  respdata2 <- data.frame(response = c(resp2), 
                          dose = rep(doses, ncol(resp)),
                          model = modelfact2)
  pp <- pp +
    ggplot2::geom_point(ggplot2::aes(x=.data$dose, y=.data$response), size=1.8, data=respdata2) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  # plot also the EDps
  if (plotED) {
    if (is.null(p)) stop("'p' must be provided if plotED is TRUE.")
    
    dfED <- function(p, color) {
      EDp <- as.numeric(ED(ModsObj, p))
      names(EDp) <- mod_nams
      EDpResp <- as.numeric(diag(getResp(ModsObj, doses = EDp)))
      data.frame(model = factor(mod_nams, levels = mod_nams), EDp = EDp, y = EDpResp,
                 y0 = rep(placEff, nM), color = color, p = p)
    }
    plotED <- function(pp, dt) {
      pp +
        ggplot2::geom_segment(data = dt, aes(x = 0, y = y, xend = EDp, yend = y), color = dt$color[1], 
                              linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE) +
        ggplot2::geom_segment(data = dt, ggplot2::aes(x = EDp, y = y0, xend = EDp, yend = y), color = dt$color[1], 
                              linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE) +
        ggplot2::geom_point(data = dt, ggplot2::aes(x = EDp, y = y), color = dt$color[1], size = 1, inherit.aes = FALSE) + 
        ggplot2::geom_text(data = dt, ggplot2::aes(x = EDp, y = y0, label = paste0("ED", round(p*100))), 
                           color = dt$color[1],  vjust = 1.3,  fontface = "bold", size = 3, inherit.aes = FALSE)
    }
    
    dt <- dfED(p, "#E69F00")
    if (!is.null(pLB)) {
      dtLB <- dfED(pLB, "#0072B2")
      pp <- plotED(pp, dtLB)
    }
    if (!is.null(pUB)) {
      dtUB <- dfED(pUB, "#0072B2")
      pp <- plotED(pp, dtUB)
    }
    
    pp <- plotED(pp, dt)
  }
  pp
}

#' Plot dose-response models
#'
#' @param Delta Delta: The target effect size use for the target dose (TD)
#' (Delta should be > 0).
#' @param x Object of class Mods with type Mods
#' @param plotTD \samp{plotTD} is a logical determining, whether the TD should
#' be plotted. \samp{Delta} is the target effect to estimate for the TD.
#' 
#' @rdname Mods
#' @method plot Mods
#' @export
plot.Mods <- function(x, nPoints = 200, superpose = FALSE, xlab = "Dose",
                      ylab = "Model means", modNams = NULL, plotTD = FALSE, Delta, ...){
  plotModels(x, nPoints = nPoints, superpose = superpose, xlab = xlab,
             ylab = ylab, modNams = modNams, plotTD = plotTD, Delta, ...)
}


#' Calculate dose estimates for a fitted dose-response model (via [fitMod()], [bFitMod()]) 
#'  or [maFitMod()]) or a [Mods()] object
#'
#' @description The TD (target dose) is defined as the dose that achieves a target effect of Delta over placebo (if
#' there are multiple such doses, the smallest is chosen):
#'
#' \deqn{TD_\Delta = \min \{x|f(x) > f(0)+\Delta\}}{TD = min {x|f(x) > f(0)+Delta}}
#'
#' If a decreasing trend is beneficial the definition of the TD is
#'
#' \deqn{TD_\Delta = \min \{x|f(x) < f(0)-\Delta\}}{TD = min {x|f(x) < f(0)-Delta}}
#'
#' When \eqn{\Delta}{Delta} is the clinical relevance threshold, then the TD is similar to the usual definition of the
#' minimum effective dose (MED).
#'
#' The ED (effective dose) is defined as the dose that achieves a certain percentage p of the full effect size (within
#' the observed dose-range!) over placebo (if there are multiple such doses, the smallest is chosen).
#' \deqn{ED_p=\min\{x|f(x) > f(0) + p(f(dmax)-f(0))}{ EDp=min{x|f(x) > f(0) + p(f(dmax)-f(0))}}
#'
#' Note that this definition of the EDp is different from traditional definition based on the Emax model,
#' where the EDp is defined relative to the *asymptotic* maximum effect (rather than the maximum effect in the observed dose-range).
#'
#' ED or TD calculation for bootstrap model averaging (maFit) objects is based on first calculating the pointwise median dose-response curve estimate. Then calculating the dose estimate based on this curve.
#'
#' @name Target doses
#' @rdname targdose
#' @aliases ED
#' @param object An object of class c(Mods, fullMod), DRMod, bFitMod or maFit
#' @param Delta,p
#' Delta: The target effect size use for the target dose (TD) (Delta should be > 0).
#'
#' p: The percentage of the dose to use for the effective dose.
#' @param TDtype,EDtype character that determines, whether the dose should be treated as a continuous
#' variable when calculating the TD/ED or whether the TD/ED should be calculated based on a grid of doses specified in \samp{doses}
#' @param direction Direction to be used in defining the TD. This depends on whether an increasing
#' or decreasing of the response variable is beneficial. In case of ED calculation only needed for maFit objects.
#' @param doses Dose levels to be used if \samp{TDtype} or \samp{EDtype} are
#' equal to \samp{"discrete"}. Needs to include placebo, and may not exceed the dose range of the model(s) provided in \samp{object}.
#'
#' @return Returns the dose estimate
#'
#' @author Bjoern Bornkamp
#' @seealso [Mods()], [drmodels()],
#' [fitMod()], [bFitMod()]
#'
#' @examples
#' ## example for creating a "full-model" candidate set placebo response
#' ## and maxEff already fixed in Mods call
#' doses <- c(0, 10, 25, 50, 100, 150)
#' fmodels <- Mods(linear = NULL, emax = 25,
#'                 logistic = c(50, 10.88111), exponential = 85,
#'                 betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)),
#'                 linInt = rbind(c(0, 1, 1, 1, 1),
#'                                c(0, 0, 1, 1, 0.8)),
#'                 doses=doses, placEff = 0, maxEff = 0.4,
#'                 addArgs=list(scal=200))
#' ## calculate doses giving an improvement of 0.3 over placebo
#' TD(fmodels, Delta=0.3)
#' ## discrete version
#' TD(fmodels, Delta=0.3, TDtype = "discrete", doses=doses)
#' ## doses giving 50% of the maximum effect
#' ED(fmodels, p=0.5)
#' ED(fmodels, p=0.5, EDtype = "discrete", doses=doses)
#' plot(fmodels, plotTD = TRUE, Delta = 0.3)
#' @export
TD <- function(object, Delta, TDtype = c("continuous", "discrete"),
               direction = c("increasing", "decreasing"), doses = NULL){
  ## calculate target doses for Mods or DRMod object, return in a numeric
  if(missing(Delta))
    stop("need \"Delta\" to calculate TD")
  if(Delta <= 0)
    stop("\"Delta\" needs to be > 0")
  modNams <- tds <- NULL
  if(inherits(object, "Mods")){
    off <- attr(object, "off")
    scal <- attr(object, "scal")
    nodes <- attr(object, "doses")
    maxD <- max(attr(object, "doses"))
    TDtype <- match.arg(TDtype)
    if(TDtype == "discrete" & any(doses > maxD))
      stop("Doses provided may not exceed the observed dose range")
    ## loop through list
    for(nam in names(object)){
      par <- object[[nam]]
      if(is.matrix(par)){
        for(i in 1:nrow(par)){
          td <- calcTD(nam, par[i,], Delta, TDtype, direction, doses, off, scal, nodes)
          modNams <- c(modNams, paste(nam, i, sep=""))
          tds <- c(tds, td)
        }
      } else { # single model
        td <- calcTD(nam, par, Delta, TDtype, direction, doses, off, scal, nodes)
        modNams <- c(modNams, nam)
        tds <- c(tds, td)
      }
    }
    names(tds) <- modNams
    return(tds)
  }
  if(inherits(object, "DRMod")){ # if fmodel is a DRMod object
    nam <- attr(object, "model")
    par <- sepCoef(object)$DRpars
    scal <- attr(object, "scal")
    off <- attr(object, "off")
    nodes <- attr(object, "nodes")
    doseNam <- attr(object, "doseRespNam")[1]
    maxD <- max(attr(object,"data")[[doseNam]])
    TDtype <- match.arg(TDtype)
    if(TDtype == "discrete" & any(doses > maxD))
      stop("Doses provided may not exceed the observed dose range")
    if(attr(object, "placAdj")){
      par <- c(0, par)
      if(nam == "linInt")
        nodes <- c(0, nodes)
    }
    td <- calcTD(nam, par, Delta, TDtype, direction, doses, off, scal, nodes)
    names(td) <- NULL
    return(td)
  }
  if(inherits(object, "bFitMod")){ # if fmodel is a bFitMod object
    nam <- attr(object, "model")
    scal <- attr(object, "scal")
    off <- attr(object, "off")
    nodes <- attr(object, "nodes")
    doseNam <- attr(object, "doseRespNam")[1]
    maxD <- max(attr(object,"data")[[doseNam]])
    TDtype <- match.arg(TDtype)
    if(TDtype == "discrete" & any(doses > maxD))
      stop("Doses provided may not exceed the observed dose range")
    if(attr(object, "placAdj")){
      if(nam == "linInt")
        nodes <- c(0, nodes)
    }
    td <- apply(object$samples, 1, function(x){
      if(attr(object, "placAdj")){
        par <- c(0, x)
      } else {
        par <- x
      }
      calcTD(nam, par, Delta, TDtype, direction, doses, off, scal, nodes)
    })
    return(td)
  }
  if(inherits(object, "maFit")){
    direction <- match.arg(direction, c("increasing", "decreasing"))
    TDtype <- match.arg(TDtype)
    maxD <-  max(object$args$dose)
    if(TDtype == "discrete"){
      if(is.null(doses))
        stop("For TDtype = \"discrete\" need the possible doses in doses argument")
      if(doses[1] != 0)
        stop("need placebo dose for TD calculation")
      if(any(doses > maxD))
        stop("Doses provided may not exceed the observed dose range")
      doseSeq <- doses
    } else { # TDtype == "continuous"
      doseSeq <- seq(0, maxD, length=501) 
    }
    pred_med <- predict(object, doseSeq = doseSeq, summaryFct = stats::median)
    
    if(direction == "decreasing")
      pred_med <- -pred_med
    
    ind <- which(pred_med > pred_med[1] + Delta)
    
    if (length(ind)>0) {
      return(min(doseSeq[ind]))
    } else {
      return(NA)
    }
  }
}

#' #' Calculate effective dose for a dose-response model
#'
#' @rdname targdose
#' @export
ED <- function(object, p, EDtype = c("continuous", "discrete"),
               direction = c("increasing", "decreasing"), doses = NULL){
  ## calculate target doses for Mods or DRMod object, return in a numeric
  if(missing(p))
    stop("need \"p\" to calculate ED")
  if(p <= 0 | p >= 1)
    stop("\"p\" needs to be in (0,1)")
  modNams <- eds <- NULL
  if(inherits(object, "Mods")){
    off <- attr(object, "off")
    scal <- attr(object, "scal")
    nodes <- attr(object, "doses")
    maxD <- max(attr(object, "doses"))
    ## loop through list
    for(nam in names(object)){
      par <- object[[nam]]
      if(is.matrix(par)){
        for(i in 1:nrow(par)){
          ed <- calcED(nam, par[i,], p, maxD, EDtype, doses, off, scal, nodes)
          modNams <- c(modNams, paste(nam, i, sep=""))
          eds <- c(eds, ed)
        }
      } else { # single model
        ed <- calcED(nam, par, p, maxD, EDtype, doses, off, scal, nodes)
        modNams <- c(modNams, nam)
        eds <- c(eds, ed)
      }
    }
    names(eds) <- modNams
    return(eds)
  }
  if(inherits(object, "DRMod")){ # if fmodel is a DRMod object
    nam <- attr(object, "model")
    par <- sepCoef(object)$DRpars
    doseNam <- attr(object, "doseRespNam")[1]
    maxD <- max(attr(object,"data")[[doseNam]])
    scal <- attr(object, "scal")
    off <- attr(object, "off")
    nodes <- attr(object, "nodes")
    if(attr(object, "placAdj")){
      par <- c(0, par)
      if(nam == "linInt")
        nodes <- c(0, nodes)
    }
    ed <- calcED(nam, par, p, maxD, EDtype, doses, off, scal, nodes)
    names(ed) <- NULL
    return(ed)
  }
  if(inherits(object, "bFitMod")){ # if fmodel is a bFitMod object
    nam <- attr(object, "model")
    scal <- attr(object, "scal")
    off <- attr(object, "off")
    nodes <- attr(object, "nodes")
    if(attr(object, "placAdj")){
      if(nam == "linInt")
        nodes <- c(0, nodes)
    }
    doseNam <- attr(object, "doseRespNam")[1]
    maxD <- max(attr(object,"data")[[doseNam]])
    ed <- apply(object$samples, 1, function(x){
      if(attr(object, "placAdj")){
        par <- c(0, x)
      } else {
        par <- x
      }
      calcED(nam, par, p, maxD, EDtype, doses, off, scal, nodes)
    })
    return(ed)
  }
  if(inherits(object, "maFit")){
    EDtype <- match.arg(EDtype)
    maxD <- max(object$args$dose)
    if(EDtype == "discrete"){
      if(is.null(doses))
        stop("For EDtype = \"discrete\" need the possible doses in doses argument")
      if(!any(doses == 0))
        stop("need placebo dose for ED calculation")
      if(any(doses > maxD))
        stop("Doses provided may not exceed the observed dose range.")
    }
    if(missing(direction)){
      stop("Need to provide direction of dose-response (\"increasing\" or \"decreasing\") for objects of class maFitMod.")
    } else {
      direction <- match.arg(direction, c("increasing", "decreasing"))
    }
    if(EDtype == "discrete"){
      doseSeq <- unique(c(sort(doses), maxD)) 
    } else { # EDtype == "continuous"
      doseSeq <- seq(0, maxD, length=501) 
    }
    pred_med <- predict(object, doseSeq = doseSeq, summaryFct = stats::median)
    
    if(direction == "decreasing")
      pred_med <- -pred_med
    
    difs <- (pred_med - pred_med[1])
    ind <- which(difs > p*max(difs))
    if(length(ind) == 0)
      return(NA)
    
    edose <- min(doseSeq[ind])
    if (EDtype == "continuous" | edose %in% doses) {## don't return maxD if it was not in originally provided doses for discrete type
      return(edose)
    } else {
      return(NA)
    }
  }
}
