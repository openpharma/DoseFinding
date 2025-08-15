## optimal designs for model-fitting

#' Function to calculate optimal designs
#'
#' Given a set of models (with full parameter values and model probabilities) the \samp{optDesign} function calculates
#' the optimal design for estimating the dose-response model parameters (D-optimal) or the design for estimating the
#' target dose (TD-optimal design) \insertCite{@see @dette2008}{DoseFinding}, or a mixture of these two
#' criteria. The design can be plotted (together with the candidate models) using \samp{plot.design}. \samp{calcCrit}
#' calculates the design criterion for a discrete set of design(s). \samp{rndDesign} provides efficient rounding for the
#' calculated continous design to a finite sample size.
#'
#' Let \eqn{M_m}{M_m} denote the Fisher information matrix under model m (up to
#' proportionality). \eqn{M_m}{M_m} is given by \eqn{\sum a_i w_i }{\sum a_i
#' w_i g_i^Tg_i}\eqn{ g_i^Tg_i}{\sum a_i w_i g_i^Tg_i}, where \eqn{a_i}{a_i} is
#' the allocation weight to dose i, \eqn{w_i}{w_i} the weight for dose i specified via \samp{weights} and \eqn{g_i}{g_i}
#' the gradient vector of model m evaluated at dose i.
#'
#' For \samp{designCrit = "Dopt"} the code minimizes the design criterion
#'
#' \deqn{-\sum_{m}p_m/k_m \log(\det(M_m))}{-sum_m p_m/k_m log(det(M_m))} where \eqn{p_m}{p_m} is the probability for
#' model m and \eqn{k_m}{k_m} is the number of parameters for model m. When \samp{standDopt = FALSE} the \eqn{k_m}{k_m}
#' are all assumed to be equal to one.
#'
#' For \samp{designCrit = "TD"} the code minimizes the design criterion
#'
#' \deqn{\sum_{m}p_m \log(v_m)}{sum_m p_m log(v_m)} where \eqn{p_m}{p_m} is the probability for model m and
#' \eqn{v_m}{v_m} is proportional to the asymptotic
#' variance of the TD estimate and given by \eqn{b_m'M_m^{-}b_m}{b_m'Minv_m
#' b_m} \insertCite{@see @dette2008, p. 1227 for details}{DoseFinding}.
#'
#' For \samp{designCrit = "Dopt&TD"} the code minimizes the design criterion
#' \deqn{\sum_{m}p_m(-0.5\log(\det(M_m))/k_m+0.5\log(v_m))}{sum_m
#' p_m(-0.5log(det(M_m))/k_m+0.5log(v_m))}
#'
#' Again, for \samp{standDopt = FALSE} the \eqn{k_m}{k_m} are all assumed to be equal to one.
#'
#' For details on the \samp{rndDesign} function, see \insertCite{pukelsheim1993;textual}{DoseFinding}, Chapter 12.
#'
#' @aliases optDesign plot.DRdesign calcCrit rndDesign
#' @param models An object of class \samp{c(Mods, fullMod)}, see the [Mods()] function for details. When an TD
#'   optimal design should be calculated, the TD needs to exist for all models. If a D-optimal design should be
#'   calculated, you need at least as many doses as there are parameters in the specified models.
#' @param probs Vector of model probabilities for the models specified in \samp{models}, assumed in the same order as
#'   specified in models
#' @param doses Optional argument. If this argument is missing the doses attribute in the \samp{c(Mods, fullMod)} object
#'   specified in \samp{models} is used.
#' @param designCrit Determines which type of design to calculate. "TD&Dopt" uses both optimality criteria with equal
#'   weight.
#' @param Delta Target effect needed for calculating "TD" and "TD&Dopt" type designs.
#' @param standDopt Logical determining, whether the D-optimality criterion (specifically the log-determinant) should be
#'   standardized by the number of parameters in the model or not (only of interest if type = "Dopt" or type =
#'   "TD&Dopt"). This is of interest, when there is more than one model class in the candidate model set (traditionally
#'   standardization this is done in the optimal design literature).
#' @param weights Vector of weights associated with the response at the doses. Needs to be of the same length as the
#'   \samp{doses}.  This can be used to calculate designs for heteroscedastic or for generalized linear model
#'   situations.
#' @param nold,n When calculating an optimal design at an interim analysis, \samp{nold} specifies the vector of sample
#'   sizes already allocated to the different doses, and \samp{n} gives sample size for the next cohort.
#'
#'   For \samp{optimizer = "exact"} one always needs to specify the total sample size via \samp{n}.
#' @param control List containing control parameters passed down to numerical optimization algorithms
#'   ([optim()], [nlminb()] or solnp function).\cr
#'
#'   For \samp{type = "exact"} this should be a list with possible entries \samp{maxvls1} and \samp{maxvls2},
#'   determining the maximum number of designs allowed for passing to the criterion function (default
#'   \samp{maxvls2=1e5}) and for creating the initial unrestricted matrix of designs (default \samp{maxvls1=1e6}). In
#'   addition there can be an entry \samp{groupSize} in case the patients are allocated a minimum group size is
#'   required.
#' @param optimizer Algorithm used for calculating the optimal design. Options "Nelder-Mead" and "nlminb" use the
#'   [optim()] and [nlminb()] function and use a trigonometric transformation to turn the
#'   constrained optimization problem into an unconstrained one \insertCite{@see @atkinson2007 pages 130,131}{DoseFinding}.
#'
#'   Option "solnp" uses the solnp function from the Rsolnp package, which implements an optimizer for non-linear
#'   optimization under general constraints.
#'
#'   Option "exact" tries all given combinations of \samp{n} patients to the given dose groups (subject to the bounds
#'   specified via \samp{lowbnd} and \samp{uppbnd}) and reports the best design. When patients are only allowed to be
#'   allocated in groups of a certain \samp{groupSize}, this can be adjusted via the control argument.
#'   \samp{n/groupSize} and \samp{length(doses)} should be rather small for this approach to be feasible.
#'
#'   When the number of doses is small (<8) usually \samp{"Nelder-Mead"} and \samp{"nlminb"} are best suited
#'   (\samp{"nlminb"} is usually a bit faster but less stable than \samp{"Nelder-Mead"}). For a larger number of doses
#'   \samp{"solnp"} is the most reliable option (but also slowest) (\samp{"Nelder-Mead"} and \samp{"nlminb"} often
#'   fail). When the sample size is small \samp{"exact"} provides the optimal solution rather quickly.
#' @param lowbnd,uppbnd Vectors of the same length as dose vector specifying upper and lower limits for the allocation
#'   weights. This option is only available when using the "solnp" and "exact" optimizers.
#' @param userCrit User defined design criterion, should be a function that given a vector of allocation weights and the
#'   doses returns the criterion function. When specified \samp{models} does not need to be handed over.
#'
#'   The first argument of \samp{userCrit} should be the vector of design weights, while the second argument should be
#'   the \samp{doses} argument (see example below). Additional arguments to \samp{userCrit} can be passed via ...
#' @param ...  For function \samp{optDesign} these are additional arguments passed to \samp{userCrit}.\cr \cr For
#'   function \samp{plot.design} these are additional parameters passed to [plot.Mods()].\cr
#' @note In some cases (particularly when the number of doses is large, e.g. 7 or larger) it might be necessary to allow
#'   a larger number of iterations in the algorithm (via the argument \samp{control}), particularly for the Nelder-Mead
#'   algorithm. Alternatively one can use the solnp optimizer that is usually the most reliable, but not fastest option.
#' @author Bjoern Bornkamp
#' @seealso [Mods()], [drmodels()]
#' @references 
#' \insertRef{atkinson2007}{DoseFinding}
#' 
#' \insertRef{dette2008}{DoseFinding}
#' 
#' \insertRef{pinheiro2017}{DoseFinding}
#'
#' \insertRef{pukelsheim1993}{DoseFinding}
#' @examples
#'
#' ## calculate designs for Emax model
#' doses <- c(0, 10, 100)
#' emodel <- Mods(emax = 15, doses=doses, placEff = 0, maxEff = 1)
#' optDesign(emodel, probs = 1)
#' ## TD-optimal design
#' optDesign(emodel, probs = 1, designCrit = "TD", Delta=0.5)
#' ## 50-50 mixture of Dopt and TD
#' optDesign(emodel, probs = 1, designCrit = "Dopt&TD", Delta=0.5)
#' ## use dose levels different from the ones specified in emodel object
#' des <- optDesign(emodel, probs = 1, doses = c(0, 5, 20, 100))
#' ## plot models overlaid by design
#' plot(des, emodel)
#' ## round des to a sample size of exactly 90 patients
#' rndDesign(des, n=90) ## using the round function would lead to 91 patients
#'
#' ## illustrating different optimizers (see Note above for more comparison)
#' optDesign(emodel, probs=1, optimizer="Nelder-Mead")
#' optDesign(emodel, probs=1, optimizer="nlminb")
#' ## optimizer solnp (the default) can deal with lower and upper bounds:
#' optDesign(emodel, probs=1, designCrit = "TD", Delta=0.5,
#'           optimizer="solnp", lowbnd = rep(0.2,3))
#' ## exact design using enumeration of all possibilites
#' optDesign(emodel, probs=1, optimizer="exact", n = 30)
#' ## also allows to fix minimum groupSize
#' optDesign(emodel, probs=1, designCrit = "TD", Delta=0.5,
#'           optimizer="exact", n = 30,  control = list(groupSize=5))
#'
#'
#' ## optimal design at interim analysis
#' ## assume there are already 10 patients on each dose and there are 30
#' ## left to randomize, this calculates the optimal increment design
#' optDesign(emodel, 1, designCrit = "TD", Delta=0.5,
#'           nold = c(10, 10, 10), n=30)
#'
#' ## use a larger candidate model set
#' doses <- c(0, 10, 25, 50, 100, 150)
#' fmods <- Mods(linear = NULL, emax = 25, exponential = 85,
#'              linlog = NULL, logistic = c(50, 10.8811),
#'              doses = doses, addArgs=list(off=1),
#'              placEff=0, maxEff=0.4)
#' probs <- rep(1/5, 5) # assume uniform prior
#' desDopt <- optDesign(fmods, probs, optimizer = "nlminb")
#' desTD <- optDesign(fmods, probs, designCrit = "TD", Delta = 0.2,
#'                    optimizer = "nlminb")
#' desMix <- optDesign(fmods, probs, designCrit = "Dopt&TD", Delta = 0.2)
#' ## plot design and truth
#' plot(desMix, fmods)
#'
#' ## illustrate calcCrit function
#' ## calculate optimal design for beta model
#' doses <- c(0, 0.49, 25.2, 108.07, 150)
#' models <- Mods(betaMod = c(0.33, 2.31), doses=doses,
#'                 addArgs=list(scal=200),
#'                 placEff=0, maxEff=0.4)
#' probs <- 1
#' deswgts <- optDesign(models, probs, designCrit = "Dopt",
#'                      control=list(maxit=1000))
#' ## now compare this design to equal allocations on
#' ## 0, 10, 25, 50, 100, 150
#' doses2 <- c(0, 10, 25, 50, 100, 150)
#' design2 <- c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6)
#' crit2 <- calcCrit(design2, models, probs, doses2, designCrit = "Dopt")
#' ## ratio of determinants (returned criterion value is on log scale)
#' exp(deswgts$crit-crit2)
#'
#' ## example for calculating an optimal design for logistic regression
#' doses <- c(0, 0.35, 0.5, 0.65, 1)
#' fMod <- Mods(linear = NULL, doses=doses, placEff=-5, maxEff = 10)
#' ## now calculate weights to use in the covariance matrix
#' mu <- as.numeric(getResp(fMod, doses=doses))
#' mu <- 1/(1+exp(-mu))
#' weights <- mu*(1-mu)
#' des <- optDesign(fMod, 1, doses, weights = weights)
#'
#' ## one can also specify a user defined criterion function
#' ## here D-optimality for cubic polynomial
#' CubeCrit <- function(w, doses){
#'   X <- cbind(1, doses, doses^2, doses^3)
#'   CVinv <- crossprod(X*w)
#'   -log(det(CVinv))
#' }
#' optDesign(doses = c(0,0.05,0.2,0.6,1),
#'           designCrit = "userCrit", userCrit = CubeCrit,
#'           optimizer = "nlminb")
#' @export
optDesign <- function(models, probs, doses,
                      designCrit = c("Dopt", "TD", "Dopt&TD", "userCrit"),
                      Delta, standDopt = TRUE, weights,
                      nold = rep(0, length(doses)),  n,
                      control=list(), 
                      optimizer = c("solnp", "Nelder-Mead", "nlminb", "exact"),
                      lowbnd = rep(0, length(doses)), uppbnd = rep(1, length(doses)),
                      userCrit, ...){
  if(!missing(models)){
    if(!inherits(models, "Mods"))
      stop("\"models\" needs to be of class Mods")
    direction <- attr(models, "direction")
    off <- attr(models, "off")
    scal <- attr(models, "scal")
    if(missing(doses))
      doses <- attr(models, "doses")
  } else {
    if(missing(userCrit))
      stop("either \"models\" or \"userCrit\" need to be specified")
    if(missing(doses))
      stop("For userCrit one always needs to specify doses")
  }
  ## check arguments
  designCrit <- match.arg(designCrit)
  optimizer <- match.arg(optimizer)
  if(missing(n)){
    if(optimizer == "exact")
      stop("need to specify sample size via n argument")
    if(any(nold > 0))
      stop("need to specify sample size for next cohort via n argument")
    n <- 1 ## value is arbitrary in this case
  } else {
    if(length(n) > 1)
      stop("n needs to be of length 1")
  }
  if(missing(Delta)){
    if(designCrit %in% c("TD", "Dopt&TD"))
      stop("need to specify target difference \"Delta\"")
  } else {
    if(Delta <= 0)
      stop("\"Delta\" needs to be > 0, if curve decreases use \"direction = decreasing\"")
  }
  if(missing(weights)){
    weights <- rep(1, length(doses))
  } else {
    if(length(weights) != length(doses))
      stop("weights and doses need to be of equal length")
  }
  if(length(lowbnd) != length(doses))
    stop("lowbnd needs to be of same length as doses")
  if(length(uppbnd) != length(doses))
    stop("uppbnd needs to be of same length as doses")
  if(any(lowbnd > 0) | any(uppbnd < 1)){
    if(optimizer != "solnp" & optimizer != "exact")
      stop("only optimizers solnp or exact can handle additional constraints on allocations")
  }
  if(sum(lowbnd) > 1)
    stop("Infeasible lower bound specified (\"sum(lowbnd) > 1\"!)")
  if(sum(uppbnd) < 1)
    stop("Infeasible upper bound specified (\"sum(lowbnd) < 1\"!)")
  if(!is.logical(standDopt))
    stop("standDopt needs to contain a logical value")
  standInt <- as.integer(standDopt) # use standardized or non-stand. D-optimality
  nD <- length(doses)
  if(designCrit == "TD" | designCrit == "Dopt&TD"){ # check whether TD exists in (0,max(dose))
    if(length(unique(direction)) > 1)
      stop("need to provide either \"increasing\" or \"decreasing\" as direction to optDesign, when TD optimal designs should be calculated")
    direction <- unique(direction)
    tdMods <- TD(models, Delta, "continuous", direction)
    tdMods[tdMods > max(doses)] <- NA
    if(any(is.na(tdMods)))
      stop("TD does not exist for ",
           paste(names(tdMods)[is.na(tdMods)], collapse=", " ), " model(s)")
  }
  if(designCrit == "Dopt" | designCrit == "Dopt&TD"){ # check whether Fisher matrix can be singular
    np <- nPars(names(models))
    if(max(np) > length(doses))
      stop("need at least as many dose levels as there are parameters to calculate Dopt design.")
  } 

  ## use transformation for Nelder-Mead and nlminb
  if(is.element(optimizer, c("Nelder-Mead", "nlminb"))){ 
    transform <- transTrig
  } else {
    transform <- idtrans
  }
  
  if(designCrit != "userCrit"){ # prepare criterion function
    ## check arguments
    if(abs(sum(probs)-1) > sqrt(.Machine$double.eps)){
      stop("probs need to sum to 1")
    }
    ## prepare criterion function
    lst <- calcGrads(models, doses, weights,
                     Delta, off, scal, direction, designCrit)
    ## check for invalid values (NA, NaN and +-Inf)
    checkInvalid <- function(x)
      any(is.na(x)|(is.nan(x)|!is.finite(x)))
    grInv <- checkInvalid(lst$modgrads)
    MvInv <- ifelse(designCrit != "Dopt", checkInvalid(lst$TDgrad), FALSE)
    if(grInv | MvInv)
      stop("NA, NaN or +-Inf in gradient or bvec")
    ## prepare arguments before passing to C
    M <- as.integer(length(probs))
    if(M != length(lst$nPar))
      stop("probs of wrong length")
    if(length(lst$modgrads) != length(doses)*sum(lst$nPar))
      stop("Gradient of wrong length.")
    if(length(nold) != nD)
      stop("Either nold or doses of wrong length.")
    nD <- as.integer(nD)
    p <- as.integer(lst$nPar)
    intdesignCrit <- match(designCrit, c("TD", "Dopt", "Dopt&TD"))
    objFunc <- function(par){
      optFunc(par, xvec=as.double(lst$modgrads),
              pvec=as.integer(p), nD=nD, probs=as.double(probs),
              M=M, n=as.double(n), nold = as.double(nold),
              bvec=as.double(lst$TDgrad), trans = transform,
              standInt = standInt,designCrit = as.integer(intdesignCrit))
    }
  } else { # user criterion
    if(missing(userCrit))
      stop("need design criterion in userCrit when specified")
    if(!is.function(userCrit))
      stop("userCrit needs to be a function")
    objFunc <- function(par){
      par2 <- do.call("transform", list(par, nD))
      userCrit((par2*n+nold)/(sum(nold)+n), doses, ...)
    }
  }

  ## perform actual optimization
  if(optimizer != "exact"){ # use callOptim function
    res <- callOptim(objFunc, optimizer, nD, control, lowbnd, uppbnd)
    if(optimizer == "Nelder-Mead" | optimizer == "nlminb"){ # transform results back
      des <- transTrig(res$par, length(doses))
      if(optimizer == "Nelder-Mead"){
        crit <- res$value
      } else {
        crit <- res$objective
      }
    }
    if(optimizer == "solnp"){ # no need to transform back
      des <- res$pars
      crit <- res$values[length(res$values)]
    }
    if(res$convergence){
      message("Message: algorithm indicates no convergence, the 'optimizerResults'
               attribute of the returned object contains more details.")
    }
  } else { # exact criterion (enumeration of all designs)
    ## enumerate possible exact designs
    con <- list(maxvls1 = 1e6, maxvls2 = 1e5, groupSize = 1)
    con[(namc <- names(control))] <- control    
    mat <- getDesMat(n, nD, lowbnd, uppbnd,
                     con$groupSize, con$maxvls1, con$maxvls2)
    designmat <- sweep(mat*n, 2, nold, "+")
    res <- sweep(designmat, 2, n+sum(nold), "/")
    ## evaluate criterion function
    if(designCrit != "userCrit"){
      critv <- calcCrit(res, models, probs, doses,
                        designCrit, Delta, standDopt,
                        weights, nold, n)
    } else {
      critv <- apply(res, 1, objFunc)
    }
    des <- mat[which.min(critv),]
    crit <- min(critv)
  }
  out <- list(crit = crit, design = des, doses = doses, n = n,
              nold = nold, designCrit = designCrit)
  attr(out, "optimizerResults") <- res
  class(out) <- "DRdesign"
  out
}

#' Calculate design criteria for set of designs
#'
#' @param design Argument for \samp{rndDesign} and \samp{calcCrit} functions: Numeric vector (or matrix) of allocation
#'   weights for the different doses. The rows of the matrices need to sum to 1. Alternatively also an object of class
#'   "DRdesign" can be used for \samp{rndDesign}. Note that there should be at least as many design points available as
#'   there are parameters in the dose-response models selected in `models` (otherwise the code returns an NA).
#'
#' @rdname optDesign
#' @export
calcCrit <- function(design, models, probs, doses, 
                     designCrit = c("Dopt", "TD", "Dopt&TD"),
                     Delta, standDopt = TRUE, weights,
                     nold = rep(0, length(doses)), n){
  if(!inherits(models, "Mods"))
    stop("\"models\" needs to be of class Mods")
  off <- attr(models, "off")
  scal <- attr(models, "scal")
  if(missing(doses))
    doses <- attr(models, "doses")  
  ## extract design
  if(inherits(design, "DRdesign"))
    design <- design$design
  if(!is.numeric(design))
    stop("design needs to be numeric")
  if(!is.matrix(design))
    design <- matrix(design, ncol = length(design))
  if(ncol(design) != length(doses))
    stop("design and doses should be of the same length")      
  if(any(abs(rowSums(design)-1) > 0.001))
    stop("design needs to sum to 1")
  if(missing(n)){
    n <- 1 # value arbitrary
  } else {
    if(length(n) > 1)
      stop("n needs to be of length 1")
  }
  if(missing(weights)){
    weights <- rep(1, length(doses))
  } else {
    if(length(weights) != length(doses))
      stop("weights and doses need to be of equal length")
  }
  designCrit <- match.arg(designCrit)
  if(missing(Delta) & substr(designCrit, 1, 3) == "TD")
    stop("need to specify clinical relevance parameter")
  direction <- attr(models, "direction")

  if(designCrit == "TD" | designCrit == "Dopt&TD"){ # check whether TD exists in (0,max(dose))
    if(length(unique(direction)) > 1)
      stop("need to provide either \"increasing\" or \"decreasing\" as direction to optDesign, when TD optimal designs should be calculated")
    direction <- unique(direction)
    tdMods <- TD(models, Delta, "continuous", direction)
    tdMods[tdMods > max(doses)] <- NA
    if(any(is.na(tdMods)))
      stop("TD does not exist for ",
           paste(names(tdMods)[is.na(tdMods)], collapse=", " ), " model(s)")
  }
  if(designCrit == "Dopt" | designCrit == "Dopt&TD"){ # check whether Fisher matrix can be singular
    np <- nPars(names(models))
    if(max(np) > length(doses))
      stop("need more dose levels to calculate Dopt design.")
  }
  if(!is.logical(standDopt))
    stop("standDopt needs to contain a logical value")
  standInt <- as.integer(standDopt)
  lst <- calcGrads(models, doses, weights, Delta, off, scal,
                   direction, designCrit)
  ## check for invalid values (NA, NaN and +-Inf)
  checkInvalid <- function(x)
    any(is.na(x)|(is.nan(x)|!is.finite(x)))
  grInv <- checkInvalid(lst$modgrads)
  MvInv <- ifelse(designCrit != "Dopt", checkInvalid(lst$TDgrad), FALSE)
  if(grInv | MvInv)
    stop("NA, NaN or +-Inf in gradient or bvec")
  ## prepare for input into C
  M <- as.integer(length(probs))
  nD <- as.integer(length(doses))
  if(M != length(lst$nPar))
    stop("Probs of wrong length")
  if(length(lst$modgrads) != length(doses)*sum(lst$nPar))
    stop("Gradient of wrong length.")
  
  if(length(nold) != nD)
    stop("Either nold or doses of wrong length.")
  p <- as.integer(lst$nPar)
  intdesignCrit <- match(designCrit, c("TD", "Dopt", "Dopt&TD"))
  res <- numeric(nrow(design))
  ## check for sufficient number of design points
  iter <- 1:nrow(design)
  design0 <- sweep(design, 2, nold, "+")
  count <- apply(design0, 1, function(x) sum(x > 0.0001))
  ind <- count < max(p[probs > 0])
  if(any(ind)){
    iter <- iter[!ind]
    res[ind] <- NA
    if(all(is.na(res)))
      warning("need at least as many dose levels in the design as parameters in the model")
  }
  for(i in iter){
    res[i] <- optFunc(design[i,], xvec=as.double(lst$modgrads),
                      pvec=as.integer(p), nD=nD, probs=as.double(probs),
                      M=M, n=as.double(n), nold = as.double(nold),
                      bvec=as.double(lst$TDgrad), trans = idtrans,
                      standInt = standInt, designCrit = as.integer(intdesignCrit))
  }
  res
}

#' @export
print.DRdesign <- function(x, digits = 5, eps = 0.001, ...){
  nam <- switch(x$designCrit,
                "TD" = "TD",
                "Dopt" = "D",
                "Dopt&TD" = "TD and D mixture",
                "userCrit" = "userCrit")
  cat("Calculated", nam, "- optimal design:\n")
  ind <- x$design > eps
  vec <- x$design[ind]
  names(vec) <- x$doses[ind]
  print(round(vec, digits = digits))
}


#' Efficiently round calculated design to a finite sample size
#'
#' @param eps Argument for \samp{rndDesign} function: Value under which elements of w will be regarded as 0.
#' @rdname optDesign
#' @export
rndDesign <- function(design, n, eps = 0.0001){

  if(missing(n))
    stop("total sample size \"n\" needs to be specified")
  n <- round(n) # ensure n is an integer (at least numerically)
  if(inherits(design, "DRdesign")){
    design <- design$design
  }
  if(!inherits(design, "numeric"))
    stop("design needs to be a numeric vector.")
  zeroind <- design < eps
  if(any(zeroind)){
    design <- design[!zeroind]/sum(design[!zeroind])
  }
  l <- sum(!zeroind)
  nn <- ceiling((n-0.5*l)*design)
  while(sum(nn)!=n){
    if(sum(nn)<n){
      indmin <- which.is.max(-nn/design)
      nn[indmin] <- nn[indmin]+1
    } else {
      indmax <- which.is.max((nn-1)/design)
      nn[indmax] <- nn[indmax]-1
    }
  }
  if(any(zeroind)){
    out <- numeric(length(design))
    out[zeroind] <- 0
    out[!zeroind] <- nn
    return(out)
  } else {
    nn
  }
}


#' Plot optimal designs
#'
#' @param x Object of class \samp{DRdesign} (for \samp{plot.design})
#' @param lwdDes,colDes Line width and color of the lines plotted for the design (in \samp{plot.design})
#'
#' @rdname optDesign
#' @method plot DRdesign
#' @export
plot.DRdesign <- function(x, models, lwdDes = 10, colDes = rgb(0,0,0,0.3), ...){
  if(missing(models))
    stop("need object of class Mods to produce plot")
  plot(models, ...)
  layoutmat <- lattice::trellis.currentLayout()
  nc <- ncol(layoutmat)
  nr <- nrow(layoutmat)
  total <- sum(layoutmat > 0)
  z <- 1
  for(i in 1:nc){
    for(j in 1:nr){
      if(z > total)
        break
      lattice::trellis.focus("panel", i, j)
      args <- lattice::trellis.panelArgs()
      miny <- min(args$y)
      maxy <- max(args$y)
      dy <- maxy-miny
      for(k in 1:length(x$doses)){
        yy <- c(0,(x$design*dy)[k])+miny
        xx <- rep(x$doses[k],2)
        lattice::panel.xyplot(xx, yy, type="l", col = colDes, lwd = lwdDes)
      }
      z <- z+1
      lattice::trellis.unfocus()
    }
  }
}
