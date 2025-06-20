## functions for calculating optimal contrasts and critical value

#' Calculate optimal contrasts
#' 
#' This function calculates a contrast vectors that are optimal for detecting
#' certain alternatives. The contrast is optimal in the sense of maximizing the
#' non-centrality parameter of the underlying contrast test statistic:
#' \deqn{\frac{c'\mu}{\sqrt{c'Sc}}}{c'mu/sqrt(c'Sc).} Here \eqn{\mu}{mu} is the
#' mean vector under the alternative and \eqn{S}{S} the covariance matrix
#' associated with the estimate of \eqn{\mu}{mu}.  The optimal contrast is
#' given by \deqn{c^{opt} \propto S^{-1}\left(\mu - \frac{\mu^{\prime}S^{-1}1}
#' {1^\prime S^{-1} 1}\right),}{c propto S^(-1) (mu - mu'S^(-1)1)/(1'S^(-1)1),} 
#' see Pinheiro et al. (2014).
#' 
#' Note that the directionality (i.e. whether in "increase" in the response
#' variable is beneficial or a "decrease", is inferred from the specified
#' \samp{models} object, see [Mods()] for details).
#' 
#' Constrained contrasts (type = "constrained") add the additional constraint
#' in the optimization that the sign of the contrast coefficient for control
#' and active treatments need to be different. The quadratic programming
#' algorithm from the quadprog package is used to calculate the contrasts.
#' 
#' 
#' @aliases optContr plot.optContr plotContr
#' @param models An object of class \samp{Mods} defining the dose-response
#' shapes for which to calculate optimal contrasts.
#' @param doses Optional argument. If this argument is missing the doses
#' attribute in the \samp{Mods} object specified in \samp{models} is used.
#' @param w,S Arguments determining the matrix S used in the formula for the
#' optimal contrasts. Exactly one of \samp{w} and \samp{S} has to be specified.
#' Note that \samp{w} and \samp{S} only have to be specified up to
#' proportionality \cr \describe{ \item{w}{ Vector specifying weights for the
#' different doses, in the formula for calculation of the optimal contrasts.
#' Specifying a weights vector is equivalent to specifying S=diag(1/w) (e.g. in
#' a homoscedastic case with unequal sample sizes, \samp{w} should be
#' proportional to the group sample sizes).  } \item{S}{ Directly specify a
#' matrix proportional to the covariance matrix to use.  } }
#' @param placAdj Logical determining, whether the contrasts should be applied
#' to placebo-adjusted estimates. If yes the returned coefficients are no
#' longer contrasts (i.e. do not sum to 0). However, the result of multiplying
#' of this "contrast" matrix with the placebo adjusted estimates, will give the
#' same results as multiplying the original contrast matrix to the unadjusted
#' estimates.
#' @param type For \samp{type = "constrained"} the contrast coefficients of the
#' zero dose group are constrained to be different from the coefficients of the
#' active treatment groups. So that a weighted sum of the active treatments is
#' compared against the zero dose group. For an increasing trend the
#' coefficient of the zero dose group is negative and all other coefficients
#' have to be positive (for a decreasing trend the other way round).
#' @return Object of class \samp{optContr}. A list containing entries contMat
#' and muMat (i.e. contrast, mean and correlation matrix).
#' @author Bjoern Bornkamp
#' @seealso [MCTtest()]
#' @references Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining
#' multiple comparisons and modeling techniques in dose-response studies,
#' *Biometrics*, **61**, 738--748
#' 
#' Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014) Model-based
#' dose finding under model uncertainty using general parametric models,
#' *Statistics in Medicine*, **33**, 1646--1661
#' @examples
#' 
#' doses <- c(0,10,25,50,100,150)
#' models <- Mods(linear = NULL, emax = 25,
#'                logistic = c(50, 10.88111), exponential= 85,
#'                betaMod=rbind(c(0.33,2.31), c(1.39,1.39)),
#'                doses = doses, addArgs = list(scal = 200))
#' contMat <- optContr(models, w = rep(50,6))
#' plot(contMat)
#' plotContr(contMat) # display contrasts using ggplot2
#' 
#' ## now we would like the "contrasts" for placebo adjusted estimates
#' dosPlac <- doses[-1]
#' ## matrix proportional to cov-matrix of plac. adj. estimates for balanced data
#' S <- diag(5)+matrix(1, 5,5)
#' ## note that we explicitly hand over the doses here
#' contMat0 <- optContr(models, doses=dosPlac, S = S, placAdj = TRUE)
#' ## -> contMat0 is no longer a contrast matrix (columns do not sum to 0)
#' colSums(contMat0$contMat)
#' ## calculate contrast matrix for unadjusted estimates from this matrix
#' ## (should be same as above)
#' aux <- rbind(-colSums(contMat0$contMat), contMat0$contMat)
#' t(t(aux)/sqrt(colSums(aux^2))) ## compare to contMat$contMat
#' 
#' ## now calculate constrained contrasts 
#' if(requireNamespace("quadprog", quietly = TRUE)){
#' optContr(models, w = rep(50,6), type = "constrained")
#' optContr(models, doses=dosPlac, S = S, placAdj = TRUE,
#'          type = "constrained")
#' }
#' @export
optContr <-  function(models, doses, w, S, placAdj = FALSE,
                      type = c("unconstrained", "constrained")){
  ## calculate optimal contrasts and critical value
  if(!(inherits(models, "Mods")))
    stop("models needs to be of class Mods")
  if(missing(doses))
    doses <- attr(models, "doses")
  scal <- attr(models, "scal")
  off <- attr(models, "off")
  nodes <- attr(models, "doses")
  direction <- unique(attr(models, "direction"))
  if(length(direction) > 1)
    stop("need to provide either \"increasing\" or \"decreasing\" as direction to optContr")
  mu <- getResp(models, doses)
  if(placAdj){ 
    mu0 <- getResp(models, 0)
    mu <- mu-matrix(mu0[1,], byrow = TRUE,
                    nrow=nrow(mu), ncol=ncol(mu))
  }
  type <- match.arg(type)
  if(type == "constrained"){
    avail <- requireNamespace("quadprog", quietly = TRUE)
    if(!avail)
      stop("Need suggested package quadprog to calculate constrained contrasts")
  }
  if(any(doses == 0) & placAdj)
    stop("If placAdj == TRUE there should be no placebo group in \"doses\"")
  ## check for n and vCov arguments 
  if(!xor(missing(w), missing(S)))
    stop("Need to specify exactly one of \"w\" or \"S\"")
  if(!missing(w)){
    if(length(w) == 1){ # assume equal weights
      S <- Sinv <- diag(length(doses))
    } else {
      if(length(w) != length(doses))
        stop("w needs to be of length 1 or of the same length as doses")
      S <- diag(1/w)
      Sinv <- diag(w)
    }
  } else { 
    if(!is.matrix(S))
      stop("S needs to be a matrix")
    Sinv <- solve(S)
  }
  contMat <- modContr(mu, Sinv=Sinv, placAdj = placAdj,
                      type = type, direction = direction)
  rownames(contMat) <- doses
  corMat <- cov2cor(t(contMat) %*% S %*% contMat)
  res <- list(contMat = contMat, muMat = mu, corMat = corMat)
  attr(res, "type") <- type
  attr(res, "placAdj") <- placAdj
  class(res) <- "optContr"
  res
}

#' @export
print.optContr <- function(x, digits = 3, ...){
  cat("Optimal contrasts\n")
  print(round(x$contMat, digits))
}

#' @export
summary.optContr <- function(object, digits = 3, ...){
  class(object) <- "summary.optContr"
  print(object, digits = digits)
}

#' @export
print.summary.optContr <- function(x, digits = 3, ...){
  cat("Optimal contrasts\n")
  cat("\n","Optimal Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Contrast Correlation Matrix:","\n", sep="")
  print(round(x$corMat, digits))  
  cat("\n")
}

#' Plot optimal contrasts
#'
#' @param x,superpose,xlab,ylab,plotType Arguments for the plot method for
#' optContr objects. plotType determines, whether the contrasts or the
#' underlying (standardized) mean matrix should be plotted.
#' @param ...  Additional arguments for plot method
#' 
#' @rdname optContr
#' @method plot optContr
#' @export
plot.optContr <- function (x, superpose = TRUE, xlab = "Dose",
                           ylab = NULL, plotType = c("contrasts", "means"), ...){
  plotType <- match.arg(plotType)
  if (is.null(ylab)) {
    if (plotType == "contrasts") {
      ylab <- "Contrast coefficients"
    } else {
      ylab <- "Normalized model means"
    }
  }
  cM <- x$contMat
  if (plotType == "means")
    cM <- t(t(x$muMat)/apply(x$muMat, 2, max))
  nD <- nrow(cM)
  nM <- ncol(cM)
  cMtr <- data.frame(resp = as.vector(cM),
                     dose = rep(as.numeric(dimnames(cM)[[1]]), nM),
                     model = factor(rep(dimnames(cM)[[2]], each = nD),
                                    levels = dimnames(cM)[[2]]))
  if(superpose){
    spL <- lattice::trellis.par.get("superpose.line")
    spL$lty <- rep(spL$lty, nM%/%length(spL$lty) + 1)[1:nM]
    spL$lwd <- rep(spL$lwd, nM%/%length(spL$lwd) + 1)[1:nM]
    spL$col <- rep(spL$col, nM%/%length(spL$col) + 1)[1:nM]
    ## number of columns in legend
    nCol <- ifelse(nM < 5, nM, min(4,ceiling(nM/min(ceiling(nM/4),3))))
    key <- list(lines = spL, transparent = TRUE, 
                text = list(levels(cMtr$model), cex = 0.9),
                columns = nCol)
    ltplot <- lattice::xyplot(resp ~ dose, data = cMtr, subscripts = TRUE,
                              groups = cMtr$model, panel = panel.superpose,
                              type = "o", xlab = xlab, ylab = ylab,
                              key = key, ...)
  } else {
    ltplot <- lattice::xyplot(resp ~ dose | model, data = cMtr, type = "o", 
                              xlab = xlab, ylab = ylab,
                              strip = function(...){
                                lattice::strip.default(..., style = 1)
                              }, ...)
  }
  print(ltplot)
}

#' Plot optimal contrasts
#'
#' @param optContrObj For function \samp{plotContr} the \samp{optContrObj}
#' should contain an object of class \samp{optContr}.
#' 
#' @rdname optContr
#' @export
plotContr <- function(optContrObj, xlab = "Dose", ylab = "Contrast coefficients"){
  if(!inherits(optContrObj, "optContr"))
    stop("\"optContrObj\" needs to be of class Mods")

  parList <- attr(optContrObj$muMat, "parList")
  mod_nams <- getModNams(parList)
  cM <- optContrObj$contMat
  nD <- nrow(cM)
  nM <- ncol(cM)
  cMtr <- data.frame(resp = as.vector(cM),
                     dose = rep(as.numeric(dimnames(cM)[[1]]), nM),
                     model = factor(rep(mod_nams, each = nD), levels=mod_nams),
                     levels = dimnames(cM)[[2]])
  ggplot2::ggplot(cMtr, ggplot2::aes(.data$dose, .data$resp, col=.data$model))+
    ggplot2::geom_line(linewidth=1.2)+
    ggplot2::geom_point()+
    ggplot2::theme_bw()+
    ggplot2::geom_point(size=1.8)+
    ggplot2::xlab(xlab)+ggplot2::ylab(ylab)+
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
}
